/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "pair_resquared_gpu.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "gpu_extra.h"
#include "lammps_gpu.h"
#include "math_extra.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace LAMMPS_GPU;


/* ---------------------------------------------------------------------- */

PairRESquaredGPU::PairRESquaredGPU(LAMMPS *lmp) : PairRESquared(lmp), gpu_mode(GPU_FORCE)
{
  reinitflag = 0;
  avec = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  if (!avec) error->all(FLERR, "Pair resquared/gpu requires atom style ellipsoid");
  suffix_flag |= Suffix::GPU;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairRESquaredGPU::~PairRESquaredGPU()
{
  re_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairRESquaredGPU::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  int nall = atom->nlocal + atom->nghost;
  int inum;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;

  if (gpu_mode != GPU_FORCE) {
    double sublo[3], subhi[3];
    if (domain->triclinic == 0) {
      sublo[0] = domain->sublo[0];
      sublo[1] = domain->sublo[1];
      sublo[2] = domain->sublo[2];
      subhi[0] = domain->subhi[0];
      subhi[1] = domain->subhi[1];
      subhi[2] = domain->subhi[2];
    } else {
      domain->bbox(domain->sublo_lamda, domain->subhi_lamda, sublo, subhi);
    }
    inum = atom->nlocal;
    firstneigh =
        re_gpu_compute_n(neighbor->ago, inum, nall, atom->x, atom->type, sublo, subhi, atom->tag,
                         atom->nspecial, atom->special, eflag, vflag, eflag_atom, vflag_atom,
                         &ilist, &numneigh, success, ellipsoid, bonus);
  } else {
    inum = list->inum;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    ilist = re_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type,
                           list->ilist, numneigh, firstneigh, eflag, vflag,
                           eflag_atom, vflag_atom, success, ellipsoid, bonus);
  }
  if (!success) error->one(FLERR, "Insufficient memory on accelerator");

  if (atom->molecular != Atom::ATOMIC && neighbor->ago == 0)
    neighbor->build_topology();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairRESquaredGPU::init_style()
{
  if (!atom->ellipsoid_flag) error->all(FLERR, "Pair resquared/gpu requires atom style ellipsoid");

  // per-type shape precalculations
  // require that atom shapes are identical within each type
  // if shape = 0 for point particle, set shape = 1 as required by Gay-Berne

  for (int i = 1; i <= atom->ntypes; i++) {
    if (!atom->shape_consistency(i, shape1[i][0], shape1[i][1], shape1[i][2]))
      error->all(FLERR, "Pair resquared/gpu requires atoms with same type have same shape");
    if (setwell[i]) {
      shape2[i][0] = shape1[i][0] * shape1[i][0];
      shape2[i][1] = shape1[i][1] * shape1[i][1];
      shape2[i][2] = shape1[i][2] * shape1[i][2];
      lshape[i] = shape1[i][0] * shape1[i][1] * shape1[i][2];
    }
  }

  // Repeat cutsq calculation because done after call to init_style
  double maxcut = -1.0;
  double cut;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) {
        cut = init_one(i, j);
        cut *= cut;
        if (cut > maxcut) maxcut = cut;
        cutsq[i][j] = cutsq[j][i] = cut;
      } else
        cutsq[i][j] = cutsq[j][i] = 0.0;
    }
  }

  double cell_size = sqrt(maxcut) + neighbor->skin;

  int maxspecial = 0;
  if (atom->molecular != Atom::ATOMIC) maxspecial = atom->maxspecial;
  int mnf = 5e-2 * neighbor->oneatom;
  int success =
      re_gpu_init(atom->ntypes + 1, shape1, well, cutsq, sigma, epsilon, form,
                  lj1, lj2, lj3, lj4, offset, force->special_lj, atom->nlocal,
                  atom->nlocal + atom->nghost, mnf, maxspecial, cell_size,
                  gpu_mode, screen);
  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE) neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

double PairRESquaredGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + re_gpu_bytes();
}
