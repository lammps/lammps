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
   Contributing author: Claude Opus 5.5 (Anthropic), under the direction of
   Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_morse_soft_gpu.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "gpu_extra.h"
#include "lammps_gpu.h"
#include "math_special.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace LAMMPS_GPU;

/* ---------------------------------------------------------------------- */

PairMorseSoftGPU::PairMorseSoftGPU(LAMMPS *lmp) :
    PairMorseSoft(lmp), gpu_mode(GPU_FORCE), update_coeffs(false), bs(nullptr), vs(nullptr)
{
  suffix_flag |= Suffix::GPU;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairMorseSoftGPU::~PairMorseSoftGPU()
{
  morse_soft_gpu_clear();
  memory->destroy(bs);
  memory->destroy(vs);
}

/* ---------------------------------------------------------------------- */

void PairMorseSoftGPU::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  if (update_coeffs) {
    setup_soft_factors();
    morse_soft_gpu_reinit(atom->ntypes + 1, cutsq, morse1, r0, alpha, d0, offset, bs, vs);
    update_coeffs = false;
  }

  int nall = atom->nlocal + atom->nghost;
  int inum;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
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
    firstneigh = morse_soft_gpu_compute_n(neighbor->ago, inum, nall, atom->x, atom->type, sublo,
                                          subhi, atom->tag, atom->nspecial, atom->special, eflag,
                                          vflag, eflag_atom, vflag_atom, &ilist, &numneigh, success,
                                          domain->prd, domain->periodicity);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    morse_soft_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type, ilist, numneigh,
                           firstneigh, eflag, vflag, eflag_atom, vflag_atom, success);
  }
  if (!success) error->one(FLERR, "Insufficient memory on accelerator");

  if (atom->molecular != Atom::ATOMIC && neighbor->ago == 0) neighbor->build_topology();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairMorseSoftGPU::init_style()
{

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
  setup_soft_factors();

  int maxspecial = 0;
  if (atom->molecular != Atom::ATOMIC) maxspecial = atom->maxspecial;
  int mnf = 5e-2 * neighbor->oneatom;
  int success = morse_soft_gpu_init(atom->ntypes + 1, cutsq, morse1, r0, alpha, d0, offset, bs, vs,
                                    force->special_lj, atom->nlocal, atom->nlocal + atom->nghost,
                                    mnf, maxspecial, cell_size, gpu_mode, screen);
  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE) neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

double PairMorseSoftGPU::init_one(int i, int j)
{
  // the parameters may have been changed, e.g. by fix adapt,
  // so copy them to the GPU again before the next compute()

  update_coeffs = true;
  return PairMorseSoft::init_one(i, j);
}

/* ----------------------------------------------------------------------
   combine the lambda dependent factors of PairMorseSoft::compute() into
   the prefactor of the exponential^3 term (bs) and the scale factor of
   the Morse potential (vs) for each pair of atom types
------------------------------------------------------------------------- */

void PairMorseSoftGPU::setup_soft_factors()
{
  const int np1 = atom->ntypes + 1;
  memory->destroy(bs);
  memory->destroy(vs);
  memory->create(bs, np1, np1, "pair:bs");
  memory->create(vs, np1, np1, "pair:vs");

  for (int i = 0; i < np1; i++) {
    for (int j = 0; j < np1; j++) {
      bs[i][j] = vs[i][j] = 0.0;
      if ((i == 0) || (j == 0) || (cutsq[i][j] == 0.0)) continue;
      const double ea = exp(alpha[i][j] * r0[i][j]);
      const double iea2 = exp(-2. * alpha[i][j] * r0[i][j]);
      const double B = -2.0 * d0[i][j] * iea2 * (ea - 1.0) / 3.0;
      const double l = lambda[i][j];
      if (l >= shift_range) {
        bs[i][j] = B * (l - 1.0) / (shift_range - 1.0);
        vs[i][j] = 1.0;
      } else {
        const double llf = MathSpecial::powint(l / shift_range, nlambda);
        bs[i][j] = B * llf;
        vs[i][j] = llf;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairMorseSoftGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + morse_soft_gpu_bytes();
}
