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
   Contributing author: Eddy BARRAUD (IFPEN/Sorbonne)
------------------------------------------------------------------------- */

#include "pair_dpd_coul_slater_long_gpu.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "gpu_extra.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"
#include "update.h"

#include "ewald_const.h"
#include "kspace.h"
#include "lammps_gpu.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace LAMMPS_GPU;
using namespace EwaldConst;


/* ---------------------------------------------------------------------- */

PairDPDCoulSlaterLongGPU::PairDPDCoulSlaterLongGPU(LAMMPS *lmp) : PairDPDCoulSlaterLong(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  reinitflag = 0;
  suffix_flag |= Suffix::GPU;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairDPDCoulSlaterLongGPU::~PairDPDCoulSlaterLongGPU()
{
  dpd_coul_slater_long_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairDPDCoulSlaterLongGPU::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  int nall = atom->nlocal + atom->nghost;
  int inum;

  double dtinvsqrt = 1.0 / sqrt(update->dt);

  bool success = true;
  int *ilist, *numneigh, **firstneigh;

  double *q = atom->q;
  dpd_coul_slater_long_gpu_get_extra_data(q);

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
    firstneigh = dpd_coul_slater_long_gpu_compute_n(
        neighbor->ago, inum, nall, atom->x, atom->type, sublo, subhi, atom->tag, atom->nspecial,
        atom->special, eflag, vflag, eflag_atom, vflag_atom, &ilist, &numneigh,
        success, atom->v, dtinvsqrt, seed, update->ntimestep, domain->boxlo, domain->prd);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    dpd_coul_slater_long_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type, ilist, numneigh, firstneigh,
                    eflag, vflag, eflag_atom, vflag_atom, success, atom->tag,
                    atom->v, dtinvsqrt, seed, update->ntimestep, atom->nlocal, domain->boxlo,
                    domain->prd);
  }
  if (!success) error->one(FLERR, "Insufficient memory on accelerator");

  if (atom->molecular != Atom::ATOMIC && neighbor->ago == 0)
    neighbor->build_topology();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairDPDCoulSlaterLongGPU::init_style()
{

  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair dpd requires ghost atoms store velocity");
  if (!atom->q_flag)
    error->all(FLERR,"Pair style coul/slater/long requires atom attribute q");
  // ensure use of KSpace long-range solver, set g_ewald
  if (force->kspace == nullptr) error->all(FLERR, "Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;

  // Repeat cutsq calculation because done after call to init_style
  double maxcut = -1.0;
  double mcut;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) {
        mcut = init_one(i, j);
        mcut *= mcut;
        if (mcut > maxcut) maxcut = mcut;
        cutsq[i][j] = cutsq[j][i] = mcut;
      } else
        cutsq[i][j] = cutsq[j][i] = 0.0;
    }
  }
  double cell_size = sqrt(maxcut) + neighbor->skin;

  int maxspecial = 0;
  if (atom->molecular != Atom::ATOMIC) maxspecial = atom->maxspecial;
  int mnf = 5e-2 * neighbor->oneatom;
  int success =
      dpd_coul_slater_long_gpu_init(atom->ntypes + 1, cutsq, a0, gamma, sigma,
                   cut_dpd, cut_dpdsq, cut_slatersq, force->special_lj, atom->nlocal,
                   atom->nlocal + atom->nghost, mnf, maxspecial, cell_size, gpu_mode, screen,
                   force->special_coul, force->qqrd2e, g_ewald, lamda);
  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE) neighbor->add_request(this, NeighConst::REQ_FULL);

}

/* ---------------------------------------------------------------------- */

double PairDPDCoulSlaterLongGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + dpd_coul_slater_long_gpu_bytes();
}
