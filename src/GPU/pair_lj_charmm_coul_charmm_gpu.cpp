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

#include "pair_lj_charmm_coul_charmm_gpu.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "gpu_extra.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>
#include "lammps_gpu.h"

using namespace LAMMPS_NS;
using namespace LAMMPS_GPU;


/* ---------------------------------------------------------------------- */

PairLJCharmmCoulCharmmGPU::PairLJCharmmCoulCharmmGPU(LAMMPS *lmp) :
    PairLJCharmmCoulCharmm(lmp), gpu_mode(GPU_FORCE)
{
  reinitflag = 0;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairLJCharmmCoulCharmmGPU::~PairLJCharmmCoulCharmmGPU()
{
  crm_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulCharmmGPU::compute(int eflag, int vflag)
{
  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0;

  int nall = atom->nlocal + atom->nghost;
  int inum;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    inum = atom->nlocal;
    firstneigh = crm_gpu_compute_n(neighbor->ago, inum, nall, atom->x, atom->type, domain->sublo,
                                   domain->subhi, atom->tag, atom->nspecial, atom->special, eflag,
                                   vflag, eflag_atom, vflag_atom, &ilist, &numneigh,
                                   success, atom->q, domain->boxlo, domain->prd,
                                   domain->periodicity);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    crm_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type, ilist, numneigh, firstneigh,
                    eflag, vflag, eflag_atom, vflag_atom, success, atom->q,
                    atom->nlocal, domain->boxlo, domain->prd);
  }
  if (!success) error->one(FLERR, "Insufficient memory on accelerator");

  if (atom->molecular != Atom::ATOMIC && neighbor->ago == 0)
    neighbor->build_topology();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCharmmCoulCharmmGPU::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR, "Pair style lj/charmm/coul/long/gpu requires atom attribute q");

  // Repeated cutsq calculation in init_one() is required for GPU package

  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) init_one(i, j);
    }
  }

  cut_lj_innersq = cut_lj_inner * cut_lj_inner;
  cut_coul_innersq = cut_coul_inner * cut_coul_inner;
  cut_ljsq = cut_lj * cut_lj;
  cut_coulsq = cut_coul * cut_coul;
  cut_bothsq = MAX(cut_ljsq, cut_coulsq);

  denom_lj =
      (cut_ljsq - cut_lj_innersq) * (cut_ljsq - cut_lj_innersq) * (cut_ljsq - cut_lj_innersq);
  denom_lj = 1.0 / denom_lj;

  denom_coul = (cut_coulsq - cut_coul_innersq) * (cut_coulsq - cut_coul_innersq) *
      (cut_coulsq - cut_coul_innersq);
  denom_coul = 1.0 / denom_coul;

  double cell_size = sqrt(cut_bothsq) + neighbor->skin;

  int maxspecial = 0;
  if (atom->molecular != Atom::ATOMIC) maxspecial = atom->maxspecial;

  bool arithmetic = true;
  for (int i = 1; i < atom->ntypes + 1; i++)
    for (int j = i + 1; j < atom->ntypes + 1; j++) {
      if (epsilon[i][j] != sqrt(epsilon[i][i] * epsilon[j][j])) arithmetic = false;
      if (sigma[i][j] != 0.5 * (sigma[i][i] + sigma[j][j])) arithmetic = false;
    }

  int mnf = 5e-2 * neighbor->oneatom;
  int success =
      crm_gpu_init(atom->ntypes + 1, cut_bothsq, lj1, lj2, lj3, lj4, force->special_lj,
                   atom->nlocal, atom->nlocal + atom->nghost, mnf, maxspecial, cell_size, gpu_mode,
                   screen, cut_ljsq, cut_coulsq, force->special_coul, force->qqrd2e, cut_lj_innersq,
                   cut_coul_innersq, denom_lj, denom_coul, epsilon, sigma, arithmetic);
  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE) neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

double PairLJCharmmCoulCharmmGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + crm_gpu_bytes();
}

