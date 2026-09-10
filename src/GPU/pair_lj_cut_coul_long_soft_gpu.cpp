/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Trung Nguyen (U Chicago)
------------------------------------------------------------------------- */

#include "pair_lj_cut_coul_long_soft_gpu.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "gpu_extra.h"
#include "kspace.h"
#include "lammps_gpu.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace LAMMPS_GPU;
using namespace EwaldConst;


/* ---------------------------------------------------------------------- */

PairLJCutCoulLongSoftGPU::PairLJCutCoulLongSoftGPU(LAMMPS *lmp) :
    PairLJCutCoulLongSoft(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  suffix_flag |= Suffix::GPU;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairLJCutCoulLongSoftGPU::~PairLJCutCoulLongSoftGPU()
{
  ljcls_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulLongSoftGPU::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

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
    firstneigh = ljcls_gpu_compute_n(neighbor->ago, inum, nall, atom->x, atom->type, sublo, subhi,
                                     atom->tag, atom->nspecial, atom->special, eflag, vflag,
                                     eflag_atom, vflag_atom, &ilist, &numneigh, success, atom->q,
                                     domain->boxlo, domain->prd, domain->periodicity);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    ljcls_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type, ilist, numneigh, firstneigh,
                     eflag, vflag, eflag_atom, vflag_atom, success, atom->q,
                     atom->nlocal, domain->boxlo, domain->prd);
  }
  if (!success) error->one(FLERR, "Insufficient memory on accelerator");

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCutCoulLongSoftGPU::init_style()
{
  cut_respa = nullptr;

  if (!atom->q_flag) error->all(FLERR, "Pair style lj/cut/coul/long/soft/gpu requires atom attribute q");

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

  cut_coulsq = cut_coul * cut_coul;

  // insure use of KSpace long-range solver, set g_ewald

  if (force->kspace == nullptr) error->all(FLERR, "Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;

  // setup force tables

  if (ncoultablebits) init_tables(cut_coul, cut_respa);

  int maxspecial = 0;
  if (atom->molecular != Atom::ATOMIC) maxspecial = atom->maxspecial;
  int mnf = 5e-2 * neighbor->oneatom;
  int success =
      ljcls_gpu_init(atom->ntypes + 1, cutsq, lj1, lj2, lj3, lj4, offset, epsilon, force->special_lj,
                    atom->nlocal, atom->nlocal + atom->nghost, mnf, maxspecial, cell_size, gpu_mode,
                    screen, cut_ljsq, cut_coulsq, force->special_coul, force->qqrd2e, g_ewald);
  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE) neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulLongSoftGPU::reinit()
{
  Pair::reinit();

  ljcls_gpu_reinit(atom->ntypes + 1, cutsq, lj1, lj2, lj3, lj4, offset, epsilon, cut_ljsq);
}

/* ---------------------------------------------------------------------- */

double PairLJCutCoulLongSoftGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + ljcls_gpu_bytes();
}
