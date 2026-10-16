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

#include "pair_lj_charmmfsw_coul_long_gpu.h"

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

PairLJCharmmfswCoulLongGPU::PairLJCharmmfswCoulLongGPU(LAMMPS *lmp) :
    PairLJCharmmfswCoulLong(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  reinitflag = 0;
  suffix_flag |= Suffix::GPU;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairLJCharmmfswCoulLongGPU::~PairLJCharmmfswCoulLongGPU()
{
  charmmfsw_long_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmfswCoulLongGPU::compute(int eflag, int vflag)
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
    firstneigh = charmmfsw_long_gpu_compute_n(
        neighbor->ago, inum, nall, atom->x, atom->type, sublo, subhi, atom->tag, atom->nspecial,
        atom->special, eflag, vflag, eflag_atom, vflag_atom, &ilist, &numneigh, success, atom->q,
        domain->boxlo, domain->prd, domain->periodicity);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    charmmfsw_long_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type, ilist, numneigh,
                               firstneigh, eflag, vflag, eflag_atom, vflag_atom, success, atom->q,
                               atom->nlocal, domain->boxlo, domain->prd);
  }
  if (!success) error->one(FLERR, "Insufficient memory on accelerator");

  if (atom->molecular != Atom::ATOMIC && neighbor->ago == 0) neighbor->build_topology();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCharmmfswCoulLongGPU::init_style()
{
  cut_respa = nullptr;

  if (!atom->q_flag)
    error->all(FLERR, "Pair style lj/charmmfsw/coul/long/gpu requires atom attribute q");

  // Repeat cutsq calculation because done after call to init_style

  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) init_one(i, j);
    }
  }

  // require cut_lj_inner < cut_lj

  if (cut_lj_inner >= cut_lj) error->all(FLERR, "Pair inner cutoff >= Pair outer cutoff");

  // same derived cutoff values as in PairLJCharmmfswCoulLong::init_style()

  cut_lj_innersq = cut_lj_inner * cut_lj_inner;
  cut_ljsq = cut_lj * cut_lj;
  cut_ljinv = 1.0 / cut_lj;
  cut_lj_innerinv = 1.0 / cut_lj_inner;
  cut_lj3 = cut_lj * cut_lj * cut_lj;
  cut_lj3inv = cut_ljinv * cut_ljinv * cut_ljinv;
  cut_lj_inner3inv = cut_lj_innerinv * cut_lj_innerinv * cut_lj_innerinv;
  cut_lj_inner3 = cut_lj_inner * cut_lj_inner * cut_lj_inner;
  cut_lj6 = cut_ljsq * cut_ljsq * cut_ljsq;
  cut_lj6inv = cut_lj3inv * cut_lj3inv;
  cut_lj_inner6inv = cut_lj_inner3inv * cut_lj_inner3inv;
  cut_lj_inner6 = cut_lj_innersq * cut_lj_innersq * cut_lj_innersq;
  cut_coulsq = cut_coul * cut_coul;
  cut_bothsq = MAX(cut_ljsq, cut_coulsq);

  denom_lj =
      (cut_ljsq - cut_lj_innersq) * (cut_ljsq - cut_lj_innersq) * (cut_ljsq - cut_lj_innersq);
  denom_lj12 = 1.0 / (cut_lj6 - cut_lj_inner6);
  denom_lj6 = 1.0 / (cut_lj3 - cut_lj_inner3);

  // constants of the force-switched LJ energy in the order the GPU kernel expects

  double fsw[8] = {cut_lj6,   denom_lj12, cut_lj6inv,       cut_lj3,
                   denom_lj6, cut_lj3inv, cut_lj_inner6inv, cut_lj_inner3inv};

  double cell_size = sqrt(cut_bothsq) + neighbor->skin;

  // ensure use of KSpace long-range solver, set g_ewald

  if (force->kspace == nullptr) error->all(FLERR, "Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;

  // setup force tables

  if (ncoultablebits) init_tables(cut_coul, cut_respa);

  int maxspecial = 0;
  if (atom->molecular != Atom::ATOMIC) maxspecial = atom->maxspecial;

  bool arithmetic = true;
  for (int i = 1; i < atom->ntypes + 1; i++)
    for (int j = i + 1; j < atom->ntypes + 1; j++) {
      if (epsilon[i][j] != sqrt(epsilon[i][i] * epsilon[j][j])) arithmetic = false;
      if (sigma[i][j] != 0.5 * (sigma[i][i] + sigma[j][j])) arithmetic = false;
    }

  int mnf = 5e-2 * neighbor->oneatom;
  int success = charmmfsw_long_gpu_init(
      atom->ntypes + 1, cut_bothsq, lj1, lj2, lj3, lj4, offset, force->special_lj, atom->nlocal,
      atom->nlocal + atom->nghost, mnf, maxspecial, cell_size, gpu_mode, screen, cut_ljsq,
      cut_coulsq, force->special_coul, force->qqrd2e, g_ewald, cut_lj_innersq, denom_lj, epsilon,
      sigma, arithmetic, fsw);
  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE) neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

double PairLJCharmmfswCoulLongGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + charmmfsw_long_gpu_bytes();
}
