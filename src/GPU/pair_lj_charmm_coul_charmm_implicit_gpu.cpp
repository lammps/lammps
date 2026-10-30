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

#include "pair_lj_charmm_coul_charmm_implicit_gpu.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "gpu_extra.h"
#include "lammps_gpu.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace LAMMPS_GPU;

/* ---------------------------------------------------------------------- */

PairLJCharmmCoulCharmmImplicitGPU::PairLJCharmmCoulCharmmImplicitGPU(LAMMPS *lmp) :
    PairLJCharmmCoulCharmmImplicit(lmp), gpu_mode(GPU_FORCE), update_coeffs(false)
{
  respa_enable = 0;
  suffix_flag |= Suffix::GPU;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairLJCharmmCoulCharmmImplicitGPU::~PairLJCharmmCoulCharmmImplicitGPU()
{
  lj_charmm_coul_charmm_implicit_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulCharmmImplicitGPU::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  if (update_coeffs) {
    lj_charmm_coul_charmm_implicit_gpu_reinit(atom->ntypes + 1, lj1, lj2, lj3, lj4);
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
    firstneigh = lj_charmm_coul_charmm_implicit_gpu_compute_n(
        neighbor->ago, inum, nall, atom->x, atom->type, sublo, subhi, atom->tag, atom->nspecial,
        atom->special, eflag, vflag, eflag_atom, vflag_atom, &ilist, &numneigh, success, atom->q,
        domain->boxlo, domain->prd, domain->periodicity);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    lj_charmm_coul_charmm_implicit_gpu_compute(
        neighbor->ago, inum, nall, atom->x, atom->type, ilist, numneigh, firstneigh, eflag, vflag,
        eflag_atom, vflag_atom, success, atom->q, atom->nlocal, domain->boxlo, domain->prd);
  }
  if (!success) error->one(FLERR, "Insufficient memory on accelerator");

  if (atom->molecular != Atom::ATOMIC && neighbor->ago == 0) neighbor->build_topology();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCharmmCoulCharmmImplicitGPU::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR, "Pair style lj/charmm/coul/charmm/implicit/gpu requires atom attribute q");

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

  // require cut_lj_inner < cut_lj, cut_coul_inner < cut_coul

  if (cut_lj_inner >= cut_lj || cut_coul_inner >= cut_coul)
    error->all(FLERR, "Pair inner cutoff >= Pair outer cutoff");

  cut_lj_innersq = cut_lj_inner * cut_lj_inner;
  cut_ljsq = cut_lj * cut_lj;
  cut_coul_innersq = cut_coul_inner * cut_coul_inner;
  cut_coulsq = cut_coul * cut_coul;
  cut_bothsq = MAX(cut_ljsq, cut_coulsq);

  denom_lj =
      (cut_ljsq - cut_lj_innersq) * (cut_ljsq - cut_lj_innersq) * (cut_ljsq - cut_lj_innersq);
  denom_coul = (cut_coulsq - cut_coul_innersq) * (cut_coulsq - cut_coul_innersq) *
      (cut_coulsq - cut_coul_innersq);

  int maxspecial = 0;
  if (atom->molecular != Atom::ATOMIC) maxspecial = atom->maxspecial;
  int mnf = 5e-2 * neighbor->oneatom;
  int success = lj_charmm_coul_charmm_implicit_gpu_init(
      atom->ntypes + 1, lj1, lj2, lj3, lj4, force->special_lj, force->special_coul, atom->nlocal,
      atom->nlocal + atom->nghost, mnf, maxspecial, cell_size, gpu_mode, screen, force->qqrd2e,
      cut_bothsq, cut_ljsq, cut_lj_innersq, cut_coulsq, cut_coul_innersq, 1.0 / denom_lj,
      1.0 / denom_coul);
  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE) neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

double PairLJCharmmCoulCharmmImplicitGPU::init_one(int i, int j)
{
  // the parameters may have been changed, e.g. by fix adapt,
  // so copy them to the GPU again before the next compute()

  update_coeffs = true;
  return PairLJCharmmCoulCharmmImplicit::init_one(i, j);
}

/* ---------------------------------------------------------------------- */

double PairLJCharmmCoulCharmmImplicitGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + lj_charmm_coul_charmm_implicit_gpu_bytes();
}
