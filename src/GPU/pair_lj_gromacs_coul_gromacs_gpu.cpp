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

#include "pair_lj_gromacs_coul_gromacs_gpu.h"

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

PairLJGromacsCoulGromacsGPU::PairLJGromacsCoulGromacsGPU(LAMMPS *lmp) :
    PairLJGromacsCoulGromacs(lmp), gpu_mode(GPU_FORCE), update_coeffs(false)
{
  respa_enable = 0;
  suffix_flag |= Suffix::GPU;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairLJGromacsCoulGromacsGPU::~PairLJGromacsCoulGromacsGPU()
{
  lj_gromacs_coul_gromacs_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairLJGromacsCoulGromacsGPU::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  if (update_coeffs) {
    lj_gromacs_coul_gromacs_gpu_reinit(atom->ntypes + 1, lj1, lj2, ljsw1, ljsw2, lj3, lj4, ljsw3,
                                       ljsw4, ljsw5);
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
    firstneigh = lj_gromacs_coul_gromacs_gpu_compute_n(
        neighbor->ago, inum, nall, atom->x, atom->type, sublo, subhi, atom->tag, atom->nspecial,
        atom->special, eflag, vflag, eflag_atom, vflag_atom, &ilist, &numneigh, success, atom->q,
        domain->boxlo, domain->prd, domain->periodicity);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    lj_gromacs_coul_gromacs_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type, ilist,
                                        numneigh, firstneigh, eflag, vflag, eflag_atom, vflag_atom,
                                        success, atom->q, atom->nlocal, domain->boxlo, domain->prd);
  }
  if (!success) error->one(FLERR, "Insufficient memory on accelerator");

  if (atom->molecular != Atom::ATOMIC && neighbor->ago == 0) neighbor->build_topology();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJGromacsCoulGromacsGPU::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR, "Pair style lj/gromacs/coul/gromacs/gpu requires atom attribute q");

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

  cut_lj_innersq = cut_lj_inner * cut_lj_inner;
  cut_ljsq = cut_lj * cut_lj;
  cut_coul_innersq = cut_coul_inner * cut_coul_inner;
  cut_coulsq = cut_coul * cut_coul;
  cut_bothsq = MAX(cut_ljsq, cut_coulsq);

  int maxspecial = 0;
  if (atom->molecular != Atom::ATOMIC) maxspecial = atom->maxspecial;
  int mnf = 5e-2 * neighbor->oneatom;
  int success = lj_gromacs_coul_gromacs_gpu_init(
      atom->ntypes + 1, lj1, lj2, ljsw1, ljsw2, lj3, lj4, ljsw3, ljsw4, ljsw5, force->special_lj,
      force->special_coul, atom->nlocal, atom->nlocal + atom->nghost, mnf, maxspecial, cell_size,
      gpu_mode, screen, force->qqrd2e, cut_bothsq, cut_ljsq, cut_lj_innersq, cut_lj_inner,
      cut_coulsq, cut_coul_innersq, cut_coul_inner, coulsw1, coulsw2, coulsw3, coulsw4, coulsw5);
  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE) neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

double PairLJGromacsCoulGromacsGPU::init_one(int i, int j)
{
  // the parameters may have been changed, e.g. by fix adapt,
  // so copy them to the GPU again before the next compute()

  update_coeffs = true;
  return PairLJGromacsCoulGromacs::init_one(i, j);
}

/* ---------------------------------------------------------------------- */

double PairLJGromacsCoulGromacsGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + lj_gromacs_coul_gromacs_gpu_bytes();
}
