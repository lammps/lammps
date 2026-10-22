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
   Contributing author: Trung Nguyen (Northwestern)
------------------------------------------------------------------------- */

#include "pair_coul_long_cs_gpu.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
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

/* ---------------------------------------------------------------------- */

PairCoulLongCSGPU::PairCoulLongCSGPU(LAMMPS *lmp) : PairCoulLongCS(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  suffix_flag |= Suffix::GPU;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairCoulLongCSGPU::~PairCoulLongCSGPU()
{
  clcs_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairCoulLongCSGPU::compute(int eflag, int vflag)
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
    firstneigh = clcs_gpu_compute_n(neighbor->ago, inum, nall, atom->x, atom->type, sublo, subhi,
                                    atom->tag, atom->nspecial, atom->special, eflag, vflag,
                                    eflag_atom, vflag_atom, &ilist, &numneigh, success, atom->q,
                                    domain->boxlo, domain->prd, domain->periodicity);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    clcs_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type, ilist, numneigh, firstneigh,
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

void PairCoulLongCSGPU::init_style()
{
  cut_respa = nullptr;

  if (!atom->q_flag) error->all(FLERR, "Pair style coul/long/cs/gpu requires atom attribute q");

  // Call init_one calculation make sure scale is correct
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) { init_one(i, j); }
    }
  }
  double cell_size = cut_coul + neighbor->skin;

  cut_coulsq = cut_coul * cut_coul;

  // ensure use of KSpace long-range solver, set g_ewald

  if (force->kspace == nullptr) error->all(FLERR, "Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;

  // setup force tables

  if (ncoultablebits) init_tables(cut_coul, cut_respa);

  int maxspecial = 0;
  if (atom->molecular != Atom::ATOMIC) maxspecial = atom->maxspecial;
  int mnf = 5e-2 * neighbor->oneatom;
  int success = clcs_gpu_init(atom->ntypes + 1, scale, atom->nlocal, atom->nlocal + atom->nghost,
                              mnf, maxspecial, cell_size, gpu_mode, screen, cut_coulsq,
                              force->special_coul, force->qqrd2e, g_ewald);

  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE) neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

void PairCoulLongCSGPU::reinit()
{
  Pair::reinit();

  clcs_gpu_reinit(atom->ntypes + 1, scale);
}

/* ---------------------------------------------------------------------- */

double PairCoulLongCSGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + clcs_gpu_bytes();
}
