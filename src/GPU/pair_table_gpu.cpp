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
   Contributing authors: Trung Dac Nguyen (ORNL)
------------------------------------------------------------------------- */

#include "pair_table_gpu.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "gpu_extra.h"
#include "lammps_gpu.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace LAMMPS_GPU;


/* ---------------------------------------------------------------------- */

PairTableGPU::PairTableGPU(LAMMPS *lmp) : PairTable(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  reinitflag = 0;
  suffix_flag |= Suffix::GPU;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairTableGPU::~PairTableGPU()
{
  table_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairTableGPU::compute(int eflag, int vflag)
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
    firstneigh =
        table_gpu_compute_n(neighbor->ago, inum, nall, atom->x, atom->type,
                            sublo, subhi, atom->tag, atom->nspecial,
                            atom->special, eflag, vflag, eflag_atom,
                            vflag_atom, &ilist, &numneigh,
                            success, domain->prd, domain->periodicity);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    table_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type, ilist, numneigh, firstneigh,
                      eflag, vflag, eflag_atom, vflag_atom, success);
  }
  if (!success) error->one(FLERR, "Insufficient memory on accelerator");

  if (atom->molecular != Atom::ATOMIC && neighbor->ago == 0)
    neighbor->build_topology();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairTableGPU::init_style()
{

  int ntypes = atom->ntypes;

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

  // pack tables and send them to device
  double ***table_coeffs = nullptr;
  double **table_data = nullptr;
  memory->create(table_coeffs, ntypes + 1, ntypes + 1, 6, "table:coeffs");

  Table *tb;
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = 1; j <= atom->ntypes; j++) {
      int n = tabindex[i][j];
      tb = &tables[n];
      table_coeffs[i][j][0] = n;
      table_coeffs[i][j][1] = tb->nshiftbits;
      table_coeffs[i][j][2] = tb->nmask;
      table_coeffs[i][j][3] = tb->innersq;
      table_coeffs[i][j][4] = tb->invdelta;
      table_coeffs[i][j][5] = tb->deltasq6;
    }

  if (tabstyle != BITMAP) {
    memory->create(table_data, ntables, 6 * tablength, "table:data");
    for (int n = 0; n < ntables; n++) {
      tb = &tables[n];
      if (tabstyle == LOOKUP) {
        for (int k = 0; k < tablength - 1; k++) {
          table_data[n][6 * k + 1] = tb->e[k];
          table_data[n][6 * k + 2] = tb->f[k];
        }
      } else if (tabstyle == LINEAR) {
        for (int k = 0; k < tablength; k++) {
          table_data[n][6 * k + 0] = tb->rsq[k];
          table_data[n][6 * k + 1] = tb->e[k];
          table_data[n][6 * k + 2] = tb->f[k];
          if (k < tablength - 1) {
            table_data[n][6 * k + 3] = tb->de[k];
            table_data[n][6 * k + 4] = tb->df[k];
          }
        }
      } else if (tabstyle == SPLINE) {
        for (int k = 0; k < tablength; k++) {
          table_data[n][6 * k + 0] = tb->rsq[k];
          table_data[n][6 * k + 1] = tb->e[k];
          table_data[n][6 * k + 2] = tb->f[k];
          table_data[n][6 * k + 3] = tb->e2[k];
          table_data[n][6 * k + 4] = tb->f2[k];
        }
      }
    }
  } else {
    int ntable = 1 << tablength;
    memory->create(table_data, ntables, 6 * ntable, "table:data");

    for (int n = 0; n < ntables; n++) {
      tb = &tables[n];
      for (int k = 0; k < ntable; k++) {
        table_data[n][6 * k + 0] = tb->rsq[k];
        table_data[n][6 * k + 1] = tb->e[k];
        table_data[n][6 * k + 2] = tb->f[k];
        table_data[n][6 * k + 3] = tb->de[k];
        table_data[n][6 * k + 4] = tb->df[k];
        table_data[n][6 * k + 5] = tb->drsq[k];
      }
    }
  }

  int maxspecial = 0;
  if (atom->molecular != Atom::ATOMIC) maxspecial = atom->maxspecial;
  int mnf = 5e-2 * neighbor->oneatom;
  int success = table_gpu_init(atom->ntypes + 1, cutsq, table_coeffs, table_data, force->special_lj,
                               atom->nlocal, atom->nlocal + atom->nghost, mnf, maxspecial,
                               cell_size, gpu_mode, screen, tabstyle, ntables, tablength);
  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE) neighbor->add_request(this, NeighConst::REQ_FULL);
  memory->destroy(table_coeffs);
  memory->destroy(table_data);
}

/* ---------------------------------------------------------------------- */

double PairTableGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + table_gpu_bytes();
}
