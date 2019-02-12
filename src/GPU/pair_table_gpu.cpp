/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Trung Dac Nguyen (ORNL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_table_gpu.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "integrate.h"
#include "memory.h"
#include "error.h"
#include "neigh_request.h"
#include "universe.h"
#include "update.h"
#include "domain.h"
#include "gpu_extra.h"

#define LOOKUP 0
#define LINEAR 1
#define SPLINE 2
#define BITMAP 3

using namespace LAMMPS_NS;

// External functions from cuda library for atom decomposition

int table_gpu_init(const int ntypes, double **cutsq,
                   double ***host_table_coeffs, double **host_table_data,
                   double *special_lj, const int nlocal, const int nall,
                   const int max_nbors, const int maxspecial,
                   const double cell_size, int &gpu_mode, FILE *screen,
                   int tabstyle, int ntables, int tablength);
void table_gpu_clear();
int ** table_gpu_compute_n(const int ago, const int inum, const int nall,
                           double **host_x, int *host_type, double *sublo,
                           double *subhi, tagint *tag, int **nspecial,
                           tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int &host_start,
                           int **ilist, int **jnum, const double cpu_time,
                           bool &success);
void table_gpu_compute(const int ago, const int inum, const int nall,
                       double **host_x, int *host_type, int *ilist, int *numj,
                       int **firstneigh, const bool eflag, const bool vflag,
                       const bool eatom, const bool vatom, int &host_start,
                       const double cpu_time, bool &success);
double table_gpu_bytes();

/* ---------------------------------------------------------------------- */

PairTableGPU::PairTableGPU(LAMMPS *lmp) : PairTable(lmp),
                                          gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  reinitflag = 0;
  cpu_time = 0.0;
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
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    inum = atom->nlocal;
    firstneigh = table_gpu_compute_n(neighbor->ago, inum, nall, atom->x,
                                     atom->type, domain->sublo, domain->subhi,
                                     atom->tag, atom->nspecial, atom->special,
                                     eflag, vflag, eflag_atom, vflag_atom,
                                     host_start, &ilist, &numneigh, cpu_time,
                                     success);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    table_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type,
                      ilist, numneigh, firstneigh, eflag, vflag, eflag_atom,
                      vflag_atom, host_start, cpu_time, success);
  }
  if (!success)
    error->one(FLERR,"Insufficient memory on accelerator");

  if (host_start<inum) {
    cpu_time = MPI_Wtime();
    cpu_compute(host_start, inum, eflag, vflag, ilist, numneigh, firstneigh);
    cpu_time = MPI_Wtime() - cpu_time;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairTableGPU::init_style()
{
  if (force->newton_pair)
    error->all(FLERR,"Cannot use newton pair with table/gpu pair style");

  int ntypes = atom->ntypes;

  // Repeat cutsq calculation because done after call to init_style
  double maxcut = -1.0;
  double cut;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) {
        cut = init_one(i,j);
        cut *= cut;
        if (cut > maxcut)
          maxcut = cut;
        cutsq[i][j] = cutsq[j][i] = cut;
      } else
        cutsq[i][j] = cutsq[j][i] = 0.0;
    }
  }
  double cell_size = sqrt(maxcut) + neighbor->skin;

  // pack tables and send them to device
  double ***table_coeffs = NULL;
  double **table_data = NULL;
  memory->create(table_coeffs, ntypes+1, ntypes+1, 6, "table:coeffs");

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
    memory->create(table_data, ntables, 6*tablength, "table:data");
    for (int n = 0; n < ntables; n++) {
      tb = &tables[n];
      if (tabstyle == LOOKUP) {
        for (int k = 0; k<tablength-1; k++) {
          table_data[n][6*k+1] = tb->e[k];
          table_data[n][6*k+2] = tb->f[k];
        }
      } else if (tabstyle == LINEAR) {
        for (int k = 0; k<tablength; k++) {
          table_data[n][6*k+0] = tb->rsq[k];
          table_data[n][6*k+1] = tb->e[k];
          table_data[n][6*k+2] = tb->f[k];
          if (k<tablength-1) {
            table_data[n][6*k+3] = tb->de[k];
            table_data[n][6*k+4] = tb->df[k];
          }
       }
      } else if (tabstyle == SPLINE) {
        for (int k = 0; k<tablength; k++) {
          table_data[n][6*k+0] = tb->rsq[k];
          table_data[n][6*k+1] = tb->e[k];
          table_data[n][6*k+2] = tb->f[k];
          table_data[n][6*k+3] = tb->e2[k];
          table_data[n][6*k+4] = tb->f2[k];
        }
      }
    }
  } else {
    int ntable = 1 << tablength;
    memory->create(table_data, ntables, 6*ntable, "table:data");

    for (int n = 0; n < ntables; n++) {
      tb = &tables[n];
      for (int k = 0; k<ntable; k++) {
        table_data[n][6*k+0] = tb->rsq[k];
        table_data[n][6*k+1] = tb->e[k];
        table_data[n][6*k+2] = tb->f[k];
        table_data[n][6*k+3] = tb->de[k];
        table_data[n][6*k+4] = tb->df[k];
        table_data[n][6*k+5] = tb->drsq[k];
      }
    }
  }

  int maxspecial=0;
  if (atom->molecular)
    maxspecial=atom->maxspecial;
  int success = table_gpu_init(atom->ntypes+1, cutsq, table_coeffs, table_data,
                               force->special_lj, atom->nlocal,
                               atom->nlocal+atom->nghost, 300, maxspecial,
                               cell_size, gpu_mode, screen, tabstyle, ntables,
                               tablength);
  GPU_EXTRA::check_flag(success,error,world);

  if (gpu_mode == GPU_FORCE) {
    int irequest = neighbor->request(this,instance_me);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
  }

  memory->destroy(table_coeffs);
  memory->destroy(table_data);
}

/* ---------------------------------------------------------------------- */

double PairTableGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + table_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairTableGPU::cpu_compute(int start, int inum, int eflag, int /* vflag */,
                               int *ilist, int *numneigh, int **firstneigh) {
  int i,j,ii,jj,jnum,itype,jtype,itable;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,factor_lj,fraction,value,a,b;
  int *jlist;
  Table *tb;

  union_int_float_t rsq_lookup;
  int tlm1 = tablength - 1;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  double *special_lj = force->special_lj;

  // loop over neighbors of my atoms

  for (ii = start; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        tb = &tables[tabindex[itype][jtype]];
        if (rsq < tb->innersq)
          error->one(FLERR,"Pair distance < table inner cutoff");

        if (tabstyle == LOOKUP) {
          itable = static_cast<int> ((rsq - tb->innersq) * tb->invdelta);
          if (itable >= tlm1)
            error->one(FLERR,"Pair distance > table outer cutoff");
          fpair = factor_lj * tb->f[itable];
        } else if (tabstyle == LINEAR) {
          itable = static_cast<int> ((rsq - tb->innersq) * tb->invdelta);
          if (itable >= tlm1)
            error->one(FLERR,"Pair distance > table outer cutoff");
          fraction = (rsq - tb->rsq[itable]) * tb->invdelta;
          value = tb->f[itable] + fraction*tb->df[itable];
          fpair = factor_lj * value;
        } else if (tabstyle == SPLINE) {
          itable = static_cast<int> ((rsq - tb->innersq) * tb->invdelta);
          if (itable >= tlm1)
            error->one(FLERR,"Pair distance > table outer cutoff");
          b = (rsq - tb->rsq[itable]) * tb->invdelta;
          a = 1.0 - b;
          value = a * tb->f[itable] + b * tb->f[itable+1] +
            ((a*a*a-a)*tb->f2[itable] + (b*b*b-b)*tb->f2[itable+1]) *
            tb->deltasq6;
          fpair = factor_lj * value;
        } else {
          rsq_lookup.f = rsq;
          itable = rsq_lookup.i & tb->nmask;
          itable >>= tb->nshiftbits;
          fraction = (rsq_lookup.f - tb->rsq[itable]) * tb->drsq[itable];
          value = tb->f[itable] + fraction*tb->df[itable];
          fpair = factor_lj * value;
        }

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;

        if (eflag) {
          if (tabstyle == LOOKUP)
            evdwl = tb->e[itable];
          else if (tabstyle == LINEAR || tabstyle == BITMAP)
            evdwl = tb->e[itable] + fraction*tb->de[itable];
          else
            evdwl = a * tb->e[itable] + b * tb->e[itable+1] +
              ((a*a*a-a)*tb->e2[itable] + (b*b*b-b)*tb->e2[itable+1]) *
              tb->deltasq6;
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally_full(i,evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }
}
