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
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_lj_sdk_gpu.h"
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

using namespace LAMMPS_NS;

// External functions from cuda library for atom decomposition

int sdk_gpu_init(const int ntypes, double **cutsq, int **cg_types,
                 double **host_lj1, double **host_lj2, double **host_lj3,
                 double **host_lj4, double **offset, double *special_lj,
                 const int nlocal, const int nall, const int max_nbors,
                 const int maxspecial, const double cell_size, int &gpu_mode,
                 FILE *screen);
void sdk_gpu_clear();
int ** sdk_gpu_compute_n(const int ago, const int inum, const int nall,
                         double **host_x, int *host_type, double *sublo,
                         double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag,
                         const bool eatom, const bool vatom, int &host_start,
                         int **ilist, int **jnum,
                         const double cpu_time, bool &success);
void sdk_gpu_compute(const int ago, const int inum, const int nall,
                     double **host_x, int *host_type, int *ilist, int *numj,
                     int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, int &host_start,
                     const double cpu_time, bool &success);
double sdk_gpu_bytes();

#include "lj_sdk_common.h"

using namespace LJSDKParms;

/* ---------------------------------------------------------------------- */

PairLJSDKGPU::PairLJSDKGPU(LAMMPS *lmp) : PairLJSDK(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  reinitflag = 0;
  cpu_time = 0.0;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairLJSDKGPU::~PairLJSDKGPU()
{
  sdk_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairLJSDKGPU::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    inum = atom->nlocal;
    firstneigh = sdk_gpu_compute_n(neighbor->ago, inum, nall, atom->x,
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
    sdk_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type,
                    ilist, numneigh, firstneigh, eflag, vflag, eflag_atom,
                    vflag_atom, host_start, cpu_time, success);
  }
  if (!success)
    error->one(FLERR,"Insufficient memory on accelerator");

  if (host_start<inum) {
    cpu_time = MPI_Wtime();
    if (evflag) {
      if (eflag) cpu_compute<1,1>(host_start, inum, ilist, numneigh, firstneigh);
      else cpu_compute<1,0>(host_start, inum, ilist, numneigh, firstneigh);
    } else cpu_compute<0,0>(host_start, inum, ilist, numneigh, firstneigh);
    cpu_time = MPI_Wtime() - cpu_time;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJSDKGPU::init_style()
{
  if (force->newton_pair)
    error->all(FLERR,"Cannot use newton pair with lj/sdk/gpu pair style");

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

  int maxspecial=0;
  if (atom->molecular)
    maxspecial=atom->maxspecial;
  int success = sdk_gpu_init(atom->ntypes+1,cutsq,lj_type,lj1,lj2,lj3,lj4,
                             offset, force->special_lj, atom->nlocal,
                             atom->nlocal+atom->nghost, 300, maxspecial,
                             cell_size, gpu_mode, screen);
  GPU_EXTRA::check_flag(success,error,world);

  if (gpu_mode == GPU_FORCE) {
    int irequest = neighbor->request(this,instance_me);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
  }
}

/* ---------------------------------------------------------------------- */

double PairLJSDKGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + sdk_gpu_bytes();
}

/* ---------------------------------------------------------------------- */
template <int EVFLAG, int EFLAG>
void PairLJSDKGPU::cpu_compute(int start, int inum, int *ilist,
                               int *numneigh, int **firstneigh)
{
  int i,j,ii,jj,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,forcelj,factor_lj;

  const double * const * const x = atom->x;
  double * const * const f = atom->f;
  const int * const type = atom->type;
  const double * const special_lj = force->special_lj;
  double fxtmp,fytmp,fztmp;
  evdwl=0.0;

  // loop over neighbors of my atoms

  for (ii = start; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fxtmp=fytmp=fztmp=0.0;

    const int itype = type[i];
    const int * const jlist = firstneigh[i];
    const int jnum = numneigh[i];

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
        r2inv = 1.0/rsq;
        const int ljt = lj_type[itype][jtype];

        if (ljt == LJ12_4) {
          const double r4inv=r2inv*r2inv;
          forcelj = r4inv*(lj1[itype][jtype]*r4inv*r4inv
                           - lj2[itype][jtype]);

          if (EFLAG)
            evdwl = r4inv*(lj3[itype][jtype]*r4inv*r4inv
                           - lj4[itype][jtype]) - offset[itype][jtype];

        } else if (ljt == LJ9_6) {
          const double r3inv = r2inv*sqrt(r2inv);
          const double r6inv = r3inv*r3inv;
          forcelj = r6inv*(lj1[itype][jtype]*r3inv
                           - lj2[itype][jtype]);
          if (EFLAG)
            evdwl = r6inv*(lj3[itype][jtype]*r3inv
                           - lj4[itype][jtype]) - offset[itype][jtype];

        } else if (ljt == LJ12_6) {
          const double r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv*(lj1[itype][jtype]*r6inv
                          - lj2[itype][jtype]);
          if (EFLAG)
            evdwl = r6inv*(lj3[itype][jtype]*r6inv
                           - lj4[itype][jtype]) - offset[itype][jtype];
        } else continue;

        fpair = factor_lj*forcelj*r2inv;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;

        if (EVFLAG) ev_tally_full(i,evdwl,0.0,fpair,delx,dely,delz);
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}
