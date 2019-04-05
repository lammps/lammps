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
   Contributing author: Trung Dac Nguyen (ORNL)
------------------------------------------------------------------------- */

#include "lmptype.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_zbl_gpu.h"
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

int zbl_gpu_init(const int ntypes, double **cutsq, double **host_sw1,
                 double **host_sw2, double **host_sw3, double **host_sw4,
                 double **host_sw5, double **host_d1a, double **host_d2a,
                 double **host_d3a, double **host_d4a, double **host_zze,
                 double cut_globalsq, double cut_innersq, double cut_inner,
                 const int inum, const int nall, const int max_nbors,
                 const int maxspecial, const double cell_size,
                 int &gpu_mode, FILE *screen);
void zbl_gpu_clear();
int ** zbl_gpu_compute_n(const int ago, const int inum,
                         const int nall, double **host_x, int *host_type,
                         double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag,
                         const bool eatom, const bool vatom, int &host_start,
                         int **ilist, int **jnum,
                         const double cpu_time, bool &success);
void zbl_gpu_compute(const int ago, const int inum, const int nall,
                     double **host_x, int *host_type, int *ilist, int *numj,
                     int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, int &host_start,
                     const double cpu_time, bool &success);
double zbl_gpu_bytes();

/* ---------------------------------------------------------------------- */

PairZBLGPU::PairZBLGPU(LAMMPS *lmp) : PairZBL(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  reinitflag = 0;
  cpu_time = 0.0;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairZBLGPU::~PairZBLGPU()
{
  zbl_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairZBLGPU::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    inum = atom->nlocal;
    firstneigh = zbl_gpu_compute_n(neighbor->ago, inum, nall,
                                   atom->x, atom->type, domain->sublo,
                                   domain->subhi, atom->tag, atom->nspecial,
                                   atom->special, eflag, vflag, eflag_atom,
                                   vflag_atom, host_start,
                                   &ilist, &numneigh, cpu_time, success);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    zbl_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type,
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

void PairZBLGPU::init_style()
{
  if (force->newton_pair)
    error->all(FLERR,"Cannot use newton pair with zbl/gpu pair style");

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

  cut_innersq = cut_inner * cut_inner;
  cut_globalsq = cut_global * cut_global;

  int maxspecial=0;
  if (atom->molecular)
    maxspecial=atom->maxspecial;
  int success = zbl_gpu_init(atom->ntypes+1, cutsq, sw1, sw2, sw3, sw4,
                             sw5, d1a, d2a, d3a, d4a, zze,
                             cut_globalsq, cut_innersq, cut_inner,
                             atom->nlocal, atom->nlocal+atom->nghost,
                             300, maxspecial, cell_size, gpu_mode, screen);
  GPU_EXTRA::check_flag(success,error,world);

  if (gpu_mode == GPU_FORCE) {
    int irequest = neighbor->request(this);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
  }
}

/* ---------------------------------------------------------------------- */

double PairZBLGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + zbl_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairZBLGPU::cpu_compute(int start, int inum, int eflag, int /* vflag */,
                             int *ilist, int *numneigh, int **firstneigh) {
  int i,j,ii,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,t,fswitch,eswitch;
  int *jlist;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;

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
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cut_globalsq) {
              r = sqrt(rsq);
        fpair = dzbldr(r, itype, jtype);

              if (rsq > cut_innersq) {
                t = r - cut_inner;
                fswitch = t*t *
                  (sw1[itype][jtype] + sw2[itype][jtype]*t);
                fpair += fswitch;
              }

        fpair *= -1.0/r;
        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;

        if (eflag) {
          evdwl = e_zbl(r, itype, jtype);
                evdwl += sw5[itype][jtype];
                if (rsq > cut_innersq) {
                  eswitch = t*t*t *
                    (sw3[itype][jtype] + sw4[itype][jtype]*t);
                  evdwl += eswitch;
                }
        }

        if (evflag) ev_tally_full(i,evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }
}
