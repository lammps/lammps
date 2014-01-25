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
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "pair_colloid_gpu.h"
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
#include "string.h"
#include "gpu_extra.h"

using namespace LAMMPS_NS;

// External functions from cuda library for atom decomposition

int colloid_gpu_init(const int ntypes, double **cutsq, double **host_lj1,
                     double **host_lj2, double **host_lj3, double **host_lj4, 
                     double **offset, double *special_lj, double **host_a12, 
                     double **host_a1, double **host_a2, double **host_d1, 
                     double **host_d2, double **host_sigma3, double **host_sigma6, 
                     int **host_form, const int nlocal, 
                     const int nall, const int max_nbors, const int maxspecial,
                     const double cell_size, int &gpu_mode, FILE *screen);
void colloid_gpu_clear();
int ** colloid_gpu_compute_n(const int ago, const int inum,
                             const int nall, double **host_x, int *host_type, 
                             double *sublo, double *subhi, tagint *tag, int **nspecial,
                             tagint **special, const bool eflag, const bool vflag,
                             const bool eatom, const bool vatom, int &host_start,
                             int **ilist, int **jnum,
                             const double cpu_time, bool &success);
void colloid_gpu_compute(const int ago, const int inum, const int nall, 
                         double **host_x, int *host_type, int *ilist, int *numj,
                         int **firstneigh, const bool eflag, const bool vflag,
                         const bool eatom, const bool vatom, int &host_start,
                         const double cpu_time, bool &success);
double colloid_gpu_bytes();

/* ---------------------------------------------------------------------- */

PairColloidGPU::PairColloidGPU(LAMMPS *lmp) : PairColloid(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  cpu_time = 0.0;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error); 
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairColloidGPU::~PairColloidGPU()
{
  colloid_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairColloidGPU::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;
  
  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;
  
  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    inum = atom->nlocal;
    firstneigh = colloid_gpu_compute_n(neighbor->ago, inum, nall,
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
    colloid_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type,
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

void PairColloidGPU::init_style()
{
  if (force->newton_pair) 
    error->all(FLERR,"Cannot use newton pair with colloid/gpu pair style");

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

  int **_form = NULL;
  int n=atom->ntypes;
  memory->create(_form,n+1,n+1,"colloid/gpu:_form");
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      if (form[i][j] == SMALL_SMALL) _form[i][j] = 0;
      else if (form[i][j] == SMALL_LARGE) _form[i][j] = 1;
      else if (form[i][j] == LARGE_LARGE) _form[i][j] = 2;
    }
  }
  int maxspecial=0;
  if (atom->molecular)
    maxspecial=atom->maxspecial;
  int success = colloid_gpu_init(atom->ntypes+1, cutsq, lj1, lj2, lj3, lj4,
                                 offset, force->special_lj, a12, a1, a2, 
                                 d1, d2, sigma3, sigma6, _form, atom->nlocal,
                                 atom->nlocal+atom->nghost, 300, maxspecial,
                                 cell_size, gpu_mode, screen);
  memory->destroy(_form);
  GPU_EXTRA::check_flag(success,error,world);

  if (gpu_mode == GPU_FORCE) {
    int irequest = neighbor->request(this);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
  }
}

/* ---------------------------------------------------------------------- */

double PairColloidGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + colloid_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairColloidGPU::cpu_compute(int start, int inum, int eflag, int vflag, 
                                 int *ilist, int *numneigh, int **firstneigh) 
{
  int i,j,ii,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double r,rsq,r2inv,r6inv,forcelj,factor_lj;
  double c1,c2,fR,dUR,dUA;
  double K[9],h[4],g[4];
  int *jlist;

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

      if (rsq >= cutsq[itype][jtype]) continue;

      switch (form[itype][jtype]) {
      case SMALL_SMALL: 
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
        fpair = factor_lj*forcelj*r2inv;
        if (eflag) 
          evdwl = r6inv*(r6inv*lj3[itype][jtype]-lj4[itype][jtype]) -
            offset[itype][jtype];
        break;
	      
      case SMALL_LARGE:
        c2 = a2[itype][jtype];
        K[1] = c2*c2;
        K[2] = rsq;
        K[0] = K[1] - rsq;
        K[4] = rsq*rsq;
        K[3] = K[1] - K[2];
        K[3] *= K[3]*K[3];
        K[6] = K[3]*K[3];
        fR = sigma3[itype][jtype]*a12[itype][jtype]*c2*K[1]/K[3];
        fpair = 4.0/15.0*fR*factor_lj * 
          (2.0*(K[1]+K[2]) * (K[1]*(5.0*K[1]+22.0*K[2])+5.0*K[4]) * 
          sigma6[itype][jtype]/K[6]-5.0) / K[0];
        if (eflag) 
          evdwl = 2.0/9.0*fR * 
            (1.0-(K[1]*(K[1]*(K[1]/3.0+3.0*K[2])+4.2*K[4])+K[2]*K[4]) *
            sigma6[itype][jtype]/K[6]) - offset[itype][jtype];
        if (rsq <= K[1]) 
          error->one(FLERR,"Overlapping small/large in pair colloid");
        break;

      case LARGE_LARGE:
        r = sqrt(rsq);
        c1 = a1[itype][jtype];
        c2 = a2[itype][jtype];
        K[0] = c1*c2;
        K[1] = c1+c2;
        K[2] = c1-c2;
        K[3] = K[1]+r;
        K[4] = K[1]-r;
        K[5] = K[2]+r;
        K[6] = K[2]-r;
        K[7] = 1.0/(K[3]*K[4]);
        K[8] = 1.0/(K[5]*K[6]);
        g[0] = pow(K[3],-7.0);
        g[1] = pow(K[4],-7.0);
        g[2] = pow(K[5],-7.0);
        g[3] = pow(K[6],-7.0);
        h[0] = ((K[3]+5.0*K[1])*K[3]+30.0*K[0])*g[0];
        h[1] = ((K[4]+5.0*K[1])*K[4]+30.0*K[0])*g[1];
        h[2] = ((K[5]+5.0*K[2])*K[5]-30.0*K[0])*g[2];
        h[3] = ((K[6]+5.0*K[2])*K[6]-30.0*K[0])*g[3];
        g[0] *= 42.0*K[0]/K[3]+6.0*K[1]+K[3];
        g[1] *= 42.0*K[0]/K[4]+6.0*K[1]+K[4];
        g[2] *= -42.0*K[0]/K[5]+6.0*K[2]+K[5];
        g[3] *= -42.0*K[0]/K[6]+6.0*K[2]+K[6];

        fR = a12[itype][jtype]*sigma6[itype][jtype]/r/37800.0;
        evdwl = fR * (h[0]-h[1]-h[2]+h[3]);
        dUR = evdwl/r + 5.0*fR*(g[0]+g[1]-g[2]-g[3]);
        dUA = -a12[itype][jtype]/3.0*r*((2.0*K[0]*K[7]+1.0)*K[7] + 
          (2.0*K[0]*K[8]-1.0)*K[8]);
        fpair = factor_lj * (dUR+dUA)/r;
        if (eflag)
          evdwl += a12[itype][jtype]/6.0 * 
            (2.0*K[0]*(K[7]+K[8])-log(K[8]/K[7])) - offset[itype][jtype];
        if (r <= K[1]) 
          error->one(FLERR,"Overlapping large/large in pair colloid");
        break;
      }
      
      if (eflag) evdwl *= factor_lj;

      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;

      if (evflag) ev_tally_full(i,evdwl,0.0,fpair,delx,dely,delz);
    }
  }
}
