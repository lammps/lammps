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

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "pair_lj_cut_coul_cut_gpu.h"
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

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

// External functions from cuda library for atom decomposition

bool ljc_gpu_init(const int ntypes, double **cutsq, double **host_lj1,
                  double **host_lj2, double **host_lj3, double **host_lj4, 
                  double **offset, double *special_lj, const int nlocal, 
                  const int nall, const int max_nbors, const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen,
                  double **host_cut_ljsq, double **host_cut_coulsq,
                  double *host_special_coul, const double qqrd2e);
void ljc_gpu_clear();
int * ljc_gpu_compute_n(const int timestep, const int ago, const int inum,
	 	        const int nall, double **host_x, int *host_type, 
                        double *boxlo, double *boxhi, int *tag, int **nspecial,
                        int **special, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        const double cpu_time, bool &success, double *host_q);
void ljc_gpu_compute(const int timestep, const int ago, const int inum,
	 	     const int nall, double **host_x, int *host_type,
                     int *ilist, int *numj, int **firstneigh,
		     const bool eflag, const bool vflag, const bool eatom,
                     const bool vatom, int &host_start, const double cpu_time,
                     bool &success, double *host_q);
double ljc_gpu_bytes();

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJCutCoulCutGPU::PairLJCutCoulCutGPU(LAMMPS *lmp) : PairLJCutCoulCut(lmp), gpu_mode(GPU_PAIR)
{
  respa_enable = 0;
  cpu_time = 0.0;
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairLJCutCoulCutGPU::~PairLJCutCoulCutGPU()
{
  ljc_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulCutGPU::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;
  
  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;
  
  bool success = true;
  
  if (gpu_mode == GPU_NEIGH) {
    inum = atom->nlocal;
    gpulist = ljc_gpu_compute_n(update->ntimestep, neighbor->ago, inum, nall,
			        atom->x, atom->type, domain->sublo,
				domain->subhi, atom->tag, atom->nspecial,
                                atom->special, eflag, vflag, eflag_atom,
                                vflag_atom, host_start, cpu_time, success,
                                atom->q);
  } else {
    inum = list->inum;
    ljc_gpu_compute(update->ntimestep, neighbor->ago, inum, nall, atom->x,
		    atom->type, list->ilist, list->numneigh, list->firstneigh,
		    eflag, vflag, eflag_atom, vflag_atom, host_start, cpu_time,
                    success, atom->q);
  }
  if (!success)
    error->one("Out of memory on GPGPU");

  if (host_start<inum) {
    cpu_time = MPI_Wtime();
    if (gpu_mode == GPU_NEIGH)
      cpu_compute(gpulist, host_start, eflag, vflag);
    else
      cpu_compute(host_start, eflag, vflag);
    cpu_time = MPI_Wtime() - cpu_time;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCutCoulCutGPU::init_style()
{
  if (!atom->q_flag)
    error->all("Pair style lj/cut/coul/cut requires atom attribute q");
  if (force->pair_match("gpu",0) == NULL)
    error->all("Cannot use pair hybrid with multiple GPU pair styles");

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
  bool init_ok = ljc_gpu_init(atom->ntypes+1, cutsq, lj1, lj2, lj3, lj4,
                              offset, force->special_lj, atom->nlocal,
                              atom->nlocal+atom->nghost, 300, maxspecial,
                              cell_size, gpu_mode, screen, cut_ljsq, cut_coulsq,
                              force->special_coul, force->qqrd2e);
  if (!init_ok)
    error->one("Insufficient memory on accelerator (or no fix gpu).\n"); 

  if (force->newton_pair) 
    error->all("Cannot use newton pair with GPU LJ pair style");

  if (gpu_mode != GPU_NEIGH) {
    int irequest = neighbor->request(this);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
  }
}

/* ---------------------------------------------------------------------- */

double PairLJCutCoulCutGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + ljc_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulCutGPU::cpu_compute(int start, int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double rsq,r2inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = start; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      if (j < nall) factor_coul = factor_lj = 1.0;
      else {
	factor_coul = special_coul[j/nall];
	factor_lj = special_lj[j/nall];
	j %= nall;
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
	r2inv = 1.0/rsq;

	if (rsq < cut_coulsq[itype][jtype])
	  forcecoul = qqrd2e * qtmp*q[j]*sqrt(r2inv);
	else forcecoul = 0.0;

	if (rsq < cut_ljsq[itype][jtype]) {
	  r6inv = r2inv*r2inv*r2inv;
	  forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	} else forcelj = 0.0;

	fpair = (factor_coul*forcecoul + factor_lj*forcelj) * r2inv;

	f[i][0] += delx*fpair;
	f[i][1] += dely*fpair;
	f[i][2] += delz*fpair;

	if (eflag) {
	  if (rsq < cut_coulsq[itype][jtype])
	    ecoul = factor_coul * qqrd2e * qtmp*q[j]*sqrt(r2inv);
	  else ecoul = 0.0;
	  if (rsq < cut_ljsq[itype][jtype]) {
	    evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
	      offset[itype][jtype];
	    evdwl *= factor_lj;
	  } else evdwl = 0.0;
	}

	if (evflag) ev_tally_full(i,evdwl,ecoul,fpair,delx,dely,delz);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulCutGPU::cpu_compute(int *nbors, int start, int eflag,
                                      int vflag)
{
  int i,j,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double rsq,r2inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;

  evdwl = ecoul = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int stride = nlocal-start;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  double qqrd2e = force->qqrd2e;

  // loop over neighbors of my atoms

  for (i = start; i < nlocal; i++) {
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    int *nbor = nbors + i - start;
    jnum = *nbor;
    nbor += stride;
    int *nbor_end = nbor + stride * jnum;

    for (; nbor<nbor_end; nbor+=stride) {
      j = *nbor;

      if (j < nall) factor_coul = factor_lj = 1.0;
      else {
	factor_coul = special_coul[j/nall];
	factor_lj = special_lj[j/nall];
	j %= nall;
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
	r2inv = 1.0/rsq;

	if (rsq < cut_coulsq[itype][jtype])
	  forcecoul = qqrd2e * qtmp*q[j]*sqrt(r2inv);
	else forcecoul = 0.0;

	if (rsq < cut_ljsq[itype][jtype]) {
	  r6inv = r2inv*r2inv*r2inv;
	  forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	} else forcelj = 0.0;

	fpair = (factor_coul*forcecoul + factor_lj*forcelj) * r2inv;

	f[i][0] += delx*fpair;
	f[i][1] += dely*fpair;
	f[i][2] += delz*fpair;

	if (eflag) {
	  if (rsq < cut_coulsq[itype][jtype])
	    ecoul = factor_coul * qqrd2e * qtmp*q[j]*sqrt(r2inv);
	  else ecoul = 0.0;
	  if (rsq < cut_ljsq[itype][jtype]) {
	    evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
	      offset[itype][jtype];
	    evdwl *= factor_lj;
	  } else evdwl = 0.0;
	}

        if (j<start) {
  	  if (evflag) ev_tally_full(i,evdwl,ecoul,fpair,delx,dely,delz);
        } else {
          if (j<nlocal) {
	    f[j][0] -= delx*fpair;
	    f[j][1] -= dely*fpair;
	    f[j][2] -= delz*fpair;
  	  }
	  if (evflag) ev_tally(i,j,nlocal,0,
			       evdwl,ecoul,fpair,delx,dely,delz);
        }
      }
    }
  }
}
