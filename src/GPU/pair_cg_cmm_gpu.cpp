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
#include "pair_cg_cmm_gpu.h"
#include "lmptype.h"
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

bool cmm_gpu_init(const int ntypes, double **cutsq, int **cg_types, 
                  double **host_lj1, double **host_lj2, double **host_lj3,
                  double **host_lj4, double **offset, double *special_lj,
                  const int nlocal, const int nall, const int max_nbors,
                  const int maxspecial, const double cell_size, int &gpu_mode,
                  FILE *screen);
void cmm_gpu_clear();
int * cmm_gpu_compute_n(const int timestep, const int ago, const int inum,
	 	        const int nall, double **host_x, int *host_type, 
                        double *boxlo, double *boxhi, int *tag, int **nspecial,
                        int **special, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        const double cpu_time, bool &success);
void cmm_gpu_compute(const int timestep, const int ago, const int inum,
	 	     const int nall, double **host_x, int *host_type,
                     int *ilist, int *numj, int **firstneigh,
		     const bool eflag, const bool vflag, const bool eatom,
                     const bool vatom, int &host_start, const double cpu_time,
                     bool &success);
double cmm_gpu_bytes();

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairCGCMMGPU::PairCGCMMGPU(LAMMPS *lmp) : PairCGCMM(lmp), gpu_mode(GPU_PAIR)
{
  respa_enable = 0;
  cpu_time = 0.0;
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairCGCMMGPU::~PairCGCMMGPU()
{
  cmm_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairCGCMMGPU::compute(int eflag, int vflag)
{
  int ntimestep = static_cast<int>(update->ntimestep % MAXSMALLINT);

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;
  
  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;
  
  bool success = true;
  
  if (gpu_mode == GPU_NEIGH) {
    inum = atom->nlocal;
    gpulist = cmm_gpu_compute_n(ntimestep, neighbor->ago, inum, nall,
			        atom->x, atom->type, domain->sublo,
				domain->subhi, atom->tag, atom->nspecial,
                                atom->special, eflag, vflag, eflag_atom,
                                vflag_atom, host_start, cpu_time, success);
  } else {
    inum = list->inum;
    cmm_gpu_compute(ntimestep, neighbor->ago, inum, nall, atom->x,
		    atom->type, list->ilist, list->numneigh, list->firstneigh,
		    eflag, vflag, eflag_atom, vflag_atom, host_start, cpu_time,
                    success);
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

void PairCGCMMGPU::init_style()
{
  cut_respa = NULL;

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
  bool init_ok = cmm_gpu_init(atom->ntypes+1,cutsq,cg_type,lj1,lj2,lj3,lj4,
                              offset, force->special_lj, atom->nlocal,
                              atom->nlocal+atom->nghost, 300, maxspecial,
                              cell_size, gpu_mode, screen);
  if (!init_ok)
    error->one("Insufficient memory on accelerator (or no fix gpu).\n"); 

  if (force->newton_pair) 
    error->all("Cannot use newton pair with GPU CGCMM pair style");

  if (gpu_mode != GPU_NEIGH) {
    int irequest = neighbor->request(this);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
  }
}

/* ---------------------------------------------------------------------- */

double PairCGCMMGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + cmm_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairCGCMMGPU::cpu_compute(int start, int eflag, int vflag) {
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,forcelj,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double *special_lj = force->special_lj;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

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

      if (j < nall) factor_lj = 1.0;
      else {
	factor_lj = special_lj[j/nall];
	j %= nall;
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        const int cgt=cg_type[itype][jtype];
        r2inv = 1.0/rsq;

	fpair = factor_lj;
	if (eflag) evdwl = factor_lj;
	if (cgt == CG_LJ12_4) {
	  const double r4inv = r2inv*r2inv;
	  fpair *= r4inv*(lj1[itype][jtype]*r4inv*r4inv
			  - lj2[itype][jtype]);
	  if (eflag) {
	    evdwl *= r4inv*(lj3[itype][jtype]*r4inv*r4inv
			    - lj4[itype][jtype]) - offset[itype][jtype];
	  }
	} else if (cgt == CG_LJ9_6) {
	  const double r3inv = r2inv*sqrt(r2inv);
	  const double r6inv = r3inv*r3inv;
	  fpair *= r6inv*(lj1[itype][jtype]*r3inv
			  - lj2[itype][jtype]);
	  if (eflag) {
	    evdwl *= r6inv*(lj3[itype][jtype]*r3inv
			    - lj4[itype][jtype]) - offset[itype][jtype];
	  }
	} else {
	  const double r6inv = r2inv*r2inv*r2inv;
	  fpair *= r6inv*(lj1[itype][jtype]*r6inv
			  - lj2[itype][jtype]);
	  if (eflag) {
	    evdwl *= r6inv*(lj3[itype][jtype]*r6inv
			    - lj4[itype][jtype]) - offset[itype][jtype];
	  }
	}
        fpair *= r2inv;

	f[i][0] += delx*fpair;
	f[i][1] += dely*fpair;
	f[i][2] += delz*fpair;

	if (evflag) ev_tally_full(i,evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairCGCMMGPU::cpu_compute(int *nbors, int start, int eflag, int vflag) {
  int i,j,itype,jtype;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int stride = nlocal-start;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,forcelj,factor_lj;
  double *special_lj = force->special_lj;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;

  // loop over neighbors of my atoms

  for (i = start; i < nlocal; i++) {
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    int *nbor = nbors + i - start;
    int jnum = *nbor;
    nbor += stride;
    int *nbor_end = nbor + stride * jnum;

    for (; nbor<nbor_end; nbor+=stride) {
      j = *nbor;
      
      if (j < nall) factor_lj = 1.0;
      else {
	factor_lj = special_lj[j/nall];
	j %= nall;
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        const int cgt=cg_type[itype][jtype];
        r2inv = 1.0/rsq;

	fpair = factor_lj;
	if (eflag) evdwl = factor_lj;
	if (cgt == CG_LJ12_4) {
	  const double r4inv = r2inv*r2inv;
	  fpair *= r4inv*(lj1[itype][jtype]*r4inv*r4inv
			  - lj2[itype][jtype]);
	  if (eflag) {
	    evdwl *= r4inv*(lj3[itype][jtype]*r4inv*r4inv
			    - lj4[itype][jtype]) - offset[itype][jtype];
	  }
	} else if (cgt == CG_LJ9_6) {
	  const double r3inv = r2inv*sqrt(r2inv);
	  const double r6inv = r3inv*r3inv;
	  fpair *= r6inv*(lj1[itype][jtype]*r3inv
			  - lj2[itype][jtype]);
	  if (eflag) {
	    evdwl *= r6inv*(lj3[itype][jtype]*r3inv
			    - lj4[itype][jtype]) - offset[itype][jtype];
	  }
	} else {
	  const double r6inv = r2inv*r2inv*r2inv;
	  fpair *= r6inv*(lj1[itype][jtype]*r6inv
			  - lj2[itype][jtype]);
	  if (eflag) {
	    evdwl *= r6inv*(lj3[itype][jtype]*r6inv
			    - lj4[itype][jtype]) - offset[itype][jtype];
	  }
	}
        fpair *= r2inv;

	f[i][0] += delx*fpair;
	f[i][1] += dely*fpair;
	f[i][2] += delz*fpair;

        if (j<start) {
  	  if (evflag) ev_tally_full(i,evdwl,0.0,fpair,delx,dely,delz);
        } else {
          if (j<nlocal) {
	    f[j][0] -= delx*fpair;
	    f[j][1] -= dely*fpair;
	    f[j][2] -= delz*fpair;
  	  }
	  if (evflag) ev_tally(i,j,nlocal,0,
			       evdwl,0.0,fpair,delx,dely,delz);
	}
      }
    }
  }
}

