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

#include "math.h"
#include "pair_dpd_tstat_omp.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "neigh_list.h"
#include "comm.h"
#include "random_mars.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define EPSILON 1.0e-20

/* ---------------------------------------------------------------------- */

PairDPDTstatOMP::PairDPDTstatOMP(LAMMPS *lmp) : PairDPDOMP(lmp)
{
  single_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairDPDTstatOMP::compute(int eflag, int vflag)
{
  int i,j;
  
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
    ev_setup_thr(eflag,vflag);
  } else evflag = vflag_fdotr = 0;

  // adjust sigma if target T is changing

  if (t_start != t_stop) {
    double delta = update->ntimestep - update->beginstep;
    delta /= update->endstep - update->beginstep;
    temperature = t_start + delta * (t_stop-t_start);
    double boltz = force->boltz;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
	sigma[i][j] = sigma[j][i] = sqrt(2.0*boltz*temperature*gamma[i][j]);
  }

  if (evflag) {
    if (eflag) {
      if (force->newton_pair) return eval_tstat<1,1,1>();
      else return eval_tstat<1,1,0>();
    } else {
      if (force->newton_pair) return eval_tstat<1,0,1>();
      else return eval_tstat<1,0,0>();
    }
  } else {
    if (force->newton_pair) return eval_tstat<0,0,1>();
    else return eval_tstat<0,0,0>();
  }
}


template <int EVFLAG, int EFLAG, int NEWTON_PAIR> 
void PairDPDTstatOMP::eval_tstat() 
{

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
    int i,j,ii,jj,inum,jnum,itype,jtype,tid;
    double xtmp,ytmp,ztmp,delx,dely,delz,fpair;
    double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
    double rsq,r,rinv,dot,wd,randnum,factor_dpd;
    int *ilist,*jlist,*numneigh,**firstneigh;

    int nlocal = atom->nlocal;
    int nall = nlocal + atom->nghost;
    int nthreads = comm->nthreads;

    double **x = atom->x;
    double **v = atom->v;
    int *type = atom->type;
    double *special_lj = force->special_lj;
    double dtinvsqrt = 1.0/sqrt(update->dt);
    double fxtmp,fytmp,fztmp;

    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
  
    // loop over neighbors of my atoms
    int iifrom, iito;
    double **f = loop_setup_thr(atom->f,iifrom,iito,tid,inum,nall,nthreads);
    RanMars &rng = *random[tid];
    for (ii = iifrom; ii < iito; ++ii) {

      i = ilist[ii];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      vxtmp = v[i][0];
      vytmp = v[i][1];
      vztmp = v[i][2];
      itype = type[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];
      fxtmp=fytmp=fztmp=0.0;

      for (jj = 0; jj < jnum; jj++) {
	j = jlist[jj];

	if (j < nall) factor_dpd = 1.0;
	else {
	  factor_dpd = special_lj[j/nall];
	  j %= nall;
	}

	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;
	jtype = type[j];

	if ((rsq < cutsq[itype][jtype]) && (rsq > EPSILON)) { // r can be 0.0 in DPD systems
	  r = sqrt(rsq);
	  rinv = 1.0/r;
	  delvx = vxtmp - v[j][0];
	  delvy = vytmp - v[j][1];
	  delvz = vztmp - v[j][2];
	  dot = delx*delvx + dely*delvy + delz*delvz;
	  wd = 1.0 - r/cut[itype][jtype];
	  randnum = rng.gaussian();

	  // drag force = -gamma * wd^2 * (delx dot delv) / r
	  // random force = sigma * wd * rnd * dtinvsqrt;

	  fpair = -gamma[itype][jtype]*wd*wd*dot*rinv;
	  fpair += sigma[itype][jtype]*wd*randnum*dtinvsqrt;
	  fpair *= factor_dpd*rinv;	

	  fxtmp += delx*fpair;
	  fytmp += dely*fpair;
	  fztmp += delz*fpair;
	  if (NEWTON_PAIR || j < nlocal) {
	    f[j][0] -= delx*fpair;
	    f[j][1] -= dely*fpair;
	    f[j][2] -= delz*fpair;
	  }

	  if (EVFLAG) ev_tally_thr(i,j,nlocal,NEWTON_PAIR,
				   0.0,0.0,fpair,delx,dely,delz,tid);
	}
      }
      f[i][0] += fxtmp;
      f[i][1] += fytmp;
      f[i][2] += fztmp;
    }

    // reduce per thread forces into global force array.
    force_reduce_thr(atom->f,nall,nthreads,tid);
  }
  if (EVFLAG) ev_reduce_thr();
  if (vflag_fdotr) virial_compute();
}

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void PairDPDTstatOMP::settings(int narg, char **arg)
{
  if (narg != 4) error->all("Illegal pair_style command");

  t_start = force->numeric(arg[0]);
  t_stop = force->numeric(arg[1]);
  cut_global = force->numeric(arg[2]);
  seed = force->inumeric(arg[3]);

  temperature = t_start;

  // initialize Marsaglia RNG with processor-unique seed

  if (seed <= 0) error->all("Illegal pair_style command");

  // we need to set up individual RNGs for each thread.
  int nthreads = comm->nthreads;
  int tid;
  if (random) {
    for (tid=0; tid < nthreads; ++tid) 
      delete random[tid];
    delete[] random;
  }
  random = new RanMars*[nthreads];
  for (tid=0; tid < nthreads; ++tid) 
    random[tid] = new RanMars(lmp,seed + nthreads*comm->me + tid);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
	if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairDPDTstatOMP::coeff(int narg, char **arg)
{
  if (narg < 3 || narg > 4) error->all("Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double a0_one = 0.0;
  double gamma_one = force->numeric(arg[2]);

  double cut_one = cut_global;
  if (narg == 4) cut_one = force->numeric(arg[3]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      a0[i][j] = a0_one;
      gamma[i][j] = gamma_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDPDTstatOMP::write_restart_settings(FILE *fp)
{
  fwrite(&t_start,sizeof(double),1,fp);
  fwrite(&t_stop,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&seed,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDPDTstatOMP::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&t_start,sizeof(double),1,fp);
    fread(&t_stop,sizeof(double),1,fp);
    fread(&cut_global,sizeof(double),1,fp);
    fread(&seed,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&t_start,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&t_stop,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&seed,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);

  // initialize Marsaglia RNG with processor-unique seed
  // same seed that pair_style command initially specified

  // we need to set up individual RNGs for each thread.
  int nthreads = comm->nthreads;
  int tid;
  if (random) {
    for (tid=0; tid < nthreads; ++tid) 
      delete random[tid];
    delete[] random;
  }
  random = new RanMars*[nthreads];
  for (tid=0; tid < nthreads; ++tid) 
    random[tid] = new RanMars(lmp,seed + nthreads*comm->me + tid);

}

