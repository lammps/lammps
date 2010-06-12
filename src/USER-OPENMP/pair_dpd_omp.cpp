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
   Contributing author: Kurt Smith (U Pittsburgh)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "pair_dpd_omp.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define EPSILON 1.0e-20

/* ---------------------------------------------------------------------- */

PairDPDOMP::PairDPDOMP(LAMMPS *lmp) : PairOMP(lmp)
{
  random = NULL;
}

/* ---------------------------------------------------------------------- */

PairDPDOMP::~PairDPDOMP()
{
  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);

    memory->destroy_2d_double_array(cut);
    memory->destroy_2d_double_array(a0);
    memory->destroy_2d_double_array(gamma);
    memory->destroy_2d_double_array(sigma);
    allocated = 0;
  }

  if (random) {
    for (int i=0; i < comm->nthreads; ++i)
      delete random[i];
    delete[] random;
    random=NULL;
  }
}

/* ---------------------------------------------------------------------- */

void PairDPDOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
    ev_setup_thr(eflag,vflag);
  } else evflag = vflag_fdotr = 0;

  if (evflag) {
    if (eflag) {
      if (force->newton_pair) return eval<1,1,1>();
      else return eval<1,1,0>();
    } else {
      if (force->newton_pair) return eval<1,0,1>();
      else return eval<1,0,0>();
    }
  } else {
    if (force->newton_pair) return eval<0,0,1>();
    else return eval<0,0,0>();
  }
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR> 
void PairDPDOMP::eval() 
{

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
    int i,j,ii,jj,inum,jnum,itype,jtype,tid;
    double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
    double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
    double rsq,r,rinv,dot,wd,randnum,factor_dpd;
    int *ilist,*jlist,*numneigh,**firstneigh;
    evdwl=0.0;

    double **x = atom->x;
    double **v = atom->v;
    int *type = atom->type;
    int nlocal = atom->nlocal;
    int nall = nlocal + atom->nghost;
    int nthreads = comm->nthreads;
    double *special_lj = force->special_lj;
    double dtinvsqrt = 1.0/sqrt(update->dt);
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
  
    // loop over neighbors of my atoms
    int iifrom, iito;
    double **f = loop_setup_thr(atom->f,iifrom,iito,tid,inum,nall,nthreads);
    RanMars *rng = random[tid];
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
	  randnum = rng->gaussian();

	  // conservative force = a0 * wd
	  // drag force = -gamma * wd^2 * (delx dot delv) / r
	  // random force = sigma * wd * rnd * dtinvsqrt;

	  fpair = a0[itype][jtype]*wd;
	  fpair -= gamma[itype][jtype]*wd*wd*dot*rinv;
	  fpair += sigma[itype][jtype]*wd*randnum*dtinvsqrt;
	  fpair *= factor_dpd*rinv;	

	  f[i][0] += delx*fpair;
	  f[i][1] += dely*fpair;
	  f[i][2] += delz*fpair;
	  if (NEWTON_PAIR || j < nlocal) {
	    f[j][0] -= delx*fpair;
	    f[j][1] -= dely*fpair;
	    f[j][2] -= delz*fpair;
	  }

	  if (EFLAG) {
	    // unshifted eng of conservative term:
	    // evdwl = -a0[itype][jtype]*r * (1.0-0.5*r/cut[itype][jtype]);
	    // eng shifted to 0.0 at cutoff
	    evdwl = 0.5*a0[itype][jtype]*cut[itype][jtype] * wd*wd;
	    evdwl *= factor_dpd;
	  }

	  if (EVFLAG) ev_tally_thr(i,j,nlocal,NEWTON_PAIR,evdwl,0.0,
				   fpair,delx,dely,delz,tid);
	}
      }
    }

    // reduce per thread forces into global force array.
    force_reduce_thr(atom->f,nall,nthreads,tid);
  }
  ev_reduce_thr();
  if (vflag_fdotr) virial_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays 
------------------------------------------------------------------------- */

void PairDPDOMP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");

  cut = memory->create_2d_double_array(n+1,n+1,"pair:cut");
  a0 = memory->create_2d_double_array(n+1,n+1,"pair:a0");
  gamma = memory->create_2d_double_array(n+1,n+1,"pair:gamma");
  sigma = memory->create_2d_double_array(n+1,n+1,"pair:sigma");
}

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void PairDPDOMP::settings(int narg, char **arg)
{
  if (narg != 3) error->all("Illegal pair_style command");

  temperature = force->numeric(arg[0]);
  cut_global = force->numeric(arg[1]);
  seed = force->inumeric(arg[2]);

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

void PairDPDOMP::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5) error->all("Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double a0_one = force->numeric(arg[2]);
  double gamma_one = force->numeric(arg[3]);

  double cut_one = cut_global;
  if (narg == 5) cut_one = force->numeric(arg[4]);

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
   init specific to this pair style
------------------------------------------------------------------------- */

void PairDPDOMP::init_style()
{
  if (comm->ghost_velocity == 0)
    error->all("Pair dpd requires ghost atoms store velocity");

  // if newton off, forces between atoms ij will be double computed
  // using different random numbers

  if (force->newton_pair == 0 && comm->me == 0) error->warning(
      "Pair dpd needs newton pair on for momentum conservation");

  int irequest = neighbor->request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairDPDOMP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all("All pair coeffs are not set");

  sigma[i][j] = sqrt(2.0*force->boltz*temperature*gamma[i][j]);
     
  cut[j][i] = cut[i][j];
  a0[j][i] = a0[i][j];
  gamma[j][i] = gamma[i][j];
  sigma[j][i] = sigma[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDPDOMP::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
	fwrite(&a0[i][j],sizeof(double),1,fp);
	fwrite(&gamma[i][j],sizeof(double),1,fp);
	fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDPDOMP::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
	if (me == 0) {
	  fread(&a0[i][j],sizeof(double),1,fp);
	  fread(&gamma[i][j],sizeof(double),1,fp);
	  fread(&cut[i][j],sizeof(double),1,fp);
	}
	MPI_Bcast(&a0[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&gamma[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDPDOMP::write_restart_settings(FILE *fp)
{
  fwrite(&temperature,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&seed,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDPDOMP::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&temperature,sizeof(double),1,fp);
    fread(&cut_global,sizeof(double),1,fp);
    fread(&seed,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&temperature,1,MPI_DOUBLE,0,world);
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

/* ---------------------------------------------------------------------- */

double PairDPDOMP::single(int i, int j, int itype, int jtype, double rsq,
		       double factor_coul, double factor_dpd, double &fforce)
{
  double r,rinv,wd,phi;

  if (rsq < EPSILON) {
    fforce = 0.0;
    return 0.0;
  }
  r = sqrt(rsq);

  rinv = 1.0/r;
  wd = 1.0 - r/cut[itype][jtype];
  fforce = a0[itype][jtype]*wd * factor_dpd*rinv;
  
  phi = -a0[itype][jtype] * r * (1.0 - 0.5*r/cut[itype][jtype]);
  return factor_dpd*phi;
}
