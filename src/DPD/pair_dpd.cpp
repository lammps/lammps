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
#include "pair_dpd.h"
#include "atom.h"
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

#define EPSILON 1.0e-10

/* ---------------------------------------------------------------------- */

PairDPD::PairDPD(LAMMPS *lmp) : Pair(lmp)
{
  random = NULL;
}

/* ---------------------------------------------------------------------- */

PairDPD::~PairDPD()
{
  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);

    memory->destroy_2d_double_array(cut);
    memory->destroy_2d_double_array(a0);
    memory->destroy_2d_double_array(gamma);
    memory->destroy_2d_double_array(sigma);
  }

  if (random) delete random;
}

/* ---------------------------------------------------------------------- */

void PairDPD::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  double rsq,r,rinv,dot,wd,randnum,fforce,factor_dpd,phi;
  int *ilist,*jlist,*numneigh,**firstneigh;

  eng_vdwl = 0.0;
  if (vflag) for (i = 0; i < 6; i++) virial[i] = 0.0;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double dtinvsqrt = 1.0/sqrt(update->dt);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
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

      if (rsq < cutsq[itype][jtype]) {
	r = sqrt(rsq);
	if (r < EPSILON) continue;     // r can be 0.0 in DPD systems
	rinv = 1.0/r;
	delvx = vxtmp - v[j][0];
	delvy = vytmp - v[j][1];
	delvz = vztmp - v[j][2];
	dot = delx*delvx + dely*delvy + delz*delvz;
	wd = 1.0 - r/cut[itype][jtype];
	randnum = random->gaussian();

	// conservative force = a0 * wd
	// drag force = -gamma * wd^2 * (delx dot delv) / r
	// random force = sigma * wd * rnd * dtinvsqrt;

	fforce = a0[itype][jtype]*wd;
	fforce -= gamma[itype][jtype]*wd*wd*dot*rinv;
	fforce += sigma[itype][jtype]*wd*randnum*dtinvsqrt;
	fforce *= factor_dpd*rinv;	

	f[i][0] += delx*fforce;
	f[i][1] += dely*fforce;
	f[i][2] += delz*fforce;
	if (newton_pair || j < nlocal) {
	  f[j][0] -= delx*fforce;
	  f[j][1] -= dely*fforce;
	  f[j][2] -= delz*fforce;
	}

	if (eflag) {
	  phi = a0[itype][jtype] * r * (1.0 - 0.5*r/cut[itype][jtype]);
	  if (newton_pair || j < nlocal) eng_vdwl += factor_dpd*phi;
	  else eng_vdwl += 0.5*factor_dpd*phi;
	}

	if (vflag == 1) {
	  if (newton_pair == 0 && j >= nlocal) fforce *= 0.5;
	  virial[0] += delx*delx*fforce;
	  virial[1] += dely*dely*fforce;
	  virial[2] += delz*delz*fforce;
	  virial[3] += delx*dely*fforce;
	  virial[4] += delx*delz*fforce;
	  virial[5] += dely*delz*fforce;
	}
      }
    }
  }
  if (vflag == 2) virial_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays 
------------------------------------------------------------------------- */

void PairDPD::allocate()
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

void PairDPD::settings(int narg, char **arg)
{
  if (narg != 3) error->all("Illegal pair_style command");

  temperature = atof(arg[0]);
  cut_global = atof(arg[1]);
  seed = atoi(arg[2]);

  // initialize Marsaglia RNG with processor-unique seed

  if (seed <= 0 || seed > 900000000)
    error->all("Illegal fix pair_style command");
  if (random) delete random;
  random = new RanMars(lmp,seed + comm->me);

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

void PairDPD::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5) error->all("Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double a0_one = atof(arg[2]);
  double gamma_one = atof(arg[3]);

  double cut_one = cut_global;
  if (narg == 5) cut_one = atof(arg[4]);

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

void PairDPD::init_style()
{
  // check that atom style is dpd or hybrid with dpd
  // else compute() will not have ghost atom velocities

  if (atom->style_match("dpd") == 0)
    error->all("Pair style dpd requires atom style dpd");

  // if newton off, forces between atoms ij will be double computed
  // using different random numbers

  if (force->newton_pair == 0 && comm->me == 0) error->warning(
      "DPD potential needs newton pair on for momentum conservation");

  int irequest = neighbor->request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairDPD::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all("All pair coeffs are not set");

  sigma[i][j] = sqrt(2.0*temperature*gamma[i][j]);
     
  cut[j][i] = cut[i][j];
  a0[j][i] = a0[i][j];
  gamma[j][i] = gamma[i][j];
  sigma[j][i] = sigma[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDPD::write_restart(FILE *fp)
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

void PairDPD::read_restart(FILE *fp)
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

void PairDPD::write_restart_settings(FILE *fp)
{
  fwrite(&temperature,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&seed,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDPD::read_restart_settings(FILE *fp)
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

  if (random) delete random;
  random = new RanMars(lmp,seed + comm->me);
}

/* ---------------------------------------------------------------------- */

void PairDPD::single(int i, int j, int itype, int jtype, double rsq,
		     double factor_coul, double factor_dpd, int eflag,
		     One &one)
{
  double r,rinv,wd,phi;

  r = sqrt(rsq);
  if (r < EPSILON) {
    one.fforce = 0.0;
    if (eflag) one.eng_vdwl = one.eng_coul = 0.0;
    return;
  }

  rinv = 1.0/r;
  wd = 1.0 - r/cut[itype][jtype];

  one.fforce = a0[itype][jtype]*wd * factor_dpd*rinv;
  
  if (eflag) {
    phi = a0[itype][jtype] * r * (1.0 - 0.5*r/cut[itype][jtype]);
    one.eng_vdwl = factor_dpd*phi;
    one.eng_coul = 0.0;
  }
}
