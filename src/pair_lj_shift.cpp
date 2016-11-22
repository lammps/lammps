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
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_lj_shift.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

//#define MIN(a,b) ((a) < (b) ? (a) : (b))
//#define MAX(a,b) ((a) > (b) ? (a) : (b))


/* ---------------------------------------------------------------------- */

PairLJShift::PairLJShift(LAMMPS *lmp) : Pair(lmp) {

respa_enable = 1;
writedata = 1;


}

/* ---------------------------------------------------------------------- */

PairLJShift::~PairLJShift()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(k0);
    memory->destroy(shift);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairLJShift::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,forcelj,factor_lj,philj,rinv,r14inv,tangh;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0; 
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  //printf ("evflag=%i,vflag_fdotr=%i\n",evflag,vflag_fdotr);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

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
        rinv = 1.0/sqrt(rsq);
	r2inv = 1.0/rsq;
	r14inv = r2inv*r2inv*r2inv*r2inv*r2inv*r2inv*r2inv;
	tangh=tanh(k0[itype][jtype]*(sqrt(rsq)-shift[itype][jtype]));
	forcelj = lj1[itype][jtype]*r14inv + lj2[itype][jtype]*sqrt(rsq)*(1.0 - tangh*tangh);
	fpair = factor_lj*forcelj*r2inv;

	f[i][0] += delx*fpair;
	f[i][1] += dely*fpair;
	f[i][2] += delz*fpair;
	if (newton_pair || j < nlocal) {
	  f[j][0] -= delx*fpair;
	  f[j][1] -= dely*fpair;
	  f[j][2] -= delz*fpair;
	}

	if(eflag){
//	tangh=tanh(k0[itype][jtype]*(sqrt(rsq)-shift[itype][jtype]));
	philj = sigma[itype][jtype]*r14inv + (0.5)*epsilon[itype][jtype]*(1.0-tangh);
        evdwl=philj-offset[itype][jtype];
	}
/*	if (vflag == 1) {
          if (newton_pair == 0 && j >= nlocal) fpair *= 0.5;
          virial[0] += delx*delx*fpair;
          virial[1] += dely*dely*fpair;
          virial[2] += delz*delz*fpair;
          virial[3] += delx*dely*fpair;
          virial[4] += delx*delz*fpair;
          virial[5] += dely*delz*fpair;
        }
*/
	if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
	

      }
    }
  }
	if (vflag_fdotr)  virial_fdotr_compute();
//	virial_compute();
 // printf("Eng_vdwl=%f\n",evdwl);
}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJShift::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

//  setflag = memory->create(n+1,n+1,"pair:setflag");
  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;
/*
  cutsq = memory->create(n+1,n+1,"pair:cutsq");

  cut = memory->create(n+1,n+1,"pair:cut");
  epsilon = memory->create(n+1,n+1,"pair:epsilon");
  sigma = memory->create(n+1,n+1,"pair:sigma");
  k0 = memory->create(n+1,n+1,"pair:k0");
  shift = memory->create(n+1,n+1,"pair:shift");
  lj1 = memory->create(n+1,n+1,"pair:lj1");
  lj2 = memory->create(n+1,n+1,"pair:lj2");
  offset = memory->create(n+1,n+1,"pair:offset");
*/


  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(k0,n+1,n+1,"pair:k0");
  memory->create(shift,n+1,n+1,"pair:shift");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(offset,n+1,n+1,"pair:offset");

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJShift::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

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

void PairLJShift::coeff(int narg, char **arg)
{
  //printf("narg=%i \n",narg);
  if (narg < 5 || narg > 6) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();


  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);
  
  double sigma_one = force->numeric(FLERR,arg[2]);
  double epsilon_one =  force->numeric(FLERR,arg[3]);
  double k0_one = force->numeric(FLERR,arg[4]);
  double shift_one = force->numeric(FLERR,arg[5]);

  // double epsilon_one = force->numeric(arg[2]);
 // double sigma_one = force->numeric(arg[3]);

  printf("ilo=%i,ihi=%i,jlo=%i,jhi=%i\n",ilo,ihi,jlo,jhi);
  printf("sigma_one=%f,epsilon_one=%f,k0_one=%f,shift_one=%f\n",sigma_one,epsilon_one,k0_one,shift_one);
  //lamda = atof(arg[6]);
  double cut_one = cut_global;
  //if (narg == 5) cut_one = atof(arg[4]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      k0[i][j] = k0_one;
      shift[i][j] = shift_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJShift::init_style()
{
  // request regular or rRESPA neighbor lists

  int irequest;

  if (update->whichflag == 0 && strcmp(update->integrate_style,"respa") == 0) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;

    if (respa == 0) irequest = neighbor->request(this);
    else if (respa == 1) {
      irequest = neighbor->request(this);
      neighbor->requests[irequest]->id = 1;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respainner = 1;
      irequest = neighbor->request(this);
      neighbor->requests[irequest]->id = 3;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respaouter = 1;
    } else {
      irequest = neighbor->request(this);
      neighbor->requests[irequest]->id = 1;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respainner = 1;
      irequest = neighbor->request(this);
      neighbor->requests[irequest]->id = 2;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respamiddle = 1;
      irequest = neighbor->request(this);
      neighbor->requests[irequest]->id = 3;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respaouter = 1;
    }

  } else irequest = neighbor->request(this);
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   regular or rRESPA
------------------------------------------------------------------------- */

void PairLJShift::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  else if (id == 1) listinner = ptr;
  else if (id == 2) listmiddle = ptr;
  else if (id == 3) listouter = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJShift::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
			       sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  lj1[i][j] = 14.0 * pow(sigma[i][j],14.0);
  lj2[i][j] = 0.5 *  epsilon[i][j] * k0[i][j];

  if (offset_flag) {
     offset[i][j] = 0;
  } else offset[i][j] = 0.0;

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];

  
  epsilon[j][i] = epsilon[i][j];
  sigma[j][i] = sigma[i][j];
  k0[j][i] = k0[i][j];
  shift[j][i] = shift[i][j];
  cut[j][i] = cut[i][j];
  
  
  offset[j][i] = offset[i][j];

  // set & error check interior rRESPA cutoff

  if (strcmp(update->integrate_style,"respa") == 0) {
    if (((Respa *) update->integrate)->level_inner >= 0) {
      cut_respa = ((Respa *) update->integrate)->cutoff;
      if (cut[i][j] < cut_respa[3])
	error->all(FLERR,"Pair cutoff < Respa interior cutoff");
    }
  } else cut_respa = NULL;

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);

    double PI = 4.0*atan(1.0);
    double sig2 = sigma[i][j]*sigma[i][j];
    double sig6 = sig2*sig2*sig2;
    double rc3 = cut[i][j]*cut[i][j]*cut[i][j];
    double rc6 = rc3*rc3;
    double rc9 = rc3*rc6;
    etail_ij = 0;
    ptail_ij = 0;
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJShift::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
	fwrite(&epsilon[i][j],sizeof(double),1,fp);
	fwrite(&sigma[i][j],sizeof(double),1,fp);
	fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJShift::read_restart(FILE *fp)
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
	  fread(&epsilon[i][j],sizeof(double),1,fp);
	  fread(&sigma[i][j],sizeof(double),1,fp);
	  fread(&cut[i][j],sizeof(double),1,fp);
	}
	MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJShift::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJShift::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairLJShift::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,  double &fpair)

{
  double r2inv,r6inv,forcelj,philj,r14inv,tangh,rinv;

	rinv = 1.0/sqrt(rsq);
	r2inv = 1.0/rsq;
	r14inv = r2inv*r2inv*r2inv*r2inv*r2inv*r2inv*r2inv;
	tangh=tanh(k0[itype][jtype]*(sqrt(rsq)-shift[itype][jtype]));
	forcelj = lj1[itype][jtype]*r14inv + lj2[itype][jtype]*sqrt(rsq)*(1.0 - tangh*tangh);
	fpair = factor_lj*forcelj*r2inv;


/*
tangh=tanh(sqrt(rsq)*10.0-13.5);
forcelj = 14.0*r14inv + 5.0*sqrt(rsq)*(1.0 - tangh*tangh);
fpair = factor_lj*forcelj*r2inv;
*/

  tangh=tanh(k0[itype][jtype]*(sqrt(rsq)-shift[itype][jtype]));
  philj = sigma[itype][jtype]*r14inv + (0.5)*epsilon[itype][jtype]*(1-tangh) - offset[itype][jtype];
/*  
  tangh=tanh(10*(sqrt(rsq)-1.35));
philj = r14inv + (0.5)*(1-tangh) - offset[itype][jtype];
*/
  return factor_lj*philj;


}
