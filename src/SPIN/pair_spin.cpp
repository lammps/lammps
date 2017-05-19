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

#include <math.h>
#include <stdlib.h>
#include "pair_spin.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "math_const.h"
#include "error.h"
#include "update.h"
#include <string.h>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairSpin::PairSpin(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairSpin::~PairSpin()
{
  if (allocated) {
    memory->destroy(setflag);
    
    memory->destroy(cut_spin_exchange);
    memory->destroy(cut_spin_dipolar);
    memory->destroy(J_1);
    memory->destroy(J_2);
    memory->destroy(J_2);
    memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairSpin::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;  
  double evdwl,ecoul;
  double xtmp,ytmp,ztmp,fmix,fmiy,fmiz,fmjx,fmjy,fmjz,omx,omy,omz;
  double cut, Jex, ra;
  double rsq,rd,delx,dely,delz;
  int *ilist,*jlist,*numneigh,**firstneigh;  

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;
  
  double **x = atom->x;
  double **fm = atom->fm;
  double *mumag = atom->mumag;
  double **sp = atom->sp;	
  int *type = atom->type;  
  int nlocal = atom->nlocal;  
  int newton_pair = force->newton_pair;
  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // Pair spin computations
  // Loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];  
   
    //printf(":::::::Loop for atom numb. %d,  jnum=%d ::::::: \n",i,jnum);
     
    //printf("Test print real atom: i=%d, sx=%g, sy=%g, sz=%g \n",i,sp[i][0],sp[i][1],sp[i][2]);
    //printf("Test print real atom: i=%d, rx=%g, ry=%g, rz=%g \n",i,x[i][0],x[i][1],x[i][2]);

    int testcount=0;     

    //Exchange interaction
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      fmix = fmiy = fmiz = 0.0;
      fmjx = fmjy = fmjz = 0.0;
     
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;  //square or inter-atomic distance
      rd = sqrt(rsq); //Inter-atomic distance
      cut = cut_spin_exchange_global;
      
      if (rd <= cut) {
          itype = type[i];
          jtype = type[j];

          ra = (rd/J_3[itype][jtype])*(rd/J_3[itype][jtype]); 
          Jex = 4.0*J_1[itype][jtype]*ra;
          Jex *= (1.0-J_2[itype][jtype]*ra);
          Jex *= exp(-ra);
          Jex *= mumag[ii]*mumag[jj];
      
          fmix = Jex*sp[j][0];
          fmiy = Jex*sp[j][1];
          fmiz = Jex*sp[j][2];
          
          fmjx = Jex*sp[i][0];
          fmjy = Jex*sp[i][1];
          fmjz = Jex*sp[i][2];

 
          
          //printf("Neighb pair: i=%d, j=%d \n",i,j);
          //printf("Test print ghost/real neib atom: i=%d, j=%d, sx=%g, sy=%g, sz=%g \n",i,j,sp[j][0],sp[j][1],sp[j][2]);
          //printf("Test g/r neib pair: i=%d, j=%d, rx=%g, ry=%g, rz=%g \n",i,j,x[j][0],x[j][1],x[j][2]);
          //printf("Atom i: %d of type %d, Atom j: %d of type %d, \n",i,itype,j,jtype);
          //printf("Exchange pair (%d,%d), Jij=%g, rij=%g \n",i,j,Jex,rd);
          testcount++;
	  }


      //printf("Test print ghost/real atom: j=%d, sx=%g, sy=%g, sz=%g \n",j,sp[j][0],sp[j][1],sp[j][2]);
      //printf("Test print ghost/real atom: j=%d, rx=%g, ry=%g, rz=%g \n",j,x[j][0],x[j][1],x[j][2]);
     
      fm[i][0] += fmix;	 
      fm[i][1] += fmiy;	  	  
      fm[i][2] += fmiz;
      
      if (newton_pair || j < nlocal) { 
      fm[j][0] += fmjx;	 
      fm[j][1] += fmjy;	  	  
      fm[j][2] += fmjz;
      }

      //printf("Val fm %d: [%g,%g,%g] \n",i,fm[i][0],fm[i][1],fm[i][2]);
      //printf("Val fm %d: [%g,%g,%g] \n",j,fm[j][0],fm[j][1],fm[j][2]);

      }


     // printf("Test count %d \n",testcount);
     
  }
  //printf("New pair val: %d \n",newton_pair);
  //printf("vals exchange: Jx=%g, Jy=%g, Jz=%g \n",fm[0][0],fm[0][1],fm[0][2]);  
  //printf("test exchange. 1;i=0, fx=%g, fy=%g, fz=%g \n",fm[0][0],fm[0][1],fm[0][2]);
  //printf("::::::::::::::::::::::::: End loop ::::::::::::::::::::::::\n");
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSpin::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;
      
  memory->create(cut_spin_exchange,n+1,n+1,"pair:cut_spin_exchange");
  memory->create(cut_spin_dipolar,n+1,n+1,"pair:cut_spin_dipolar");
  memory->create(J_1,n+1,n+1,"pair:J_1");
  memory->create(J_2,n+1,n+1,"pair:J_2");  
  memory->create(J_3,n+1,n+1,"pair:J_3");
  
  memory->create(cutsq,n+1,n+1,"pair:cutsq");  
  
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSpin::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2)
    error->all(FLERR,"Incorrect number of args in pair_style pair/spin command");

  if (strcmp(update->unit_style,"electron") == 0)
    error->all(FLERR,"Cannot (yet) use 'electron' units with spins");
    
  cut_spin_exchange_global = force->numeric(FLERR,arg[0]);
    
  if (narg == 1) cut_spin_dipolar_global = cut_spin_exchange_global;
  else cut_spin_dipolar_global = force->numeric(FLERR,arg[1]);
  
  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_spin_exchange[i][j] = cut_spin_exchange_global;
          cut_spin_dipolar[i][j] = cut_spin_dipolar_global;
        }
  }
   
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type spin pairs (only one for now)
------------------------------------------------------------------------- */

void PairSpin::coeff(int narg, char **arg)
{

  if (narg != 5)
    error->all(FLERR,"Incorrect number of args for pair spin coefficients");
  if (!allocated) allocate();	
  
  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);
    
  double J1 = force->numeric(FLERR,arg[2]);
  double J2 = force->numeric(FLERR,arg[3]);  
  double J3 = force->numeric(FLERR,arg[4]); 
  
  double hbar = force->hplanck/MY_2PI;
  J1 /= hbar;
    
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      J_1[i][j] = J1;
      J_2[i][j] = J2;
      J_3[i][j] = J3;
      setflag[i][j] = 1;
      count++;
    }
  }
  if (count == 0) error->all(FLERR,"Incorrect args for spinpair coefficients");  
 
  //Simple (Anti)Ferromagnetic exchange for now. 
  //Check if Jex [][] still works for Ferrimagnetic exchange
  
}



/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSpin::init_style()
{
  if (!atom->sp_flag || !atom->mumag_flag)
    error->all(FLERR,"Pair spin requires atom attributes sp, mumag");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSpin::init_one(int i, int j)
{
   
   if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cut_spin_exchange_global;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSpin::write_restart(FILE *fp)
{
  write_restart_settings(fp);
  
  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&J_1[i][j],sizeof(double),1,fp);
        fwrite(&J_2[i][j],sizeof(double),1,fp);
        fwrite(&J_3[i][j],sizeof(double),1,fp);
        fwrite(&cut_spin_exchange[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSpin::read_restart(FILE *fp)
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
          fread(&J_1[i][j],sizeof(double),1,fp);
          fread(&J_2[i][j],sizeof(double),1,fp);
          fread(&J_2[i][j],sizeof(double),1,fp);
          fread(&cut_spin_exchange[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&J_1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&J_2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&J_3[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_spin_exchange[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}
 

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSpin::write_restart_settings(FILE *fp)
{
  fwrite(&cut_spin_exchange_global,sizeof(double),1,fp);
  fwrite(&cut_spin_dipolar_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSpin::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_spin_exchange_global,sizeof(double),1,fp);
    fread(&cut_spin_dipolar_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_spin_exchange_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_spin_dipolar_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world); 
}
