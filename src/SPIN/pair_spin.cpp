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

//Add. lib. for full neighb. list
#include "neigh_request.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairSpin::PairSpin(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  exch_flag = 0; 
  dmi_flag = 0;
  me_flag = 0;
}

/* ---------------------------------------------------------------------- */

PairSpin::~PairSpin()
{
  if (allocated) {
    memory->destroy(setflag);
    
    memory->destroy(cut_spin_exchange);
    memory->destroy(J_1);
    memory->destroy(J_2);
    memory->destroy(J_2);  
    
    memory->destroy(cut_spin_dmi);
    memory->destroy(DM);
    memory->destroy(v_dmx);
    memory->destroy(v_dmy);
    memory->destroy(v_dmz);

    memory->destroy(cut_spin_me);
    memory->destroy(ME);
    memory->destroy(v_mex);
    memory->destroy(v_mey);
    memory->destroy(v_mez);
    
    memory->destroy(spi);
    memory->destroy(spj);
    memory->destroy(fmi);
    memory->destroy(fmj);

    memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairSpin::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;  
  double evdwl,ecoul;
  double xtmp,ytmp,ztmp;
  double fmix,fmiy,fmiz,fmjx,fmjy,fmjz;
  double cut_ex_2,cut_dmi_2,cut_me_2;
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
  // Loop over neighbors of my itoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i]; 
    spi[0] = sp[i][0]; 
    spi[1] = sp[i][1]; 
    spi[2] = sp[i][2]; 
   
    //Loop on Neighbors
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      spj[0] = sp[j][0]; 
      spj[1] = sp[j][1]; 
      spj[2] = sp[j][2]; 

      fmi[0] = fmi[1] = fmi[2] = 0.0;
      fmj[0] = fmj[1] = fmj[2] = 0.0;
     
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;  //square or inter-atomic distance
      itype = type[i];
      jtype = type[j];

      //Exchange interaction
      if (exch_flag) {
        cut_ex_2 = cut_spin_exchange[itype][jtype]*cut_spin_exchange[itype][jtype];
        if (rsq <= cut_ex_2) {
          compute_exchange(i,j,rsq,fmi,fmj,spi,spj);   
        }
      }
      //DM interaction
      if (dmi_flag){
        cut_dmi_2 = cut_spin_dmi[itype][jtype]*cut_spin_dmi[itype][jtype];
        if (rsq <= cut_dmi_2){
          compute_dmi(i,j,fmi,fmj,spi,spj);
        } 
      }
      //ME interaction
      if (me_flag){
        cut_me_2 = cut_spin_me[itype][jtype]*cut_spin_me[itype][jtype];
        if (rsq <= cut_me_2){
          compute_me(i,j,fmi,fmj,spi,spj);
        } 
      }

      fm[i][0] += fmi[0];	 
      fm[i][1] += fmi[1];	  	  
      fm[i][2] += fmi[2];

      if (newton_pair || j < nlocal) {
          fm[j][0] += fmj[0];	 
          fm[j][1] += fmj[1];	  	  
          fm[j][2] += fmj[2];
          }
      }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */
void PairSpin::compute_exchange(int i, int j, double rsq, double *fmi,  double *fmj, double *spi, double *spj)
{
  int *type = atom->type;  
  int itype, jtype;
  double dmix,dmiy,dmiz;	
  double Jex, ra;
  itype = type[i];
  jtype = type[j];
          
  ra = rsq/J_3[itype][jtype]/J_3[itype][jtype]; 
  Jex = 4.0*J_1[itype][jtype]*ra;
  Jex *= (1.0-J_2[itype][jtype]*ra);
  Jex *= exp(-ra);

  fmi[0] += Jex*spj[0];
  fmi[1] += Jex*spj[1];
  fmi[2] += Jex*spj[2];
          
  fmj[0] += Jex*spi[0];
  fmj[1] += Jex*spi[1];
  fmj[2] += Jex*spi[2];
  
}

/* ---------------------------------------------------------------------- */
void PairSpin::compute_dmi(int i, int j, double *fmi,  double *fmj, double *spi, double *spj)
{
  int *type = atom->type;  
  int itype, jtype;
  double dmix,dmiy,dmiz;	
  itype = type[i];
  jtype = type[j];
          
  dmix = DM[itype][jtype]*v_dmx[itype][jtype];
  dmiy = DM[itype][jtype]*v_dmy[itype][jtype];
  dmiz = DM[itype][jtype]*v_dmz[itype][jtype];

  fmi[0] += spj[1]*dmiz-spj[2]*dmiy;
  fmi[1] += spj[2]*dmix-spj[0]*dmiz;
  fmi[2] += spj[0]*dmiy-spj[1]*dmix;

  fmj[0] -= spi[1]*dmiz-spi[2]*dmiy;
  fmj[1] -= spi[2]*dmix-spi[0]*dmiz;
  fmj[2] -= spi[0]*dmiy-spi[1]*dmix;
}

/* ---------------------------------------------------------------------- */
void PairSpin::compute_me(int i, int j, double *fmi,  double *fmj, double *spi, double *spj)
{
  int *type = atom->type;  
  int itype, jtype;
  itype = type[i];
  jtype = type[j];
  double **sp = atom->sp;
  double **x = atom->x;
  double meix,meiy,meiz;
  double rx, ry, rz, inorm;

  rx = x[j][0] - x[i][0];
  ry = x[j][1] - x[i][1];
  rz = x[j][2] - x[i][2];
  inorm = 1.0/sqrt(rx*rx+ry*ry+rz*rz);
  rx *= inorm;
  ry *= inorm; 
  rz *= inorm; 

  meix = v_mey[itype][jtype]*rz - v_mez[itype][jtype]*ry; 
  meiy = v_mez[itype][jtype]*rx - v_mex[itype][jtype]*rz; 
  meiz = v_mex[itype][jtype]*ry - v_mey[itype][jtype]*rx; 

  meix *= ME[itype][jtype]; 
  meiy *= ME[itype][jtype]; 
  meiz *= ME[itype][jtype]; 

  fmi[0] += spj[1]*meiz - spj[2]*meiy;
  fmi[1] += spj[2]*meix - spj[0]*meiz;
  fmi[2] += spj[0]*meiy - spj[1]*meix;
          
  fmj[0] -= spi[1]*meiz - spi[2]*meiy;
  fmj[1] -= spi[2]*meix - spi[0]*meiz;
  fmj[2] -= spi[0]*meiy - spi[1]*meix;

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
  memory->create(J_1,n+1,n+1,"pair:J_1");
  memory->create(J_2,n+1,n+1,"pair:J_2");  
  memory->create(J_3,n+1,n+1,"pair:J_3");

  memory->create(cut_spin_dmi,n+1,n+1,"pair:cut_spin_dmi");
  memory->create(DM,n+1,n+1,"pair:DM");
  memory->create(v_dmx,n+1,n+1,"pair:DM_vector_x");
  memory->create(v_dmy,n+1,n+1,"pair:DM_vector_y");
  memory->create(v_dmz,n+1,n+1,"pair:DM_vector_z");
 
  memory->create(cut_spin_me,n+1,n+1,"pair:cut_spin_me");
  memory->create(ME,n+1,n+1,"pair:ME");
  memory->create(v_mex,n+1,n+1,"pair:ME_vector_x");
  memory->create(v_mey,n+1,n+1,"pair:ME_vector_y");
  memory->create(v_mez,n+1,n+1,"pair:ME_vector_z");
 
  memory->create(spi,3,"pair:spi");
  memory->create(spj,3,"pair:spj");
  memory->create(fmi,3,"pair:fmi");
  memory->create(fmj,3,"pair:fmj");
 
  memory->create(cutsq,n+1,n+1,"pair:cutsq");  
  
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSpin::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2)
    error->all(FLERR,"Incorrect number of args in pair_style pair/spin command");

  if (strcmp(update->unit_style,"metal") != 0)
    error->all(FLERR,"Spin simulations require metal unit style");
    
  cut_spin_pair_global = force->numeric(FLERR,arg[0]);
    
  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_spin_exchange[i][j] = cut_spin_pair_global;
          cut_spin_dmi[i][j] = cut_spin_pair_global;
          cut_spin_me[i][j] = cut_spin_pair_global;
        }
  }
   
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type spin pairs (only one for now)
------------------------------------------------------------------------- */

void PairSpin::coeff(int narg, char **arg)
{
  const double hbar = force->hplanck/MY_2PI;

  if (!allocated) allocate();

  if (strcmp(arg[2],"exchange")==0){
    if (narg != 7) error->all(FLERR,"Incorrect args in pair_style command");
    exch_flag = 1;    
    
    int ilo,ihi,jlo,jhi;
    force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
    force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);
    
    const double rij = force->numeric(FLERR,arg[3]);
    const double J1 = (force->numeric(FLERR,arg[4]))/hbar;
    const double J2 = force->numeric(FLERR,arg[5]);  
    const double J3 = force->numeric(FLERR,arg[6]); 
  
    int count = 0;
    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo,i); j <= jhi; j++) {
        cut_spin_exchange[i][j] = rij;   
        J_1[i][j] = J1;
        J_2[i][j] = J2;
        J_3[i][j] = J3;
        setflag[i][j] = 1;
        count++;
      }
    }
    if (count == 0) error->all(FLERR,"Incorrect args in pair_style command");  
  } else if (strcmp(arg[2],"dmi")==0) {
    if (narg != 8) error->all(FLERR,"Incorrect args in pair_style command");
    dmi_flag = 1;   

    int ilo,ihi,jlo,jhi;
    force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
    force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);
    
    const double rij = force->numeric(FLERR,arg[3]);
    const double dm = (force->numeric(FLERR,arg[4]))/hbar;
    double dmx = force->numeric(FLERR,arg[5]);  
    double dmy = force->numeric(FLERR,arg[6]); 
    double dmz = force->numeric(FLERR,arg[7]); 

    double inorm = 1.0/(dmx*dmx+dmy*dmy+dmz*dmz);
    dmx *= inorm; 
    dmy *= inorm; 
    dmz *= inorm; 
 
    int count = 0;
    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo,i); j <= jhi; j++) {
        cut_spin_dmi[i][j] = rij;
        DM[i][j] = dm;
        v_dmx[i][j] = dmx;
        v_dmy[i][j] = dmy;
        v_dmz[i][j] = dmz;
        setflag[i][j] = 1;
        count++;
      }
    }
    if (count == 0) error->all(FLERR,"Incorrect args in pair_style command"); 
  } else if (strcmp(arg[2],"me")==0) {
    if (narg != 8) error->all(FLERR,"Incorrect args in pair_style command");
    me_flag = 1;    
    int ilo,ihi,jlo,jhi;
    force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
    force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);
    
    const double rij = force->numeric(FLERR,arg[3]);
    const double me = (force->numeric(FLERR,arg[4]))/hbar;
    double mex = force->numeric(FLERR,arg[5]);  
    double mey = force->numeric(FLERR,arg[6]); 
    double mez = force->numeric(FLERR,arg[7]); 

    double inorm = 1.0/(mex*mex+mey*mey+mez*mez);
    mex *= inorm; 
    mey *= inorm; 
    mez *= inorm; 
 
    int count = 0;
    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo,i); j <= jhi; j++) {
        cut_spin_me[i][j] = rij;
        DM[i][j] = me;
        v_mex[i][j] = mex;
        v_mey[i][j] = mey;
        v_mez[i][j] = mez;
        setflag[i][j] = 1;
        count++;
      }
    }
    if (count == 0) error->all(FLERR,"Incorrect args in pair_style command"); 
  } else error->all(FLERR,"Incorrect args in pair_style command");

  // check if Jex [][] still works for Ferrimagnetic exchange  

}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSpin::init_style()
{
  if (!atom->sp_flag || !atom->mumag_flag)
    error->all(FLERR,"Pair spin requires atom attributes sp, mumag");

  neighbor->request(this,instance_me);

#define FULLNEI
#if defined FULLNEI
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
#endif
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSpin::init_one(int i, int j)
{
   
   if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cut_spin_pair_global;
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
        if (exch_flag){
          fwrite(&J_1[i][j],sizeof(double),1,fp);
          fwrite(&J_2[i][j],sizeof(double),1,fp);
          fwrite(&J_3[i][j],sizeof(double),1,fp);
          fwrite(&cut_spin_exchange[i][j],sizeof(double),1,fp);
        }
        if (dmi_flag) {
          fwrite(&DM[i][j],sizeof(double),1,fp);
          fwrite(&v_dmx[i][j],sizeof(double),1,fp);
          fwrite(&v_dmy[i][j],sizeof(double),1,fp);
          fwrite(&v_dmz[i][j],sizeof(double),1,fp);
          fwrite(&cut_spin_dmi[i][j],sizeof(double),1,fp);
        } 
        if (me_flag) {
          fwrite(&ME[i][j],sizeof(double),1,fp);
          fwrite(&v_mex[i][j],sizeof(double),1,fp);
          fwrite(&v_mey[i][j],sizeof(double),1,fp);
          fwrite(&v_mez[i][j],sizeof(double),1,fp);
          fwrite(&cut_spin_me[i][j],sizeof(double),1,fp);
        } 
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
  fwrite(&cut_spin_pair_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSpin::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_spin_pair_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_spin_pair_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world); 
}
