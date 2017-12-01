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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "pair_hybrid.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "math_const.h"
#include "memory.h"
#include "pair_spin_soc_landau.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairSpinSocLandau::PairSpinSocLandau(LAMMPS *lmp) : Pair(lmp)
{
  hbar = force->hplanck/MY_2PI;

  newton_pair_spin = 0; // no newton pair for now => to be corrected
 // newton_pair = 0;

  single_enable = 0;
  soc_neel_flag = 0; 

  no_virial_fdotr_compute = 1;
}

/* ---------------------------------------------------------------------- */

PairSpinSocLandau::~PairSpinSocLandau()
{
  if (allocated) {
    memory->destroy(setflag);
    
    memory->destroy(cut_soc_landau);
    memory->destroy(K1);
    memory->destroy(K1_mech);
    memory->destroy(K2);
    memory->destroy(K3);  

    memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairSpinSocLandau::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;  
  double evdwl, ecoul;
  double xi[3], rij[3];
  double spi[3], spj[3];
  double fi[3], fj[3];
  double fmi[3], fmj[3];
  double cut_soc_landau_2, cut_soc_global2;
  double rsq, rd, inorm;
  int *ilist,*jlist,*numneigh,**firstneigh;  

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;
  cut_soc_global2 = cut_soc_global*cut_soc_global;
  
  double **x = atom->x;
  double **f = atom->f;
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

  // pair spin computations
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xi[0] = x[i][0];
    xi[1] = x[i][1];
    xi[2] = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i]; 
    spi[0] = sp[i][0]; 
    spi[1] = sp[i][1]; 
    spi[2] = sp[i][2];
  
    // loop on neighbors
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      spj[0] = sp[j][0]; 
      spj[1] = sp[j][1]; 
      spj[2] = sp[j][2]; 

      evdwl = 0.0;

      fi[0] = fi[1] = fi[2] = 0.0;
      fj[0] = fj[1] = fj[2] = 0.0;
      fmi[0] = fmi[1] = fmi[2] = 0.0;
      fmj[0] = fmj[1] = fmj[2] = 0.0;
      rij[0] = rij[1] = rij[2] = 0.0;
     
      rij[0] = x[j][0] - xi[0];
      rij[1] = x[j][1] - xi[1];
      rij[2] = x[j][2] - xi[2];
      rsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]; 
      inorm = 1.0/sqrt(rsq);
      rij[0] *= inorm;
      rij[1] *= inorm;
      rij[2] *= inorm;

      itype = type[i];
      jtype = type[j];

      // compute mag. and mech. components of soc
      cut_soc_landau_2 = cut_soc_landau[itype][jtype]*cut_soc_landau[itype][jtype];
      if (rsq <= cut_soc_landau_2) {
        compute_soc_landau(i,j,rsq,rij,fmi,fmj,spi,spj);   
        compute_soc_mech_landau(i,j,rsq,rij,fi,fj,spi,spj);
      }

      f[i][0] += fi[0];	 
      f[i][1] += fi[1];	  	  
      f[i][2] += fi[2];
      fm[i][0] += fmi[0];	 
      fm[i][1] += fmi[1];	  	  
      fm[i][2] += fmi[2];

//      if (newton_pair || j < nlocal) {  =>  to be corrected
      if (newton_pair_spin) {
	f[j][0] += fj[0];	 
        f[j][1] += fj[1];	  	  
        f[j][2] += fj[2];
        fm[j][0] += fmj[0];	 
        fm[j][1] += fmj[1];	  	  
        fm[j][2] += fmj[2];
      }
 
      if (eflag) {
	if (rsq <= cut_soc_landau_2) {
	  evdwl -= spi[0]*fmi[0];
	  evdwl -= spi[1]*fmi[1];
	  evdwl -= spi[2]*fmi[2];
	  evdwl *= hbar;
	} else evdwl = 0.0;
      }

      if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
	  evdwl,ecoul,fi[0],fi[1],fi[2],rij[0],rij[1],rij[2]);
    }
  }  

  if (vflag_fdotr) virial_fdotr_compute();
  
}

/* ---------------------------------------------------------------------- */

void PairSpinSocLandau::compute_soc_landau(int i, int j, double rsq, double rij[3], double fmi[3],  double fmj[3], double spi[3], double spj[3])
{
  int *type = atom->type;  
  int itype, jtype;
  double Kij, Kij_3, ra, scalar;
  itype = type[i];
  jtype = type[j];
          
  ra = rsq/K3[itype][jtype]/K3[itype][jtype]; 
  Kij = 4.0*K1[itype][jtype]*ra;
  Kij *= (1.0-K2[itype][jtype]*ra);
  Kij *= exp(-ra);

  scalar = rij[0]*spj[0]+rij[1]*spj[1]+rij[2]*spj[2];
  Kij_3 = Kij/3.0;

  fmi[0] += Kij*scalar*rij[0]-Kij_3*spj[0];
  fmi[1] += Kij*scalar*rij[1]-Kij_3*spj[1];
  fmi[2] += Kij*scalar*rij[2]-Kij_3*spj[2];
          
  fmj[0] -= Kij*scalar*rij[0]+Kij_3*spi[0];
  fmj[1] -= Kij*scalar*rij[1]+Kij_3*spi[1];
  fmj[2] -= Kij*scalar*rij[2]+Kij_3*spi[2];

}

/* ---------------------------------------------------------------------- */

void PairSpinSocLandau::compute_soc_mech_landau(int i, int j, double rsq, double rij[3], double fi[3],  double fj[3], double spi[3], double spj[3])
{
  int *type = atom->type;  
  int itype, jtype;
  double scalar_si_sj, scalar_rij_si, scalar_rij_sj;
  double K_mech, Kij, dKij, ra, rr, drij, iK3;
  double t1, t2, t3;
  itype = type[i];
  jtype = type[j];

  scalar_si_sj = spi[0]*spj[0]+spi[1]*spj[1]+spi[2]*spj[2];
  scalar_rij_si = rij[0]*spi[0]+rij[1]*spi[1]+rij[2]*spi[2];
  scalar_rij_sj = rij[0]*spj[0]+rij[1]*spj[1]+rij[2]*spj[2];

  K_mech = K1_mech[itype][jtype];        
  iK3 = 1.0/(K3[itype][jtype]*K3[itype][jtype]);

  drij = sqrt(rsq);
  ra = rsq*iK3; 
  rr = drij*iK3;

  Kij *= (1.0-K2[itype][jtype]*ra);
  Kij *= 4.0*K_mech*ra*exp(-ra);

  dKij = 1.0-ra-K2[itype][jtype]*ra*(2.0-ra);
  dKij *= 8.0*K_mech*rr*exp(-ra);

  t1 = (dKij-2.0*Kij/drij)*scalar_rij_si*scalar_rij_sj;
  t1 -= scalar_si_sj*dKij/3.0;
  t2 = scalar_rij_sj*Kij/drij;
  t3 = scalar_rij_si*Kij/drij;

  fi[0] += t1*rij[0]+t2*spi[0]+t3*spj[0];
  fi[1] += t1*rij[1]+t2*spi[1]+t3*spj[1];
  fi[2] += t1*rij[2]+t2*spi[2]+t3*spj[2];
          
  fj[0] -= t1*rij[0]-t2*spi[0]-t3*spj[0];
  fj[1] -= t1*rij[1]-t2*spi[1]-t3*spj[1];
  fj[2] -= t1*rij[2]-t2*spi[2]-t3*spj[2];

}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSpinSocLandau::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;
      
  memory->create(cut_soc_landau,n+1,n+1,"pair/spin/soc/landau:cut_soc_landau");
  memory->create(K1,n+1,n+1,"pair/spin/soc/landau:K1");
  memory->create(K1_mech,n+1,n+1,"pair/spin/soc/landau:K1_mech");
  memory->create(K2,n+1,n+1,"pair/spin/soc/landau:K2");  
  memory->create(K3,n+1,n+1,"pair/spin/soc/landau:K3");
 
  memory->create(cutsq,n+1,n+1,"pair/spin/soc/landau:cutsq");  
  
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSpinSocLandau::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2)
    error->all(FLERR,"Incorrect number of args in pair_style pair/spin command");

  if (strcmp(update->unit_style,"metal") != 0)
    error->all(FLERR,"Spin simulations require metal unit style");
    
  cut_soc_global = force->numeric(FLERR,arg[0]);
    
  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_soc_landau[i][j] = cut_soc_global;
        }
  }
   
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type spin pairs (only one for now)
------------------------------------------------------------------------- */

void PairSpinSocLandau::coeff(int narg, char **arg)
{
  const double hbar = force->hplanck/MY_2PI;

  if (!allocated) allocate();
  
  // set mech_flag to 1 if magneto-mech simulation
  //no longer correct: can be hybrid without magneto-mech
  if (strstr(force->pair_style,"pair/spin")) {
    mech_flag = 0;
  } else if (strstr(force->pair_style,"hybrid/overlay")) {
    mech_flag = 1;
  } else error->all(FLERR,"Incorrect args in pair_style command");


  if (strcmp(arg[2],"neel")==0){
    if (narg != 7) error->all(FLERR,"Incorrect args in pair_style command");
    soc_neel_flag = 1;    
    
    int ilo,ihi,jlo,jhi;
    force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
    force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);
    
    const double rij = force->numeric(FLERR,arg[3]);
    const double k1 = (force->numeric(FLERR,arg[4]));
    const double k2 = force->numeric(FLERR,arg[5]);  
    const double k3 = force->numeric(FLERR,arg[6]); 
  
    int count = 0;
    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo,i); j <= jhi; j++) {
        cut_soc_landau[i][j] = rij;   
        K1[i][j] = k1/hbar;
	if (mech_flag) {
	  K1_mech[i][j] = k1;
	} else {
	  K1_mech[i][j] = 0.0;
	}
        K2[i][j] = k2;
        K3[i][j] = k3;
        setflag[i][j] = 1;
        count++;
      }
    }
    if (count == 0) error->all(FLERR,"Incorrect args in pair_style command");  
  } else error->all(FLERR,"Incorrect args in pair_style command");

}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSpinSocLandau::init_style()
{
  if (!atom->sp_flag || !atom->mumag_flag)
    error->all(FLERR,"Pair spin requires atom attributes sp, mumag");

  neighbor->request(this,instance_me);

  // check this half/full request  =>  to be corrected
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

double PairSpinSocLandau::init_one(int i, int j)
{
   
   if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cut_soc_global;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSpinSocLandau::write_restart(FILE *fp)
{
  write_restart_settings(fp);
  
  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        if (soc_neel_flag){
          fwrite(&K1[i][j],sizeof(double),1,fp);
          fwrite(&K1_mech[i][j],sizeof(double),1,fp);
          fwrite(&K2[i][j],sizeof(double),1,fp);
          fwrite(&K3[i][j],sizeof(double),1,fp);
          fwrite(&cut_soc_landau[i][j],sizeof(double),1,fp);
        }
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSpinSocLandau::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&K1[i][j],sizeof(double),1,fp);
          fread(&K1_mech[i][j],sizeof(double),1,fp);
          fread(&K2[i][j],sizeof(double),1,fp);
          fread(&K2[i][j],sizeof(double),1,fp);
          fread(&cut_soc_landau[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&K1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&K1_mech[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&K2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&K3[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_soc_landau[i][j],1,MPI_DOUBLE,0,world);
      }
    }
  }
}

 
/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSpinSocLandau::write_restart_settings(FILE *fp)
{
  fwrite(&cut_soc_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSpinSocLandau::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_soc_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_soc_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world); 
}
