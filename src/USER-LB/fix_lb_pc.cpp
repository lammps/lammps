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
   Contributing authors: Frances Mackay, Santtu Ollila, Colin Denniston (UWO)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "string.h"
#include "fix_lb_pc.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "comm.h"
#include "domain.h"
#include "fix_lb_fluid.h"
#include "modify.h"
#include "mpi.h"
#include "group.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixLbPC::FixLbPC(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 3)
    error->all(FLERR,"Illegal fix lb/pc command");

  time_integrate = 1;

  // perform initial allocation of atom-based array
  // register with Atom class

  force_old = NULL;
  up = NULL;
  up_old = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  Gamma_MD = new double[atom->ntypes+1];

 int groupbit_lb_fluid = 0;
  for(int ifix=0; ifix<modify->nfix; ifix++)
    if(strcmp(modify->fix[ifix]->style,"lb/fluid")==0){
      fix_lb_fluid = (FixLbFluid *)modify->fix[ifix];
      groupbit_lb_fluid = group->bitmask[modify->fix[ifix]->igroup];
    }

  if(groupbit_lb_fluid == 0)
    error->all(FLERR,"the lb/fluid fix must also be used if using the lb/pc fix");     

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  for(int j=0; j<nlocal; j++){
    if((mask[j] & groupbit) && !(mask[j] & groupbit_lb_fluid))
      error->one(FLERR,"can only use the lb/pc fix for an atom if also using the lb/fluid fix for that atom");
  }
     
}

/* ---------------------------------------------------------------------- */

FixLbPC::~FixLbPC() {

  atom->delete_callback(id,0);
  
  memory->destroy(force_old);
  memory->destroy(up);
  memory->destroy(up_old);

  delete [] Gamma_MD;
}


/* ---------------------------------------------------------------------- */

int FixLbPC::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLbPC::init()
{
  double *Gamma = fix_lb_fluid->Gamma;
  double dm_lb = fix_lb_fluid->dm_lb;
  double dt_lb = fix_lb_fluid->dt_lb;

  MPI_Comm_rank(world,&me);
  
  dtv = update->dt;
  dtf = update->dt * force->ftm2v;

  for(int i=0; i<=atom->ntypes; i++)
    Gamma_MD[i] = Gamma[i]*dm_lb/dt_lb;

}

/* ---------------------------------------------------------------------- */
void FixLbPC::initial_integrate(int vflag) {

  double dtfm;
  
  double **x = atom->x;
  double dx[3];
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  compute_up();

  for(int i=0; i<nlocal; i++){
    up_old[i][0] = up[i][0];
    up_old[i][1] = up[i][1];
    up_old[i][2] = up[i][2];
    force_old[i][0] = f[i][0];
    force_old[i][1] = f[i][1];
    force_old[i][2] = f[i][2];
  }

   
  if(rmass){
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	dtfm = dtf/rmass[i];
	expminusdttimesgamma = exp(-dtv*Gamma_MD[type[i]]/rmass[i]);	
	
	dx[0] = dtv*v[i][0] + 0.5*(f[i][0]*force->ftm2v - Gamma_MD[type[i]]*(v[i][0]-up[i][0]))*dtv*dtv/rmass[i];
	dx[1] = dtv*v[i][1] + 0.5*(f[i][1]*force->ftm2v - Gamma_MD[type[i]]*(v[i][1]-up[i][1]))*dtv*dtv/rmass[i];
	dx[2] = dtv*v[i][2] + 0.5*(f[i][2]*force->ftm2v - Gamma_MD[type[i]]*(v[i][2]-up[i][2]))*dtv*dtv/rmass[i];

	x[i][0] += dx[0];
	x[i][1] += dx[1];
	x[i][2] += dx[2];
	
	// Approximation for v
	if(Gamma_MD[type[i]] == 0.0){
	  v[i][0] += f[i][0]*dtfm;
	  v[i][1] += f[i][1]*dtfm;
	  v[i][2] += f[i][2]*dtfm;
	}else{
 	  v[i][0] = (v[i][0]-up[i][0]-f[i][0]*force->ftm2v/Gamma_MD[type[i]])*expminusdttimesgamma + 
 	    f[i][0]*force->ftm2v/Gamma_MD[type[i]] + up[i][0];
 	  v[i][1] = (v[i][1]-up[i][1]-f[i][1]*force->ftm2v/Gamma_MD[type[i]])*expminusdttimesgamma + 
 	    f[i][1]*force->ftm2v/Gamma_MD[type[i]] + up[i][1];
 	  v[i][2] = (v[i][2]-up[i][2]-f[i][2]*force->ftm2v/Gamma_MD[type[i]])*expminusdttimesgamma + 
 	    f[i][2]*force->ftm2v/Gamma_MD[type[i]] + up[i][2];
	}	  
      }
    }

  } else {
    // this does NOT take varying masses into account
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	dtfm = dtf/mass[type[i]];
	expminusdttimesgamma = exp(-dtv*Gamma_MD[type[i]]/mass[type[i]]);	

	dx[0] = dtv*v[i][0] + 0.5*(f[i][0]*force->ftm2v - Gamma_MD[type[i]]*(v[i][0]-up[i][0]))*dtv*dtv/mass[type[i]];
	dx[1] = dtv*v[i][1] + 0.5*(f[i][1]*force->ftm2v - Gamma_MD[type[i]]*(v[i][1]-up[i][1]))*dtv*dtv/mass[type[i]];
	dx[2] = dtv*v[i][2] + 0.5*(f[i][2]*force->ftm2v - Gamma_MD[type[i]]*(v[i][2]-up[i][2]))*dtv*dtv/mass[type[i]];

	x[i][0] += dx[0];
	x[i][1] += dx[1];
	x[i][2] += dx[2];
	
	// Approximation for v
	if(Gamma_MD[type[i]] == 0.0){
	  v[i][0] += f[i][0]*dtfm;
	  v[i][1] += f[i][1]*dtfm;
	  v[i][2] += f[i][2]*dtfm;
	}else{
 	  v[i][0] = (v[i][0]-up[i][0]-f[i][0]*force->ftm2v/Gamma_MD[type[i]])*expminusdttimesgamma + 
 	    f[i][0]*force->ftm2v/Gamma_MD[type[i]] + up[i][0];
 	  v[i][1] = (v[i][1]-up[i][1]-f[i][1]*force->ftm2v/Gamma_MD[type[i]])*expminusdttimesgamma + 
 	    f[i][1]*force->ftm2v/Gamma_MD[type[i]] + up[i][1];
 	  v[i][2] = (v[i][2]-up[i][2]-f[i][2]*force->ftm2v/Gamma_MD[type[i]])*expminusdttimesgamma + 
 	    f[i][2]*force->ftm2v/Gamma_MD[type[i]] + up[i][2];
	}
      }
    }
  }
}

/* ---------------------------------------------------------------------- */
void FixLbPC::final_integrate()
{
  double dtfm;
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  compute_up();
  
  if(rmass){
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	dtfm = dtf/rmass[i];
	expminusdttimesgamma = exp(-dtv*Gamma_MD[type[i]]/rmass[i]);
	DMDcoeff = (dtv - rmass[i]*(1.0-expminusdttimesgamma)/Gamma_MD[type[i]]);
	
	if(Gamma_MD[type[i]] == 0.0){
	  v[i][0] += 0.5*(f[i][0] - force_old[i][0])*dtfm;
	  v[i][1] += 0.5*(f[i][1] - force_old[i][1])*dtfm;
	  v[i][2] += 0.5*(f[i][2] - force_old[i][2])*dtfm;
	}else{
	  v[i][0] += DMDcoeff*((f[i][0] - force_old[i][0])*force->ftm2v/Gamma_MD[type[i]] + up[i][0] - up_old[i][0])/dtv;
	  v[i][1] += DMDcoeff*((f[i][1] - force_old[i][1])*force->ftm2v/Gamma_MD[type[i]] + up[i][1] - up_old[i][1])/dtv;
	  v[i][2] += DMDcoeff*((f[i][2] - force_old[i][2])*force->ftm2v/Gamma_MD[type[i]] + up[i][2] - up_old[i][2])/dtv;
	}
	

      }
    }
  } else {
    // this does NOT take varying masses into account
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	dtfm = dtf/mass[type[i]];
	expminusdttimesgamma = exp(-dtv*Gamma_MD[type[i]]/mass[type[i]]);
	DMDcoeff = (dtv - mass[type[i]]*(1.0-expminusdttimesgamma)/Gamma_MD[type[i]]);
       
	if(Gamma_MD[type[i]] == 0.0){
	  v[i][0] += 0.5*(f[i][0] - force_old[i][0])*dtfm;
	  v[i][1] += 0.5*(f[i][1] - force_old[i][1])*dtfm;
	  v[i][2] += 0.5*(f[i][2] - force_old[i][2])*dtfm;
	}else{
	  v[i][0] += DMDcoeff*((f[i][0] - force_old[i][0])*force->ftm2v/Gamma_MD[type[i]] + up[i][0] - up_old[i][0])/dtv;
	  v[i][1] += DMDcoeff*((f[i][1] - force_old[i][1])*force->ftm2v/Gamma_MD[type[i]] + up[i][1] - up_old[i][1])/dtv;
	  v[i][2] += DMDcoeff*((f[i][2] - force_old[i][2])*force->ftm2v/Gamma_MD[type[i]] + up[i][2] - up_old[i][2])/dtv;
	}	

      }
    }
  }

}


/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixLbPC::grow_arrays(int nmax)
{

  memory->grow(force_old,nmax,3,"FixLbPC:force_old");
  memory->grow(up_old,nmax,3,"FixLbPC:up_old");
  memory->grow(up,nmax,3,"FixLbPC:up");

}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixLbPC::copy_arrays(int i, int j, int delflag)
{

  force_old[j][0] = force_old[i][0];
  force_old[j][1] = force_old[i][1];
  force_old[j][2] = force_old[i][2];
  up_old[j][0] = up_old[i][0];
  up_old[j][1] = up_old[i][1];
  up_old[j][2] = up_old[i][2];
  up[j][0] = up[i][0];
  up[j][1] = up[i][1];
  up[j][2] = up[i][2];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixLbPC::pack_exchange(int i, double *buf)
{

  buf[0] = force_old[i][0];
  buf[1] = force_old[i][1];
  buf[2] = force_old[i][2];
  buf[3] = up_old[i][0];
  buf[4] = up_old[i][1];
  buf[5] = up_old[i][2];
  buf[6] = up[i][0];
  buf[7] = up[i][1];
  buf[8] = up[i][2];

  return 9;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixLbPC::unpack_exchange(int nlocal, double *buf)
{

  force_old[nlocal][0] = buf[0];
  force_old[nlocal][1] = buf[1];
  force_old[nlocal][2] = buf[2];
  up_old[nlocal][0] = buf[3];
  up_old[nlocal][1] = buf[4];
  up_old[nlocal][2] = buf[5];
  up[nlocal][0] = buf[6];
  up[nlocal][1] = buf[7];
  up[nlocal][2] = buf[8];

  return 9;
}

/* ---------------------------------------------------------------------- */
 void FixLbPC::compute_up(void)
 {
   int *mask = atom->mask;
   int nlocal = atom->nlocal;
   double **x = atom->x;
   int i,k;
   int ix,iy,iz;
   int ixp,iyp,izp;
   double dx1,dy1,dz1;
   int isten,ii,jj,kk;
   double r,rsq,weightx,weighty,weightz;
   double ****u_lb = fix_lb_fluid->u_lb;
   int subNbx = fix_lb_fluid->subNbx;
   int subNby = fix_lb_fluid->subNby;
   int subNbz = fix_lb_fluid->subNbz;
   double dx_lb = fix_lb_fluid->dx_lb;
   double dt_lb = fix_lb_fluid->dt_lb;
   double FfP[64];
   int trilinear_stencil = fix_lb_fluid->trilinear_stencil;
   
   for(i=0; i<nlocal; i++){
    if(mask[i] & groupbit){

      //Calculate nearest leftmost grid point.
      //Since array indices from 1 to subNb-2 correspond to the
      // local subprocessor domain (not indices from 0), use the 
      // ceiling value.
      ix = (int)ceil((x[i][0]-domain->sublo[0])/dx_lb);
      iy = (int)ceil((x[i][1]-domain->sublo[1])/dx_lb);
      iz = (int)ceil((x[i][2]-domain->sublo[2])/dx_lb);
      
      //Calculate distances to the nearest points.
      dx1 = x[i][0] - (domain->sublo[0] + (ix-1)*dx_lb);
      dy1 = x[i][1] - (domain->sublo[1] + (iy-1)*dx_lb);
      dz1 = x[i][2] - (domain->sublo[2] + (iz-1)*dx_lb);

      // Need to convert these to lattice units:
      dx1 = dx1/dx_lb;
      dy1 = dy1/dx_lb;
      dz1 = dz1/dx_lb;
   
      up[i][0]=0.0; up[i][1]=0.0; up[i][2]=0.0;
      if(trilinear_stencil==0){
	isten=0;
	for(ii=-1; ii<3; ii++){
	  rsq=(-dx1+ii)*(-dx1+ii);
	  
	  if(rsq>=4)
	    weightx=0.0;
	  else{
	    r=sqrt(rsq);
	    if(rsq>1){
	      weightx=(5.0-2.0*r-sqrt(-7.0+12.0*r-4.0*rsq))/8.;
	    } else{
	      weightx=(3.0-2.0*r+sqrt(1.0+4.0*r-4.0*rsq))/8.;
	    }
	  }
	  for(jj=-1; jj<3; jj++){
	    rsq=(-dy1+jj)*(-dy1+jj);
	    if(rsq>=4)
	      weighty=0.0;
	    else{
	      r=sqrt(rsq);
	      if(rsq>1){
		weighty=(5.0-2.0*r-sqrt(-7.0+12.0*r-4.0*rsq))/8.;
	      } else{
		weighty=(3.0-2.0*r+sqrt(1.0+4.0*r-4.0*rsq))/8.;
	      }
	    }
	    for(kk=-1; kk<3; kk++){
	      rsq=(-dz1+kk)*(-dz1+kk);
	      if(rsq>=4)
		weightz=0.0;
	      else{
		r=sqrt(rsq);
		if(rsq>1){
		  weightz=(5.0-2.0*r-sqrt(-7.0+12.0*r-4.0*rsq))/8.;
		} else{
		  weightz=(3.0-2.0*r+sqrt(1.0+4.0*r-4.0*rsq))/8.;
		}
	      }
	      ixp = ix+ii;
	      iyp = iy+jj;
	      izp = iz+kk;
	      
	      
	      if(ixp==-1) ixp=subNbx+2;
	      if(iyp==-1) iyp=subNby+2;
	      if(izp==-1) izp=subNbz+2;
	      
	      FfP[isten] = weightx*weighty*weightz;
	      // interpolated velocity based on delta function.
	      for(k=0; k<3; k++){
		up[i][k] += u_lb[ixp][iyp][izp][k]*FfP[isten];
	      }
	    }
	  }
	}
      }else{
	FfP[0] = (1.-dx1)*(1.-dy1)*(1.-dz1);
	FfP[1] = (1.-dx1)*(1.-dy1)*dz1;
	FfP[2] = (1.-dx1)*dy1*(1.-dz1);
	FfP[3] = (1.-dx1)*dy1*dz1;
	FfP[4] = dx1*(1.-dy1)*(1.-dz1);
	FfP[5] = dx1*(1.-dy1)*dz1;
	FfP[6] = dx1*dy1*(1.-dz1);
	FfP[7] = dx1*dy1*dz1;
	
	ixp = (ix+1);
	iyp = (iy+1);
	izp = (iz+1);
	
	for (k=0; k<3; k++) { 	// tri-linearly interpolated velocity at node
	  up[i][k] = u_lb[ix][iy][iz][k]*FfP[0]
	    + u_lb[ix][iy][izp][k]*FfP[1]
	    + u_lb[ix][iyp][iz][k]*FfP[2]
	    + u_lb[ix][iyp][izp][k]*FfP[3]
	    + u_lb[ixp][iy][iz][k]*FfP[4]
	    + u_lb[ixp][iy][izp][k]*FfP[5]
	    + u_lb[ixp][iyp][iz][k]*FfP[6]
	    + u_lb[ixp][iyp][izp][k]*FfP[7];	  
	}
      }    
      for(k=0; k<3; k++)
	up[i][k] = up[i][k]*dx_lb/dt_lb;

    }
  }
 }
