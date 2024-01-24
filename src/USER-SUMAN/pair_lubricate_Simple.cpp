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
   Contributing authors: Randy Schunk (SNL)
   Amit Kumar and Michael Bybee (UIUC)
   Dave Heine (Corning), polydispersity
   ------------------------------------------------------------------------- */
/* Modified by Ranga on 27/05/2017 from the GRM formulation
   given by Kim and Karilla, Microhydrodynamics  */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_lubricate_Simple.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "fix_deform.h"
#include "memory.h"
#include "random_mars.h"
#include "fix_wall.h"
#include "input.h"
#include "variable.h"
#include "math_const.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

// same as fix_deform.cpp

enum{NO_REMAP,X_REMAP,V_REMAP};


// same as fix_wall.cpp

enum{EDGE,CONSTANT,VARIABLE};

/* ---------------------------------------------------------------------- */

PairLubricateSimple::PairLubricateSimple(LAMMPS *lmp) : PairLubricate(lmp)
{
  no_virial_fdotr_compute = 1;
}

/* ---------------------------------------------------------------------- */

PairLubricateSimple::~PairLubricateSimple()
{
  if (copymode) return;
}

/* ---------------------------------------------------------------------- */

void PairLubricateSimple::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz,tx,ty,tz;
  double rsq,r,h_sep,beta0,beta1,radi,radj;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3;
  double vt1,vt2,vt3;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double lamda[3],vstream[3];
  
  double vxmu2f = force->vxmu2f;
  
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;
  
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    radi = radius[i];
    
    // Drag contribution to force and torque due to isotropic terms
    // Drag contribution to stress from isotropic RS0
    
    if (flagfld) {

      domain->x2lamda(x[i],lamda);
      double *h_rate = domain->h_rate;
      double *h_ratelo = domain->h_ratelo;
      vstream[0] = h_rate[0]*lamda[0] + h_rate[5]*lamda[1] +
	h_rate[4]*lamda[2] + h_ratelo[0];
      vstream[1] = h_rate[1]*lamda[1] + h_rate[3]*lamda[2] + h_ratelo[1];
      vstream[2] = h_rate[2]*lamda[2] + h_ratelo[2];
      
      
      R0  = 6.0*MY_PI*mu;
      RT0 = 8.0*MY_PI*mu;
      RS0 = 20.0/3.0*MY_PI*mu;
      
      f[i][0] += vxmu2f*R0*radi*(vstream[0]-v[i][0]);
      f[i][1] += vxmu2f*R0*radi*(vstream[1]-v[i][1]);
      f[i][2] += vxmu2f*R0*radi*(vstream[2]-v[i][2]);
      
      const double radi3 = radi*radi*radi;
      
      torque[i][0] -= vxmu2f*RT0*radi3*(omega[i][0]+0.5*h_rate[3]/domain->zprd);
      torque[i][1] -= vxmu2f*RT0*radi3*(omega[i][1]-0.5*h_rate[4]/domain->zprd);
      torque[i][2] -= vxmu2f*RT0*radi3*(omega[i][2]+0.5*h_rate[5]/domain->yprd);
      
      // Ef = (grad(vstream) + (grad(vstream))^T) / 2
      // set Ef from h_rate in strain units
      if(shearing){
	Ef[0][0] = h_rate[0]/domain->xprd;
	Ef[1][1] = h_rate[1]/domain->yprd;
	Ef[2][2] = h_rate[2]/domain->zprd;
	Ef[0][1] = Ef[1][0] = 0.5 * h_rate[5]/domain->yprd;
	Ef[0][2] = Ef[2][0] = 0.5 * h_rate[4]/domain->zprd;
	Ef[1][2] = Ef[2][1] = 0.5 * h_rate[3]/domain->zprd;
	
	if (vflag_either) {
	  double vRS0 = -vxmu2f* RS0*radi3;
	  v_tally_tensor(i,i,nlocal,newton_pair,
	   		 vRS0*Ef[0][0],vRS0*Ef[1][1],vRS0*Ef[2][2],
	   		 vRS0*Ef[0][1],vRS0*Ef[0][2],vRS0*Ef[1][2]);
	}
      }

    }

    if (!flagHI)continue;
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      radj = atom->radius[j];
      double radsum = radi+radj; //Sum of two particle's radii
      
      if (rsq < cutsq[itype][jtype] and rsq < (1.45*1.45*radsum*radsum) ) {
	r = sqrt(rsq);

	// scalar resistances XA and YA
	
	h_sep = r - radsum;
	
	// if less than the minimum gap use the minimum gap instead
	
	if (r < cut_inner[itype][jtype]*radsum) //Scaled inner cutoff by particle radii
	  h_sep = cut_inner[itype][jtype]*radsum - radsum;
	
	double hinv,lhinv,XA11=0.0,YA11=0.0,YB11=0.0,YC12=0.0,ws[3],wd[3],nx,ny,nz;
	ws[0] = (omega[i][0]+omega[j][0])/2.0;
	ws[1] = (omega[i][1]+omega[j][1])/2.0;
	ws[2] = (omega[i][2]+omega[j][2])/2.0;
	
	wd[0] = (omega[i][0]-omega[j][0])/2.0;
	wd[1] = (omega[i][1]-omega[j][1])/2.0;
	wd[2] = (omega[i][2]-omega[j][2])/2.0;
	
	nx=-delx/r;ny=-dely/r;nz=-delz/r;
	
	
	hinv=(radi+radj)/(2.0*h_sep);
	lhinv=log(hinv);
	
	if(lhinv < 0) error->all(FLERR,"Using pair lubricate with cutoff problem: log(1/h) is negative");
	
	beta0=radj/radi;
	beta1=1.0+beta0;
	
	double b0p2,b1p2,b1p3,mupradi;
	b0p2=beta0*beta0;
	b1p2=beta1*beta1;
	b1p3=beta1*b1p2;
	mupradi=mu*MY_PI*radi;
	
	// scalar resistances
	if (flaglog) {
	  XA11=6.0*mupradi*( hinv*2.0*b0p2 + lhinv*beta0*(1.0+ 7.0*beta0+ b0p2 )/5.0 )/b1p3;
	  YA11=1.6*mupradi*lhinv*beta0*(2.0+beta0+2.0*b0p2)/b1p3;
	  YB11=-0.8*mupradi*radi*lhinv*beta0*(4.0+beta0)/b1p2;
	  YC12=0.8*mupradi*radi*radi*lhinv*b0p2/beta1*(1.0-4.0/beta0); //YC12*(1-4/beta0)
	} 
	else  XA11=12.0*mupradi*hinv*b0p2/b1p3;
	
	
	
	/*
	  if (flaglog) {
	  XA11=6.0*MY_PI*mu*radi*( hinv*2.0*pow(beta0,2.0) + lhinv*beta0*(1.0+ 7.0*beta0+ beta0*beta0)/5.0 )/pow(beta1,3.0);
	  YA11=1.6*MY_PI*mu*radi*lhinv*beta0*(2.0+beta0+2.0*beta0*beta0)/pow(beta1,3.0);
	  YB11=-0.8*MY_PI*mu*radi*radi*lhinv*beta0*(4.0+beta0)/(beta1*beta1);
	  YC12=0.8*MY_PI*mu*pow(radi,3.0)*lhinv*beta0*beta0/beta1*(1.0-4.0/beta0); //YC12*(1-4/beta0)
	  } 
	  else  XA11=12*MY_PI*mu*radi*hinv*pow(beta0,2.0)/pow(beta1,3.0);
	*/
	// Relative velocity components U^2-U^1
	vr1 = v[j][0] - v[i][0];
	vr2 = v[j][1] - v[i][1];
	vr3 = v[j][2] - v[i][2];
	
	// normal component (vr.n)n
	
	vnnr = vr1*nx + vr2*ny + vr3*nz;
	vn1 = vnnr*nx;
	vn2 = vnnr*ny;
	vn3 = vnnr*nz;
	
	// tangential component vr - (vr.n)n
	
	vt1 = vr1 - vn1;
	vt2 = vr2 - vn2;
	vt3 = vr3 - vn3;
	
	// force due to squeeze type motion
	//f=XA11 nn dot vr
	fx  = XA11*vn1;
	fy  = XA11*vn2;
	fz  = XA11*vn3;
	
	
	// force due to all shear kind of motions
	
	if (flaglog) {
	  double ybfac=(1.0-beta0*(1.0+4.0*beta0)/(4.0+beta0))*YB11;
	  //f+=YA11*vt1-(r1+r2)*YA11*ws cross n + ybfac* wd cross n
	  fx = fx + YA11*(vt1-(radi+radj)*(ws[1]*nz-ws[2]*ny))+ybfac*(wd[1]*nz-wd[2]*ny);
	  fy = fy + YA11*(vt2-(radi+radj)*(ws[2]*nx-ws[0]*nz))+ybfac*(wd[2]*nx-wd[0]*nz);
	  fz = fz + YA11*(vt3-(radi+radj)*(ws[0]*ny-ws[1]*nx))+ybfac*(wd[0]*ny-wd[1]*nx);
	}
	
	// scale forces for appropriate units
	
	fx *= vxmu2f;
	fy *= vxmu2f;
	fz *= vxmu2f;
	
	// add to total force
	
	f[i][0] += fx;
	f[i][1] += fy;
	f[i][2] += fz;

	if(newton_pair || j < nlocal) {
	  f[j][0] -= fx;
	  f[j][1] -= fy;
	  f[j][2] -= fz;
	}
	
	// torque due to this force
	
	if (flaglog) {
	  double wsdotn,wddotn;
	  wsdotn=ws[0]*nx+ws[1]*ny+ws[2]*nz;
	  wddotn=wd[0]*nx+wd[1]*ny+wd[2]*nz;
	  
	  //squeeze+shear+pump contributions to torque
	  //YB11*(vr cross n + (r1+r2)* tang(ws)) +YC12*(1-4/beta)*tang(wd)
	  double tx = YB11*(vr2*nz -vr3*ny + (radi+radj)*(ws[0]-wsdotn*nx));
	  double ty = YB11*(vr3*nx -vr1*nz + (radi+radj)*(ws[1]-wsdotn*ny));
	  double tz = YB11*(vr1*ny -vr2*nx + (radi+radj)*(ws[2]-wsdotn*nz));

	  double ttx = YC12*(wd[0]-wddotn*nx);
	  double tty = YC12*(wd[1]-wddotn*ny);
	  double ttz = YC12*(wd[2]-wddotn*nz);
	  
	  torque[i][0] += vxmu2f * (tx + ttx);
	  torque[i][1] += vxmu2f * (ty + tty);
	  torque[i][2] += vxmu2f * (tz + ttz);

	  if(newton_pair || j < nlocal) {
	    torque[j][0] += vxmu2f * (tx - ttx);
	    torque[j][1] += vxmu2f * (ty - tty);
	    torque[j][2] += vxmu2f * (tz - ttz);
	  }
	}
	
	// set j = nlocal so that only I gets tallied
	
	if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,0.0,0.0,fx,fy,fz,delx,dely,delz);
	//
	//Add up stresslet for particle i
	//v_tally_tensor(i,nlocal,nlocal,newton_pair,fx*delx,fy*dely,fz*delz,
	//0.5*(fx*dely+fy*delx),0.5*(fx*delz+fz*delx),0.5*(fy*delz+fz*dely));
	
      }
    }

  }
  
}

/* ----------------------------------------------------------------------
   init specific to this pair style
   ------------------------------------------------------------------------- */

void PairLubricateSimple::init_style()
{
  if (force->newton_pair == 1)
    error->all(FLERR,"Pair lubricate/poly requires newton pair off");
  if (comm->ghost_velocity == 0)
    error->all(FLERR,
	       "Pair lubricate/poly requires ghost atoms store velocity");
  if (!atom->sphere_flag)
    error->all(FLERR,"Pair lubricate/poly requires atom style sphere");
  
  // ensure all particles are finite-size
  // for pair hybrid, should limit test to types using the pair style
  
  double *radius = atom->radius;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++)
    if (radius[i] == 0.0)
      error->one(FLERR,"Pair lubricate/poly requires extended particles");
  
  // CHRIS:: IMPORTANT :: support for half list was added
  //neighbor->add_request(this, NeighConst::REQ_FULL);
  neighbor->add_request(this);
  
  // set the isotropic constants that depend on the volume fraction
  // vol_T = total volume
  
  // check for fix deform, if exists it must use "remap v"
  // If box will change volume, set appropriate flag so that volume
  // and v.f. corrections are re-calculated at every step.
  
  // Ranga: Volume fraction correction unnecessary for our purposes
  // if available volume is different from box volume
  // due to walls, set volume appropriately; if walls will
  // move, set appropriate flag so that volume and v.f. corrections
  // are re-calculated at every step.
  
  shearing = flagdeform = flagwall = 0;
  for (int i = 0; i < modify->nfix; i++){
    if (strcmp(modify->fix[i]->style,"deform") == 0 || strcmp(modify->fix[i]->style,"deform/kk") == 0) {
      shearing = flagdeform = 1;
      if (((FixDeform *) modify->fix[i])->remapflag != V_REMAP)
	error->all(FLERR,"Using pair lubricate with inconsistent "
		   "fix deform remap option");
    }
    if (strstr(modify->fix[i]->style,"wall") != NULL) {
      if (flagwall)
	error->all(FLERR,
		   "Cannot use multiple fix wall commands with "
		   "pair lubricate/poly");
      flagwall = 1; // Walls exist
      wallfix = (FixWall *) modify->fix[i];
      if (wallfix->xflag) flagwall = 2; // Moving walls exist
    }
    
    if (strstr(modify->fix[i]->style,"wall") != NULL){
      flagwall = 1; // Walls exist
      if (((FixWall *) modify->fix[i])->xflag ) {
	flagwall = 2; // Moving walls exist
	wallfix = (FixWall *) modify->fix[i];
      }
    }
  }
  
  // check for fix deform, if exists it must use "remap v"
  
  shearing = 0; // then why set it above??
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"deform") == 0 || strcmp(modify->fix[i]->style,"deform/kk") == 0) {
      shearing = 1;
      if (((FixDeform *) modify->fix[i])->remapflag != V_REMAP)
	error->all(FLERR,"Using pair lubricate/poly with inconsistent "
		   "fix deform remap option");
    }
  
  // set Ef = 0 since used whether shearing or not
  
  Ef[0][0] = Ef[0][1] = Ef[0][2] = 0.0;
  Ef[1][0] = Ef[1][1] = Ef[1][2] = 0.0;
  Ef[2][0] = Ef[2][1] = Ef[2][2] = 0.0;
}

/* ----------------------------------------------------------------------
   allocate all arrays; overiding here so Kokkos can overide
   ------------------------------------------------------------------------- */

void PairLubricateSimple::allocate()
{
  printf("Inside PairLubricateSimple::allocate()\n");
  PairLubricate::allocate();
}
