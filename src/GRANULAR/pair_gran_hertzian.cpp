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
   Contributing authors: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "string.h"
#include "pair_gran_hertzian.h"
#include "atom.h"
#include "force.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

PairGranHertzian::PairGranHertzian(LAMMPS *lmp) : PairGranHistory(lmp)
{
  history = 1;
}

/* ---------------------------------------------------------------------- */

void PairGranHertzian::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz;
  double radi,radj,radsum,rsq,r,rinv;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double xmeff,damp,ccel,tor1,tor2,tor3;
  double fn,fs,fs1,fs2,fs3;
  double shrmag,rsht,rhertz;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *shear,*allshear,**firstshear;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = list->listgranhistory->firstneigh;
  firstshear = list->listgranhistory->firstdouble;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    touch = firsttouch[i];
    allshear = firstshear[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;

      if (rsq >= radsum*radsum) {

	// unset non-touching neighbors

        touch[jj] = 0;
	shear = &allshear[3*jj];
        shear[0] = 0.0;
        shear[1] = 0.0;
        shear[2] = 0.0;

      } else {
	r = sqrt(rsq);

	// relative translational velocity

	vr1 = v[i][0] - v[j][0];
	vr2 = v[i][1] - v[j][1];
	vr3 = v[i][2] - v[j][2];

	vr1 *= dt;
	vr2 *= dt;
	vr3 *= dt;

	// normal component

	vnnr = vr1*delx + vr2*dely + vr3*delz;
	vn1 = delx*vnnr / rsq;
	vn2 = dely*vnnr / rsq;
	vn3 = delz*vnnr / rsq;

	// tangential component

	vt1 = vr1 - vn1;
	vt2 = vr2 - vn2;
	vt3 = vr3 - vn3;

	// relative rotational velocity

	wr1 = radi*omega[i][0] + radj*omega[j][0];
	wr2 = radi*omega[i][1] + radj*omega[j][1];
	wr3 = radi*omega[i][2] + radj*omega[j][2];

	wr1 *= dt/r;
	wr2 *= dt/r;
	wr3 *= dt/r;

	// normal damping term
	// this definition of DAMP includes the extra 1/r term

	xmeff = rmass[i]*rmass[j] / (rmass[i]+rmass[j]);
	if (mask[i] & freeze_group_bit) xmeff = rmass[j];
	if (mask[j] & freeze_group_bit) xmeff = rmass[i];
	damp = xmeff*gamman_dl*vnnr/rsq;
	ccel = xkk*(radsum-r)/r - damp;
	rhertz = sqrt(radsum - r);
	ccel = rhertz * ccel;

	// relative velocities

	vtr1 = vt1 - (delz*wr2-dely*wr3);
	vtr2 = vt2 - (delx*wr3-delz*wr1);
	vtr3 = vt3 - (dely*wr1-delx*wr2);
	vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
	vrel = sqrt(vrel);

	// shear history effects
	// shrmag = magnitude of shear

	touch[jj] = 1;
	shear = &allshear[3*jj];
        shear[0] += vtr1;
        shear[1] += vtr2;
        shear[2] += vtr3;
        shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] + 
		      shear[2]*shear[2]);

	// rotate shear displacements correctly

	rsht = shear[0]*delx + shear[1]*dely + shear[2]*delz;
	rsht /= rsq;
        shear[0] -= rsht*delx;
        shear[1] -= rsht*dely;
        shear[2] -= rsht*delz;

	// tangential forces

        fs1 = -rhertz * (xkkt*shear[0] + xmeff*gammas_dl*vtr1);
        fs2 = -rhertz * (xkkt*shear[1] + xmeff*gammas_dl*vtr2);
        fs3 = -rhertz * (xkkt*shear[2] + xmeff*gammas_dl*vtr3);

	// force normalization
	// rescale frictional displacements and forces if needed

	fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
	fn = xmu * fabs(ccel*r);

	if (fs > fn) {
	  if (shrmag != 0.0) {
	    shear[0] = (fn/fs) * (shear[0] + xmeff*gammas_dl*vtr1/xkkt) -
	      xmeff*gammas_dl*vtr1/xkkt;
	    shear[1] = (fn/fs) * (shear[1] + xmeff*gammas_dl*vtr2/xkkt) -
	      xmeff*gammas_dl*vtr2/xkkt;
	    shear[2] = (fn/fs) * (shear[2] + xmeff*gammas_dl*vtr3/xkkt) -
	      xmeff*gammas_dl*vtr3/xkkt;
	    fs1 *= fn/fs;
	    fs2 *= fn/fs;
	    fs3 *= fn/fs;
	  } else {
	    fs1 = 0.0;
	    fs2 = 0.0;
	    fs3 = 0.0;
	  }
	}

	// forces & torques

	fx = delx*ccel + fs1;
	fy = dely*ccel + fs2;
	fz = delz*ccel + fs3;
	f[i][0] += fx;
	f[i][1] += fy;
	f[i][2] += fz;

	rinv = 1/r;
	tor1 = rinv * (dely*fs3 - delz*fs2);
	tor2 = rinv * (delz*fs1 - delx*fs3);
	tor3 = rinv * (delx*fs2 - dely*fs1);
	torque[i][0] -= radi*tor1;
	torque[i][1] -= radi*tor2;
	torque[i][2] -= radi*tor3;

	if (j < nlocal) {
	  f[j][0] -= fx;
	  f[j][1] -= fy;
	  f[j][2] -= fz;
	  torque[j][0] -= radj*tor1;
	  torque[j][1] -= radj*tor2;
	  torque[j][2] -= radj*tor3;
	}

	if (evflag) ev_tally_xyz(i,j,nlocal,0,
				 0.0,0.0,fx,fy,fz,delx,dely,delz);
      }
    }
  }
}
