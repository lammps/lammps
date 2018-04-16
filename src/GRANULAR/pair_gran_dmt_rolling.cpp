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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pair_gran_dmt_rolling.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "fix.h"
#include "fix_neigh_history.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define TWOTHIRDS 0.6666666666666666
#define EPSILON 1e-10

enum {TSUJI, BRILLIANTOV};
enum {INDEP, BRILLROLL};

/* ---------------------------------------------------------------------- */

PairGranDMTRolling::PairGranDMTRolling(LAMMPS *lmp) :
  PairGranHookeHistory(lmp, 7),
  E_one(0), G_one(0), pois(0), muS_one(0), cor(0), alpha_one(0),
  Ecoh_one(0), kR_one(0), muR_one(0), etaR_one(0)
{
  int ntypes = atom->ntypes;
  memory->create(E,ntypes+1,ntypes+1,"pair:E");
  memory->create(G,ntypes+1,ntypes+1,"pair:G");
  memory->create(alpha,ntypes+1,ntypes+1,"pair:alpha");
  memory->create(gamman,ntypes+1,ntypes+1,"pair:gamman");
  memory->create(muS,ntypes+1,ntypes+1,"pair:muS");
  memory->create(Ecoh,ntypes+1,ntypes+1,"pair:Ecoh");
  memory->create(kR,ntypes+1,ntypes+1,"pair:kR");
  memory->create(muR,ntypes+1,ntypes+1,"pair:muR");
  memory->create(etaR,ntypes+1,ntypes+1,"pair:etaR");
}

/* ---------------------------------------------------------------------- */
PairGranDMTRolling::~PairGranDMTRolling()
{
  delete [] E_one;
  delete [] G_one;
  delete [] pois;
  delete [] muS_one;
  delete [] cor;
  delete [] alpha_one;
  delete [] Ecoh_one;
  delete [] kR_one;
  delete [] muR_one;
  delete [] etaR_one;
  //TODO: Make all this work with standard pair coeff type commands.
  //Also these should not be in the destructor.
  memory->destroy(E);
  memory->destroy(G);
  memory->destroy(alpha);
  memory->destroy(gamman);
  memory->destroy(muS);
  memory->destroy(Ecoh);
  memory->destroy(kR);
  memory->destroy(muR);
  memory->destroy(etaR);
}
/* ---------------------------------------------------------------------- */

void PairGranDMTRolling::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum;
  int itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz,nx,ny,nz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv,R,a;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double kn, kt, k_Q, k_R, eta_N, eta_T, eta_Q, eta_R;
  double Fhz, Fdamp, Fdmt, Fne, Fntot, Fscrit, Frcrit;
  double overlap;
  double mi,mj,meff,damp,ccel,tor1,tor2,tor3;
  double relrot1,relrot2,relrot3,vrl1,vrl2,vrl3,vrlmag,vrlmaginv;
  double rollmag, rolldotn, scalefac;
  double fr, fr1, fr2, fr3;
  double signtwist, magtwist, magtortwist, Mtcrit;
  double fs,fs1,fs2,fs3,roll1,roll2,roll3,torroll1,torroll2,torroll3;
  double tortwist1, tortwist2, tortwist3;
  double shrmag,rsht;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *shear,*allshear,**firstshear;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  int shearupdate = 1;
  if (update->setupflag) shearupdate = 0;

  // update rigid body info for owned & ghost atoms if using FixRigid masses
  // body[i] = which body atom I is in, -1 if none
  // mass_body = mass of each rigid body

  if (fix_rigid && neighbor->ago == 0){
      int tmp;
      int *body = (int *) fix_rigid->extract("body",tmp);
      double *mass_body = (double *) fix_rigid->extract("masstotal",tmp);
      if (atom->nmax > nmax) {
	  memory->destroy(mass_rigid);
	  nmax = atom->nmax;
	  memory->create(mass_rigid,nmax,"pair:mass_rigid");
      }
      int nlocal = atom->nlocal;
      for (i = 0; i < nlocal; i++)
	if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
	else mass_rigid[i] = 0.0;
      comm->forward_comm_pair(this);
  }

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = fix_history->firstflag;
  firstshear = fix_history->firstvalue;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      itype = type[i];
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
	  jtype = type[j];
	  j &= NEIGHMASK;

	  delx = xtmp - x[j][0];
	  dely = ytmp - x[j][1];
	  delz = ztmp - x[j][2];
	  rsq = delx*delx + dely*dely + delz*delz;
	  radj = radius[j];
	  radsum = radi + radj;

	  if (rsq >= radsum*radsum){
	      // unset non-touching neighbors
	      touch[jj] = 0;
	      shear = &allshear[size_history*jj];
	      for (int k = 0; k < size_history; k++)
		shear[k] = 0.0;
	  } else {
	      r = sqrt(rsq);
	      rinv = 1.0/r;
	      rsqinv = 1.0/rsq;
	      R = radi*radj/(radi+radj);
	      nx = delx*rinv;
	      ny = dely*rinv;
	      nz = delz*rinv;

	      // relative translational velocity

	      vr1 = v[i][0] - v[j][0];
	      vr2 = v[i][1] - v[j][1];
	      vr3 = v[i][2] - v[j][2];

	      // normal component

	      vnnr = vr1*nx + vr2*ny + vr3*nz; //v_R . n
	      vn1 = nx*vnnr;
	      vn2 = ny*vnnr;
	      vn3 = nz*vnnr;

	      // meff = effective mass of pair of particles
	      // if I or J part of rigid body, use body mass
	      // if I or J is frozen, meff is other particle

	      mi = rmass[i];
	      mj = rmass[j];
	      if (fix_rigid) {
		  if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
		  if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
	      }

	      meff = mi*mj / (mi+mj);
	      if (mask[i] & freeze_group_bit) meff = mj;
	      if (mask[j] & freeze_group_bit) meff = mi;

	      //****************************************
	      //Normal force = Hertzian contact + DMT + damping
	      //****************************************
	      overlap = radsum - r;
	      a = sqrt(R*overlap);
	      kn = 4.0/3.0*E[itype][jtype]*a;
	      Fhz = kn*overlap;

	      //Damping (based on Tsuji et al)
	      if (normaldamp == BRILLIANTOV) eta_N = a*meff*gamman[itype][jtype];
	      else if (normaldamp == TSUJI) eta_N=alpha[itype][jtype]*sqrt(meff*kn);

	      Fdamp = -eta_N*vnnr; //F_nd eq 23 and Zhao eq 19

	      //DMT
	      Fdmt = -4*MY_PI*Ecoh[itype][jtype]*R;

	      Fne = Fhz + Fdmt;
	      Fntot = Fne + Fdamp;

	      //****************************************
	      //Tangential force, including shear history effects
	      //****************************************

	      // tangential component
	      vt1 = vr1 - vn1;
	      vt2 = vr2 - vn2;
	      vt3 = vr3 - vn3;

	      // relative rotational velocity
	      //Luding Gran Matt 2008, v10,p235 suggests correcting radi and radj by subtracting
	      //delta/2, i.e. instead of radi, use distance to center of contact point?
	      wr1 = (radi*omega[i][0] + radj*omega[j][0]);
	      wr2 = (radi*omega[i][1] + radj*omega[j][1]);
	      wr3 = (radi*omega[i][2] + radj*omega[j][2]);

	      // relative tangential velocities
	      vtr1 = vt1 - (nz*wr2-ny*wr3);
	      vtr2 = vt2 - (nx*wr3-nz*wr1);
	      vtr3 = vt3 - (ny*wr1-nx*wr2);
	      vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
	      vrel = sqrt(vrel);

	      // shear history effects
	      touch[jj] = 1;
	      shear = &allshear[size_history*jj];
	      shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] +
			    shear[2]*shear[2]);

	      // Rotate and update shear displacements.
	      // See e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
	      if (shearupdate) {
		  rsht = shear[0]*nx + shear[1]*ny + shear[2]*nz;
		  if (fabs(rsht) < EPSILON) rsht = 0;
		  if (rsht > 0){
		      scalefac = shrmag/(shrmag - rsht); //if rhst == shrmag, contacting pair has rotated 90 deg. in one step, in which case you deserve a crash!
		      shear[0] -= rsht*nx;
		      shear[1] -= rsht*ny;
		      shear[2] -= rsht*nz;
		      //Also rescale to preserve magnitude
		      shear[0] *= scalefac;
		      shear[1] *= scalefac;
		      shear[2] *= scalefac;
		  }
		  //Update shear history
		  shear[0] += vtr1*dt;
		  shear[1] += vtr2*dt;
		  shear[2] += vtr3*dt;
	      }

	      // tangential forces = shear + tangential velocity damping
	      // following Zhao and Marshall Phys Fluids v20, p043302 (2008)
	      kt=8.0*G[itype][jtype]*a;

	      eta_T = eta_N; //Based on discussion in Marshall; eta_T can also be an independent parameter
	      fs1 = -kt*shear[0] - eta_T*vtr1; //eq 26
	      fs2 = -kt*shear[1] - eta_T*vtr2;
	      fs3 = -kt*shear[2] - eta_T*vtr3;

	      // rescale frictional displacements and forces if needed
	      Fscrit = muS[itype][jtype] * fabs(Fne);
	      // For JKR, use eq 43 of Marshall. For DMT, use Fne instead
	      shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] +
			    shear[2]*shear[2]);
	      fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
	      if (fs > Fscrit) {
		  if (shrmag != 0.0) {
		      //shear[0] = (Fcrit/fs) * (shear[0] + eta_T*vtr1/kt) - eta_T*vtr1/kt;
		      //shear[1] = (Fcrit/fs) * (shear[1] + eta_T*vtr1/kt) - eta_T*vtr1/kt;
		      //shear[2] = (Fcrit/fs) * (shear[2] + eta_T*vtr1/kt) - eta_T*vtr1/kt;
		      shear[0] = -1.0/kt*(Fscrit*fs1/fs + eta_T*vtr1); //Same as above, but simpler (check!)
		      shear[1] = -1.0/kt*(Fscrit*fs2/fs + eta_T*vtr2);
		      shear[2] = -1.0/kt*(Fscrit*fs3/fs + eta_T*vtr3);
		      fs1 *= Fscrit/fs;
		      fs2 *= Fscrit/fs;
		      fs3 *= Fscrit/fs;
		  } else fs1 = fs2 = fs3 = 0.0;
	      }

	      //****************************************
	      // Rolling force, including shear history effects
	      //****************************************

	      relrot1 = omega[i][0] - omega[j][0];
	      relrot2 = omega[i][1] - omega[j][1];
	      relrot3 = omega[i][2] - omega[j][2];

	      // rolling velocity, see eq. 31 of Wang et al, Particuology v 23, p 49 (2015)
	      // This is different from the Marshall papers, which use the Bagi/Kuhn formulation
	      // for rolling velocity (see Wang et al for why the latter is wrong)
	      vrl1 = R*(relrot2*nz - relrot3*ny); //- 0.5*((radj-radi)/radsum)*vtr1;
	      vrl2 = R*(relrot3*nx - relrot1*nz); //- 0.5*((radj-radi)/radsum)*vtr2;
	      vrl3 = R*(relrot1*ny - relrot2*nx); //- 0.5*((radj-radi)/radsum)*vtr3;
	      vrlmag = sqrt(vrl1*vrl1+vrl2*vrl2+vrl3*vrl3);
	      if (vrlmag != 0.0) vrlmaginv = 1.0/vrlmag;
	      else vrlmaginv = 0.0;

	      // Rolling displacement
	      rollmag = sqrt(shear[3]*shear[3] + shear[4]*shear[4] + shear[5]*shear[5]);
	      rolldotn = shear[3]*nx + shear[4]*ny + shear[5]*nz;

	      if (shearupdate) {
		  if (fabs(rolldotn) < EPSILON) rolldotn = 0;
		  if (rolldotn > 0){ //Rotate into tangential plane
		      scalefac = rollmag/(rollmag - rolldotn);
		      shear[3] -= rolldotn*nx;
		      shear[4] -= rolldotn*ny;
		      shear[5] -= rolldotn*nz;
		      //Also rescale to preserve magnitude
		      shear[3] *= scalefac;
		      shear[4] *= scalefac;
		      shear[5] *= scalefac;
		  }
		  shear[3] += vrl1*dt;
		  shear[4] += vrl2*dt;
		  shear[5] += vrl3*dt;
	      }

	      k_R = kR[itype][jtype];
	      if (rollingdamp == INDEP) eta_R = etaR[itype][jtype];
	      else if (rollingdamp == BRILLROLL) eta_R = muR[itype][jtype]*fabs(Fne);
	      fr1 = -k_R*shear[3] - eta_R*vrl1;
	      fr2 = -k_R*shear[4] - eta_R*vrl2;
	      fr3 = -k_R*shear[5] - eta_R*vrl3;

	      // rescale frictional displacements and forces if needed
	      Frcrit = muR[itype][jtype] * fabs(Fne);

	      fr = sqrt(fr1*fr1 + fr2*fr2 + fr3*fr3);
	      if (fr > Frcrit) {
		  if (rollmag != 0.0) {
		      shear[3] = -1.0/k_R*(Frcrit*fr1/fr + eta_R*vrl1);
		      shear[4] = -1.0/k_R*(Frcrit*fr2/fr + eta_R*vrl2);
		      shear[5] = -1.0/k_R*(Frcrit*fr3/fr + eta_R*vrl3);
		      fr1 *= Frcrit/fr;
		      fr2 *= Frcrit/fr;
		      fr3 *= Frcrit/fr;
		  } else fr1 = fr2 = fr3 = 0.0;
	      }


	      //****************************************
	      // Twisting torque, including shear history effects
	      //****************************************
	      magtwist = relrot1*nx + relrot2*ny + relrot3*nz; //Omega_T (eq 29 of Marshall)
	      shear[6] += magtwist*dt;
	      k_Q = 0.5*kt*a*a;; //eq 32
	      eta_Q = 0.5*eta_T*a*a;
	      magtortwist = -k_Q*shear[6] - eta_Q*magtwist;//M_t torque (eq 30)

	      signtwist = (magtwist > 0) - (magtwist < 0);
	      Mtcrit=TWOTHIRDS*a*Fscrit;//critical torque (eq 44)
	      if (fabs(magtortwist) > Mtcrit){
		  shear[6] = 1.0/k_Q*(Mtcrit*signtwist - eta_Q*magtwist);
		  magtortwist = -Mtcrit * signtwist; //eq 34
	      }

	      // Apply forces & torques

	      fx = nx*Fntot + fs1;
	      fy = ny*Fntot + fs2;
	      fz = nz*Fntot + fs3;

	      f[i][0] += fx;
	      f[i][1] += fy;
	      f[i][2] += fz;

	      tor1 = ny*fs3 - nz*fs2;
	      tor2 = nz*fs1 - nx*fs3;
	      tor3 = nx*fs2 - ny*fs1;

	      torque[i][0] -= radi*tor1;
	      torque[i][1] -= radi*tor2;
	      torque[i][2] -= radi*tor3;

	      tortwist1 = magtortwist * nx;
	      tortwist2 = magtortwist * ny;
	      tortwist3 = magtortwist * nz;

	      torque[i][0] += tortwist1;
	      torque[i][1] += tortwist2;
	      torque[i][2] += tortwist3;

	      torroll1 = R*(ny*fr3 - nz*fr2); //n cross fr
	      torroll2 = R*(nz*fr1 - nx*fr3);
	      torroll3 = R*(nx*fr2 - ny*fr1);

	      torque[i][0] += torroll1;
	      torque[i][1] += torroll2;
	      torque[i][2] += torroll3;

	      if (force->newton_pair || j < nlocal) {
		  f[j][0] -= fx;
		  f[j][1] -= fy;
		  f[j][2] -= fz;

		  torque[j][0] -= radj*tor1;
		  torque[j][1] -= radj*tor2;
		  torque[j][2] -= radj*tor3;

		  torque[j][0] -= tortwist1;
		  torque[j][1] -= tortwist2;
		  torque[j][2] -= tortwist3;

		  torque[j][0] -= torroll1;
		  torque[j][1] -= torroll2;
		  torque[j][2] -= torroll3;
	      }
	      if (evflag) ev_tally_xyz(i,j,nlocal,0,
				       0.0,0.0,fx,fy,fz,delx,dely,delz);
	  }
      }
  }
}

/* ----------------------------------------------------------------------
   global settings
   ------------------------------------------------------------------------- */

void PairGranDMTRolling::settings(int narg, char **arg)
{
  if (narg < 6) error->all(FLERR,"Illegal pair_style command");

  int ntypes = atom->ntypes;

  if (narg < 8*ntypes) error->all(FLERR,"Illegal pair_style command");

  E_one = new double[ntypes+1];
  G_one = new double[ntypes+1];
  pois = new double[ntypes+1];
  muS_one = new double[ntypes+1];
  cor = new double[ntypes+1];
  alpha_one = new double[ntypes+1];
  Ecoh_one = new double[ntypes+1];
  kR_one = new double[ntypes+1];
  muR_one = new double[ntypes+1];
  etaR_one = new double[ntypes+1];

  for (int i=0; i < ntypes;i++){
      E_one[i+1] = force->numeric(FLERR, arg[i]);
      G_one[i+1] = force->numeric(FLERR, arg[ntypes+i]);
      muS_one[i+1] = force->numeric(FLERR, arg[2*ntypes+i]);
      cor[i+1] = force->numeric(FLERR, arg[3*ntypes+i]);
      Ecoh_one[i+1] = force->numeric(FLERR, arg[4*ntypes+i]);
      kR_one[i+1] = force->numeric(FLERR, arg[5*ntypes+i]);
      muR_one[i+1] = force->numeric(FLERR, arg[6*ntypes+i]);
      etaR_one[i+1] = force->numeric(FLERR, arg[7*ntypes+i]);
  }

  //Defaults
  normaldamp = TSUJI;
  rollingdamp = INDEP;

  int iarg = 8*ntypes;
  while (iarg < narg){
      if (strcmp(arg[iarg],"normaldamp") == 0){
	  if (iarg+2 > narg) error->all(FLERR, "Invalid pair/gran/dmt/rolling entry");
	  if (strcmp(arg[iarg+1],"tsuji") == 0) normaldamp = TSUJI;
	  else if (strcmp(arg[iarg+1],"brilliantov") == 0) normaldamp = BRILLIANTOV;
	  else error->all(FLERR, "Invalid normal damping model for pair/gran/dmt/rolling");
	  iarg += 2;
      }
      else if (strcmp(arg[iarg],"rollingdamp") == 0){
      	  if (iarg+2 > narg) error->all(FLERR, "Invalid pair/gran/dmt/rolling entry");
      	  if (strcmp(arg[iarg+1],"independent") == 0) rollingdamp = INDEP;
      	  else if (strcmp(arg[iarg+1],"brilliantov") == 0) rollingdamp = BRILLROLL;
      	  else error->all(FLERR, "Invalid rolling damping model for pair/gran/dmt/rolling");
      	  iarg += 2;
      }
      else{
	  iarg +=1;
      }
  }

  //Derived from inputs
  for (int i=1; i <= ntypes; i++){
      pois[i] = E_one[i]/(2.0*G_one[i]) - 1.0;
      alpha_one[i] = 1.2728-4.2783*cor[i]+11.087*cor[i]*cor[i]-22.348*cor[i]*cor[i]*cor[i]+27.467*cor[i]*cor[i]*cor[i]*cor[i]-18.022*cor[i]*cor[i]*cor[i]*cor[i]*cor[i]+4.8218*cor[i]*cor[i]*cor[i]*cor[i]*cor[i]*cor[i];
      for (int j=i; j <= ntypes; j++){
	  E[i][j] = E[j][i] = 1/((1-pois[i]*pois[i])/E_one[i]+(1-pois[j]*pois[j])/E_one[j]);
	  G[i][j] = G[j][i] = 1/((2-pois[i])/G_one[i]+(2-pois[j])/G_one[j]);
	  if (normaldamp == TSUJI){
	      alpha[i][j] = alpha[j][i] = sqrt(alpha_one[i]*alpha_one[j]);
	  }
	  else if (normaldamp == BRILLIANTOV){
	      gamman[i][j] = gamman[j][i] = sqrt(cor[i]*cor[j]);
	  }
	  muS[i][j] = muS[j][i] = sqrt(muS_one[i]*muS_one[j]);
	  Ecoh[i][j] = Ecoh[j][i] = sqrt(Ecoh_one[i]*Ecoh_one[j]);
	  kR[i][j] = kR[j][i] = sqrt(kR_one[i]*kR_one[j]);
	  etaR[i][j] = etaR[j][i] = sqrt(etaR_one[i]*etaR_one[j]);
	  muR[i][j] = muR[j][i] = sqrt(muR_one[i]*muR_one[j]);
      }
  }
}

/* ---------------------------------------------------------------------- */

double PairGranDMTRolling::single(int i, int j, int itype, int jtype,
				  double rsq,
				  double factor_coul, double factor_lj,
				  double &fforce)
{
  double radi,radj,radsum;
  double r,rinv,rsqinv,delx,dely,delz, nx, ny, nz, R;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3,wr1,wr2,wr3;
  double overlap, a;
  double mi,mj,meff,damp,kn,kt;
  double Fhz,Fdamp,Fdmt,Fne,Fntot,Fscrit;
  double eta_N,eta_T;
  double vtr1,vtr2,vtr3,vrel;
  double fs1,fs2,fs3,fs;
  double shrmag;


  double *radius = atom->radius;
  radi = radius[i];
  radj = radius[j];
  radsum = radi + radj;

  if (rsq >= radsum*radsum) {
      fforce = 0.0;
      svector[0] = svector[1] = svector[2] = svector[3] = 0.0;
      return 0.0;
  }

  r = sqrt(rsq);
  rinv = 1.0/r;
  rsqinv = 1.0/rsq;
  R = radi*radj/radsum;

  // relative translational velocity

  double **v = atom->v;
  vr1 = v[i][0] - v[j][0];
  vr2 = v[i][1] - v[j][1];
  vr3 = v[i][2] - v[j][2];

  // normal component

  double **x = atom->x;
  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  delz = x[i][2] - x[j][2];

  nx = delx*rinv;
  ny = dely*rinv;
  nz = delz*rinv;


  vnnr = vr1*nx + vr2*ny + vr3*nz;
  vn1 = nx*vnnr;
  vn2 = ny*vnnr;
  vn3 = nz*vnnr;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity

  double **omega = atom->omega;
  wr1 = (radi*omega[i][0] + radj*omega[j][0]);
  wr2 = (radi*omega[i][1] + radj*omega[j][1]);
  wr3 = (radi*omega[i][2] + radj*omega[j][2]);

  // meff = effective mass of pair of particles
  // if I or J part of rigid body, use body mass
  // if I or J is frozen, meff is other particle

  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;

  mi = rmass[i];
  mj = rmass[j];
  if (fix_rigid) {
      // NOTE: ensure mass_rigid is current for owned+ghost atoms?
      if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
      if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
  }

  meff = mi*mj / (mi+mj);
  if (mask[i] & freeze_group_bit) meff = mj;
  if (mask[j] & freeze_group_bit) meff = mi;


  // normal force = Hertzian contact + normal velocity damping
  overlap = radsum - r;
  a = sqrt(R*overlap);
  kn = 4.0/3.0*E[itype][jtype]*a;
  Fhz = kn*overlap;

  //Damping (based on Tsuji et al)

  eta_N=alpha[itype][jtype]*sqrt(meff*kn);
  Fdamp = -eta_N*vnnr; //F_nd eq 23 and Zhao eq 19

  //DMT
  Fdmt = -4*MY_PI*Ecoh[itype][jtype]*R;

  Fne = Fhz + Fdmt;
  Fntot = Fne + Fdamp;

  // relative velocities

  vtr1 = vt1 - (nz*wr2-ny*wr3);
  vtr2 = vt2 - (nx*wr3-nz*wr1);
  vtr3 = vt3 - (ny*wr1-nx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // shear history effects
  // neighprev = index of found neigh on previous call
  // search entire jnum list of neighbors of I for neighbor J
  // start from neighprev, since will typically be next neighbor
  // reset neighprev to 0 as necessary

  int jnum = list->numneigh[i];
  int *jlist = list->firstneigh[i];
  double *allshear = fix_history->firstvalue[i];

  for (int jj = 0; jj < jnum; jj++) {
      neighprev++;
      if (neighprev >= jnum) neighprev = 0;
      if (jlist[neighprev] == j) break;
  }

  double *shear = &allshear[size_history*neighprev];
  shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] +
		shear[2]*shear[2]);

  // tangential forces = shear + tangential velocity damping 
  kt=8.0*G[itype][jtype]*a;

  eta_T = eta_N; 
  fs1 = -kt*shear[0] - eta_T*vtr1;
  fs2 = -kt*shear[1] - eta_T*vtr2;
  fs3 = -kt*shear[2] - eta_T*vtr3;

  // rescale frictional displacements and forces if needed

  fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
  Fscrit= muS[itype][jtype] * fabs(Fne);

  if (fs > Fscrit) {
      if (shrmag != 0.0) {
	  fs1 *= Fscrit/fs;
	  fs2 *= Fscrit/fs;
	  fs3 *= Fscrit/fs;
	  fs *= Fscrit/fs;
      } else fs1 = fs2 = fs3 = fs = 0.0;
  }

  // set all forces and return no energy

  fforce = Fntot;

  // set single_extra quantities

  svector[0] = fs1;
  svector[1] = fs2;
  svector[2] = fs3;
  svector[3] = fs;
  svector[4] = vn1;
  svector[5] = vn2;
  svector[6] = vn3;
  svector[7] = vt1;
  svector[8] = vt2;
  svector[9] = vt3;
  return 0.0;
}
