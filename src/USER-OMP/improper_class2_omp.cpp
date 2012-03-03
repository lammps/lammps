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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "math.h"
#include "improper_class2_omp.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "error.h"

#include "suffix.h"
using namespace LAMMPS_NS;

#define TOLERANCE 0.05
#define SMALL     0.001

/* ---------------------------------------------------------------------- */

ImproperClass2OMP::ImproperClass2OMP(class LAMMPS *lmp)
  : ImproperClass2(lmp), ThrOMP(lmp,THR_IMPROPER)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void ImproperClass2OMP::compute(int eflag, int vflag)
{

  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else evflag = 0;

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = neighbor->nimproperlist;

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, thr);

    if (evflag) {
      if (eflag) {
	if (force->newton_bond) eval<1,1,1>(ifrom, ito, thr);
	else eval<1,1,0>(ifrom, ito, thr);
      } else {
	if (force->newton_bond) eval<1,0,1>(ifrom, ito, thr);
	else eval<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_bond) eval<0,0,1>(ifrom, ito, thr);
      else eval<0,0,0>(ifrom, ito, thr);
    }

    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <const int EVFLAG, const int EFLAG, const int NEWTON_BOND>
void ImproperClass2OMP::eval(int nfrom, int nto, ThrData * const thr)
{
  int i1,i2,i3,i4,i,j,k,n,type;
  double eimproper;
  double delr[3][3],rmag[3],rinvmag[3],rmag2[3];
  double theta[3],costheta[3],sintheta[3];
  double cossqtheta[3],sinsqtheta[3],invstheta[3];
  double rABxrCB[3],rDBxrAB[3],rCBxrDB[3];
  double ddelr[3][4],dr[3][4][3],dinvr[3][4][3];
  double dthetadr[3][4][3],dinvsth[3][4][3];
  double dinv3r[4][3],dinvs3r[3][4][3];
  double drCBxrDB[3],rCBxdrDB[3],drDBxrAB[3],rDBxdrAB[3];
  double drABxrCB[3],rABxdrCB[3];
  double dot1,dot2,dd[3];
  double fdot[3][4][3],ftmp,invs3r[3],inv3r;
  double t,tt1,tt3,sc1;
  double dotCBDBAB,dotDBABCB,dotABCBDB;
  double chi,deltachi,d2chi,cossin2;
  double drAB[3][4][3],drCB[3][4][3],drDB[3][4][3];
  double dchi[3][4][3],dtotalchi[4][3];
  double schiABCD,chiABCD,schiCBDA,chiCBDA,schiDBAC,chiDBAC;
  double fabcd[4][3];

  eimproper = 0.0;
 
  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const int * const * const improperlist = neighbor->improperlist;
  const int nlocal = atom->nlocal;

  for (i = 0; i < 3; i++)
    for (j = 0; j < 4; j++)
      for (k = 0; k < 3; k++) {
	dthetadr[i][j][k] = 0.0;
	drAB[i][j][k] = 0.0;
	drCB[i][j][k] = 0.0;
	drDB[i][j][k] = 0.0;
      }

  for (n = nfrom; n < nto; n++) {
    i1 = improperlist[n][0];
    i2 = improperlist[n][1];
    i3 = improperlist[n][2];
    i4 = improperlist[n][3];
    type = improperlist[n][4];

    if (k0[type] == 0.0) continue;

    // difference vectors

    delr[0][0] = x[i1][0] - x[i2][0];
    delr[0][1] = x[i1][1] - x[i2][1];
    delr[0][2] = x[i1][2] - x[i2][2];
    domain->minimum_image(delr[0]);

    delr[1][0] = x[i3][0] - x[i2][0];
    delr[1][1] = x[i3][1] - x[i2][1];
    delr[1][2] = x[i3][2] - x[i2][2];
    domain->minimum_image(delr[1]);

    delr[2][0] = x[i4][0] - x[i2][0];
    delr[2][1] = x[i4][1] - x[i2][1];
    delr[2][2] = x[i4][2] - x[i2][2];
    domain->minimum_image(delr[2]);

    // bond lengths and associated values

    for (i = 0; i < 3; i++) {
      rmag2[i] = delr[i][0]*delr[i][0] + delr[i][1]*delr[i][1] + 
	delr[i][2]*delr[i][2];
      rmag[i] = sqrt(rmag2[i]);
      rinvmag[i] = 1.0/rmag[i];
    }

    // angle ABC, CBD, ABD

    costheta[0] = (delr[0][0]*delr[1][0] + delr[0][1]*delr[1][1] +  
		   delr[0][2]*delr[1][2]) / (rmag[0]*rmag[1]);
    costheta[1] = (delr[1][0]*delr[2][0] + delr[1][1]*delr[2][1] + 
		   delr[1][2]*delr[2][2]) / (rmag[1]*rmag[2]);
    costheta[2] = (delr[0][0]*delr[2][0] + delr[0][1]*delr[2][1] + 
		   delr[0][2]*delr[2][2]) / (rmag[0]*rmag[2]);

    // angle error check

    for (i = 0; i < 3; i++) {
      if (costheta[i] == -1.0) {
	int me = comm->me;
	if (screen) {
	  char str[128];
	  sprintf(str,
		  "Improper problem: %d/%d " BIGINT_FORMAT " %d %d %d %d",
		  me, thr->get_tid(),update->ntimestep,
		  atom->tag[i1],atom->tag[i2],atom->tag[i3],atom->tag[i4]);
	  error->warning(FLERR,str,0);
	  fprintf(screen,"  1st atom: %d %g %g %g\n",
		  me,x[i1][0],x[i1][1],x[i1][2]);
	  fprintf(screen,"  2nd atom: %d %g %g %g\n",
		  me,x[i2][0],x[i2][1],x[i2][2]);
	  fprintf(screen,"  3rd atom: %d %g %g %g\n",
		  me,x[i3][0],x[i3][1],x[i3][2]);
	  fprintf(screen,"  4th atom: %d %g %g %g\n",
		  me,x[i4][0],x[i4][1],x[i4][2]);
	}
      }
    }

    for (i = 0; i < 3; i++) {
      if (costheta[i] > 1.0)  costheta[i] = 1.0;
      if (costheta[i] < -1.0) costheta[i] = -1.0;
      theta[i] = acos(costheta[i]);
      cossqtheta[i] = costheta[i]*costheta[i];
      sintheta[i] = sin(theta[i]);
      invstheta[i] = 1.0/sintheta[i];
      sinsqtheta[i] = sintheta[i]*sintheta[i];
    }

    // cross & dot products

    cross(delr[0],delr[1],rABxrCB);
    cross(delr[2],delr[0],rDBxrAB);
    cross(delr[1],delr[2],rCBxrDB);

    dotCBDBAB = dot(rCBxrDB,delr[0]);
    dotDBABCB = dot(rDBxrAB,delr[1]);
    dotABCBDB = dot(rABxrCB,delr[2]);

    t = rmag[0] * rmag[1] * rmag[2];
    inv3r = 1.0/t;
    invs3r[0] = invstheta[1] * inv3r;
    invs3r[1] = invstheta[2] * inv3r;
    invs3r[2] = invstheta[0] * inv3r;

    // chi ABCD, CBDA, DBAC
    // final chi is average of three

    schiABCD = dotCBDBAB * invs3r[0];
    chiABCD = asin(schiABCD);
    schiCBDA = dotDBABCB * invs3r[1];
    chiCBDA = asin(schiCBDA);
    schiDBAC = dotABCBDB * invs3r[2];
    chiDBAC = asin(schiDBAC);

    chi = (chiABCD + chiCBDA + chiDBAC) / 3.0;
    deltachi = chi - chi0[type];
    d2chi = deltachi * deltachi;

    // energy

    if (EFLAG) eimproper = k0[type]*d2chi;

    // forces
    // define d(delr)
    // i = bond AB/CB/DB, j = atom A/B/C/D

    ddelr[0][0] = 1.0;
    ddelr[0][1] = -1.0;
    ddelr[0][2] = 0.0;
    ddelr[0][3] = 0.0;
    ddelr[1][0] = 0.0;
    ddelr[1][1] = -1.0;
    ddelr[1][2] = 1.0;
    ddelr[1][3] = 0.0;
    ddelr[2][0] = 0.0;
    ddelr[2][1] = -1.0;
    ddelr[2][2] = 0.0;
    ddelr[2][3] = 1.0;

    // compute d(|r|)/dr and d(1/|r|)/dr for each direction, bond and atom
    // define d(r) for each r
    // i = bond AB/CB/DB, j = atom A/B/C/D, k = X/Y/Z

    for (i = 0; i < 3; i++)
      for (j = 0; j < 4; j++)
	for (k = 0; k < 3; k++) {
	  dr[i][j][k] = delr[i][k] * ddelr[i][j] / rmag[i];
	  dinvr[i][j][k] = -dr[i][j][k] / rmag2[i];
	}

    // compute d(1 / (|r_AB| * |r_CB| * |r_DB|) / dr
    // i = atom A/B/C/D, j = X/Y/Z

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
	dinv3r[i][j] = rinvmag[1] * (rinvmag[2] * dinvr[0][i][j] +
				     rinvmag[0] * dinvr[2][i][j]) +
	  rinvmag[2] * rinvmag[0] * dinvr[1][i][j];

    // compute d(theta)/d(r) for 3 angles
    // angleABC

    tt1 = costheta[0] / rmag2[0];
    tt3 = costheta[0] / rmag2[1];
    sc1 = 1.0 / sqrt(1.0 - cossqtheta[0]);

    dthetadr[0][0][0] = sc1 * ((tt1 * delr[0][0]) - 
			       (delr[1][0] * rinvmag[0] * rinvmag[1]));
    dthetadr[0][0][1] = sc1 * ((tt1 * delr[0][1]) - 
			       (delr[1][1] * rinvmag[0] * rinvmag[1]));
    dthetadr[0][0][2] = sc1 * ((tt1 * delr[0][2]) - 
			       (delr[1][2] * rinvmag[0] * rinvmag[1]));
    dthetadr[0][1][0] = -sc1 * ((tt1 * delr[0][0]) - 
				(delr[1][0] * rinvmag[0] * rinvmag[1]) +
				(tt3 * delr[1][0]) - 
				(delr[0][0] * rinvmag[0] * rinvmag[1]));
    dthetadr[0][1][1] = -sc1 * ((tt1 * delr[0][1]) - 
				(delr[1][1] * rinvmag[0] * rinvmag[1]) +
				(tt3 * delr[1][1]) - 
				(delr[0][1] * rinvmag[0] * rinvmag[1]));
    dthetadr[0][1][2] = -sc1 * ((tt1 * delr[0][2]) - 
				(delr[1][2] * rinvmag[0] * rinvmag[1]) +
				(tt3 * delr[1][2]) - 
				(delr[0][2] * rinvmag[0] * rinvmag[1]));
    dthetadr[0][2][0] = sc1 * ((tt3 * delr[1][0]) - 
			       (delr[0][0] * rinvmag[0] * rinvmag[1]));
    dthetadr[0][2][1] = sc1 * ((tt3 * delr[1][1]) - 
			       (delr[0][1] * rinvmag[0] * rinvmag[1]));
    dthetadr[0][2][2] = sc1 * ((tt3 * delr[1][2]) - 
			       (delr[0][2] * rinvmag[0] * rinvmag[1]));

    // angleCBD

    tt1 = costheta[1] / rmag2[1];
    tt3 = costheta[1] / rmag2[2];
    sc1 = 1.0 / sqrt(1.0 - cossqtheta[1]);

    dthetadr[1][2][0] = sc1 * ((tt1 * delr[1][0]) - 
			       (delr[2][0] * rinvmag[1] * rinvmag[2]));
    dthetadr[1][2][1] = sc1 * ((tt1 * delr[1][1]) - 
			       (delr[2][1] * rinvmag[1] * rinvmag[2]));
    dthetadr[1][2][2] = sc1 * ((tt1 * delr[1][2]) - 
			       (delr[2][2] * rinvmag[1] * rinvmag[2]));
    dthetadr[1][1][0] = -sc1 * ((tt1 * delr[1][0]) - 
				(delr[2][0] * rinvmag[1] * rinvmag[2]) +
				(tt3 * delr[2][0]) - 
				(delr[1][0] * rinvmag[2] * rinvmag[1]));
    dthetadr[1][1][1] = -sc1 * ((tt1 * delr[1][1]) - 
				(delr[2][1] * rinvmag[1] * rinvmag[2]) +
				(tt3 * delr[2][1]) - 
				(delr[1][1] * rinvmag[2] * rinvmag[1]));
    dthetadr[1][1][2] = -sc1 * ((tt1 * delr[1][2]) - 
				(delr[2][2] * rinvmag[1] * rinvmag[2]) +
				(tt3 * delr[2][2]) - 
				(delr[1][2] * rinvmag[2] * rinvmag[1]));
    dthetadr[1][3][0] = sc1 * ((tt3 * delr[2][0]) - 
			       (delr[1][0] * rinvmag[2] * rinvmag[1]));
    dthetadr[1][3][1] = sc1 * ((tt3 * delr[2][1]) - 
			       (delr[1][1] * rinvmag[2] * rinvmag[1]));
    dthetadr[1][3][2] = sc1 * ((tt3 * delr[2][2]) - 
			       (delr[1][2] * rinvmag[2] * rinvmag[1]));

    // angleABD

    tt1 = costheta[2] / rmag2[0];
    tt3 = costheta[2] / rmag2[2];
    sc1 = 1.0 / sqrt(1.0 - cossqtheta[2]);

    dthetadr[2][0][0] = sc1 * ((tt1 * delr[0][0]) - 
			       (delr[2][0] * rinvmag[0] * rinvmag[2]));
    dthetadr[2][0][1] = sc1 * ((tt1 * delr[0][1]) - 
			       (delr[2][1] * rinvmag[0] * rinvmag[2]));
    dthetadr[2][0][2] = sc1 * ((tt1 * delr[0][2]) - 
			       (delr[2][2] * rinvmag[0] * rinvmag[2]));
    dthetadr[2][1][0] = -sc1 * ((tt1 * delr[0][0]) - 
				(delr[2][0] * rinvmag[0] * rinvmag[2]) +
				(tt3 * delr[2][0]) - 
				(delr[0][0] * rinvmag[2] * rinvmag[0]));
    dthetadr[2][1][1] = -sc1 * ((tt1 * delr[0][1]) - 
				(delr[2][1] * rinvmag[0] * rinvmag[2]) +
				(tt3 * delr[2][1]) - 
				(delr[0][1] * rinvmag[2] * rinvmag[0]));
    dthetadr[2][1][2] = -sc1 * ((tt1 * delr[0][2]) - 
				(delr[2][2] * rinvmag[0] * rinvmag[2]) +
				(tt3 * delr[2][2]) - 
				(delr[0][2] * rinvmag[2] * rinvmag[0]));
    dthetadr[2][3][0] = sc1 * ((tt3 * delr[2][0]) - 
			       (delr[0][0] * rinvmag[2] * rinvmag[0]));
    dthetadr[2][3][1] = sc1 * ((tt3 * delr[2][1]) - 
			       (delr[0][1] * rinvmag[2] * rinvmag[0]));
    dthetadr[2][3][2] = sc1 * ((tt3 * delr[2][2]) - 
			       (delr[0][2] * rinvmag[2] * rinvmag[0]));

    // compute d( 1 / sin(theta))/dr
    // i = angle, j = atom, k = direction

    for (i = 0; i < 3; i++) {
      cossin2 = -costheta[i] / sinsqtheta[i];
      for (j = 0; j < 4; j++)
	for (k = 0; k < 3; k++)
	  dinvsth[i][j][k] = cossin2 * dthetadr[i][j][k];
    }

    // compute d(1 / sin(theta) * |r_AB| * |r_CB| * |r_DB|)/dr
    // i = angle, j = atom

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++) {
	dinvs3r[0][i][j] = (invstheta[1] * dinv3r[i][j]) +
	  (inv3r * dinvsth[1][i][j]);
	dinvs3r[1][i][j] = (invstheta[2] * dinv3r[i][j]) +
	  (inv3r * dinvsth[2][i][j]);
	dinvs3r[2][i][j] = (invstheta[0] * dinv3r[i][j]) +
	  (inv3r * dinvsth[0][i][j]);
      }

    // drCB(i,j,k), etc
    // i = vector X'/Y'/Z', j = atom A/B/C/D, k = direction X/Y/Z

    for (i = 0; i < 3; i++) {
      drCB[i][1][i] = -1.0;
      drAB[i][1][i] = -1.0;
      drDB[i][1][i] = -1.0;
      drDB[i][3][i] = 1.0;
      drCB[i][2][i] = 1.0;
      drAB[i][0][i] = 1.0;
    }

    // d((r_CB x r_DB) dot r_AB)
    // r_CB x d(r_DB)
    // d(r_CB) x r_DB
    // (r_CB x d(r_DB)) + (d(r_CB) x r_DB)
    // (r_CB x d(r_DB)) + (d(r_CB) x r_DB) dot r_AB
    // d(r_AB) dot (r_CB x r_DB)

    for (i = 0; i < 3; i++)
      for (j = 0; j < 4; j++) {
	cross(delr[1],drDB[i][j],rCBxdrDB);
	cross(drCB[i][j],delr[2],drCBxrDB);
	for (k = 0; k < 3; k++) dd[k] = rCBxdrDB[k] + drCBxrDB[k];
	dot1 = dot(dd,delr[0]);
	dot2 = dot(rCBxrDB,drAB[i][j]);
	fdot[0][j][i] = dot1 + dot2;
      }

    // d((r_DB x r_AB) dot r_CB)
    // r_DB x d(r_AB)
    // d(r_DB) x r_AB
    // (r_DB x d(r_AB)) + (d(r_DB) x r_AB)
    // (r_DB x d(r_AB)) + (d(r_DB) x r_AB) dot r_CB
    // d(r_CB) dot (r_DB x r_AB)

    for (i = 0; i < 3; i++)
      for (j = 0; j < 4; j++) {
	cross(delr[2],drAB[i][j],rDBxdrAB);
	cross(drDB[i][j],delr[0],drDBxrAB);
	for (k = 0; k < 3; k++) dd[k] = rDBxdrAB[k] + drDBxrAB[k];
	dot1 = dot(dd,delr[1]);
	dot2 = dot(rDBxrAB,drCB[i][j]);
	fdot[1][j][i] = dot1 + dot2;
      }

    // d((r_AB x r_CB) dot r_DB)
    // r_AB x d(r_CB)
    // d(r_AB) x r_CB
    // (r_AB x d(r_CB)) + (d(r_AB) x r_CB)
    // (r_AB x d(r_CB)) + (d(r_AB) x r_CB) dot r_DB
    // d(r_DB) dot (r_AB x r_CB)

    for (i = 0; i < 3; i++)
      for (j = 0; j < 4; j++) {
	cross(delr[0],drCB[i][j],rABxdrCB);
	cross(drAB[i][j],delr[1],drABxrCB);
	for (k = 0; k < 3; k++) dd[k] = rABxdrCB[k] + drABxrCB[k];
	dot1 = dot(dd,delr[2]);
	dot2 = dot(rABxrCB,drDB[i][j]);
	fdot[2][j][i] = dot1 + dot2;
      }

    // force on each atom

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++) {
	ftmp = (fdot[0][i][j] * invs3r[0]) + 
	  (dinvs3r[0][i][j] * dotCBDBAB);
	dchi[0][i][j] = ftmp / cos(chiABCD);
	ftmp = (fdot[1][i][j] * invs3r[1]) + 
	  (dinvs3r[1][i][j] * dotDBABCB);
	dchi[1][i][j] = ftmp / cos(chiCBDA);
	ftmp = (fdot[2][i][j] * invs3r[2]) + 
	  (dinvs3r[2][i][j] * dotABCBDB);
	dchi[2][i][j] = ftmp / cos(chiDBAC);
	dtotalchi[i][j] = (dchi[0][i][j]+dchi[1][i][j]+dchi[2][i][j]) / 3.0;
      }

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
	fabcd[i][j] = -2.0*k0[type] * deltachi*dtotalchi[i][j];

    // apply force to each of 4 atoms

    if (NEWTON_BOND || i1 < nlocal) {
      f[i1][0] += fabcd[0][0];
      f[i1][1] += fabcd[0][1];
      f[i1][2] += fabcd[0][2];
    }

    if (NEWTON_BOND || i2 < nlocal) {
      f[i2][0] += fabcd[1][0];
      f[i2][1] += fabcd[1][1];
      f[i2][2] += fabcd[1][2];
    }

    if (NEWTON_BOND || i3 < nlocal) {
      f[i3][0] += fabcd[2][0];
      f[i3][1] += fabcd[2][1];
      f[i3][2] += fabcd[2][2];
    }

    if (NEWTON_BOND || i4 < nlocal) {
      f[i4][0] += fabcd[3][0];
      f[i4][1] += fabcd[3][1];
      f[i4][2] += fabcd[3][2];
    }

    if (EVFLAG)
      ev_tally_thr(this,i1,i2,i3,i4,nlocal,NEWTON_BOND,eimproper,
		   fabcd[0],fabcd[2],fabcd[3],
		   delr[0][0],delr[0][1],delr[0][2],
		   delr[1][0],delr[1][1],delr[1][2],
		   delr[2][0]-delr[1][0],delr[2][1]-delr[1][1],
		   delr[2][2]-delr[1][2],thr);
  }

  // compute angle-angle interactions
  angleangle_thr<EVFLAG, EFLAG, NEWTON_BOND>(nfrom,nto,thr);
}

/* ----------------------------------------------------------------------
   angle-angle interactions within improper
------------------------------------------------------------------------- */
template <int EVFLAG, int EFLAG, int NEWTON_BOND>
void ImproperClass2OMP::angleangle_thr(int nfrom, int nto, ThrData * const thr)
{
  int i1,i2,i3,i4,i,j,k,n,type;
  double eimproper;
  double delxAB,delyAB,delzAB,rABmag2,rAB;
  double delxBC,delyBC,delzBC,rBCmag2,rBC;
  double delxBD,delyBD,delzBD,rBDmag2,rBD;
  double costhABC,thetaABC,costhABD;
  double thetaABD,costhCBD,thetaCBD,dthABC,dthCBD,dthABD;
  double sc1,t1,t3,r12;
  double dthetadr[3][4][3],fabcd[4][3];

  eimproper = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const int * const * const improperlist = neighbor->improperlist;
  const int nlocal = atom->nlocal;

  for (n = nfrom; n < nto; n++) {
    i1 = improperlist[n][0];
    i2 = improperlist[n][1];
    i3 = improperlist[n][2];
    i4 = improperlist[n][3];
    type = improperlist[n][4];

    // difference vectors

    delxAB = x[i1][0] - x[i2][0];
    delyAB = x[i1][1] - x[i2][1];
    delzAB = x[i1][2] - x[i2][2];
    domain->minimum_image(delxAB,delyAB,delzAB);

    delxBC = x[i3][0] - x[i2][0];
    delyBC = x[i3][1] - x[i2][1];
    delzBC = x[i3][2] - x[i2][2];
    domain->minimum_image(delxBC,delyBC,delzBC);

    delxBD = x[i4][0] - x[i2][0];
    delyBD = x[i4][1] - x[i2][1];
    delzBD = x[i4][2] - x[i2][2];
    domain->minimum_image(delxBD,delyBD,delzBD);

    // bond lengths

    rABmag2 = delxAB*delxAB + delyAB*delyAB + delzAB*delzAB;
    rAB = sqrt(rABmag2);
    rBCmag2 = delxBC*delxBC + delyBC*delyBC + delzBC*delzBC;
    rBC = sqrt(rBCmag2);
    rBDmag2 = delxBD*delxBD + delyBD*delyBD + delzBD*delzBD;
    rBD = sqrt(rBDmag2);
        
    // angle ABC, ABD, CBD

    costhABC = (delxAB*delxBC + delyAB*delyBC + delzAB*delzBC) / (rAB * rBC);
    if (costhABC > 1.0)  costhABC = 1.0;
    if (costhABC < -1.0) costhABC = -1.0;
    thetaABC = acos(costhABC);

    costhABD = (delxAB*delxBD + delyAB*delyBD + delzAB*delzBD) / (rAB * rBD);
    if (costhABD > 1.0)  costhABD = 1.0;
    if (costhABD < -1.0) costhABD = -1.0;
    thetaABD = acos(costhABD);

    costhCBD = (delxBC*delxBD + delyBC*delyBD + delzBC*delzBD) /(rBC * rBD);
    if (costhCBD > 1.0)  costhCBD = 1.0;
    if (costhCBD < -1.0) costhCBD = -1.0;
    thetaCBD = acos(costhCBD);

    dthABC = thetaABC - aa_theta0_1[type];
    dthABD = thetaABD - aa_theta0_2[type];
    dthCBD = thetaCBD - aa_theta0_3[type];

    // energy

    if (EFLAG) eimproper = aa_k2[type] * dthABC * dthABD + 
		 aa_k1[type] * dthABC * dthCBD +
		 aa_k3[type] * dthABD * dthCBD;

    // d(theta)/d(r) array
    // angle i, atom j, coordinate k

    for (i = 0; i < 3; i++)
      for (j = 0; j < 4; j++)
	for (k = 0; k < 3; k++)
	  dthetadr[i][j][k] = 0.0;

    // angle ABC

    sc1 = sqrt(1.0/(1.0 - costhABC*costhABC));
    t1 = costhABC / rABmag2;
    t3 = costhABC / rBCmag2;
    r12 = 1.0 / (rAB * rBC);

    dthetadr[0][0][0] = sc1 * ((t1 * delxAB) - (delxBC * r12));
    dthetadr[0][0][1] = sc1 * ((t1 * delyAB) - (delyBC * r12));
    dthetadr[0][0][2] = sc1 * ((t1 * delzAB) - (delzBC * r12));
    dthetadr[0][1][0] = sc1 * ((-t1 * delxAB) + (delxBC * r12) +
			       (-t3 * delxBC) + (delxAB * r12));
    dthetadr[0][1][1] = sc1 * ((-t1 * delyAB) + (delyBC * r12) +
			       (-t3 * delyBC) + (delyAB * r12));
    dthetadr[0][1][2] = sc1 * ((-t1 * delzAB) + (delzBC * r12) +
			       (-t3 * delzBC) + (delzAB * r12));
    dthetadr[0][2][0] = sc1 * ((t3 * delxBC) - (delxAB * r12));
    dthetadr[0][2][1] = sc1 * ((t3 * delyBC) - (delyAB * r12));
    dthetadr[0][2][2] = sc1 * ((t3 * delzBC) - (delzAB * r12));

    // angle CBD

    sc1 = sqrt(1.0/(1.0 - costhCBD*costhCBD));
    t1 = costhCBD / rBCmag2;
    t3 = costhCBD / rBDmag2;
    r12 = 1.0 / (rBC * rBD);

    dthetadr[1][2][0] = sc1 * ((t1 * delxBC) - (delxBD * r12));
    dthetadr[1][2][1] = sc1 * ((t1 * delyBC) - (delyBD * r12));
    dthetadr[1][2][2] = sc1 * ((t1 * delzBC) - (delzBD * r12));
    dthetadr[1][1][0] = sc1 * ((-t1 * delxBC) + (delxBD * r12) +
			       (-t3 * delxBD) + (delxBC * r12));
    dthetadr[1][1][1] = sc1 * ((-t1 * delyBC) + (delyBD * r12) +
			       (-t3 * delyBD) + (delyBC * r12));
    dthetadr[1][1][2] = sc1 * ((-t1 * delzBC) + (delzBD * r12) +
			       (-t3 * delzBD) + (delzBC * r12));
    dthetadr[1][3][0] = sc1 * ((t3 * delxBD) - (delxBC * r12));
    dthetadr[1][3][1] = sc1 * ((t3 * delyBD) - (delyBC * r12));
    dthetadr[1][3][2] = sc1 * ((t3 * delzBD) - (delzBC * r12));

    // angle ABD

    sc1 = sqrt(1.0/(1.0 - costhABD*costhABD));
    t1 = costhABD / rABmag2;
    t3 = costhABD / rBDmag2;
    r12 = 1.0 / (rAB * rBD);

    dthetadr[2][0][0] = sc1 * ((t1 * delxAB) - (delxBD * r12));
    dthetadr[2][0][1] = sc1 * ((t1 * delyAB) - (delyBD * r12));
    dthetadr[2][0][2] = sc1 * ((t1 * delzAB) - (delzBD * r12));
    dthetadr[2][1][0] = sc1 * ((-t1 * delxAB) + (delxBD * r12) +
			       (-t3 * delxBD) + (delxAB * r12));
    dthetadr[2][1][1] = sc1 * ((-t1 * delyAB) + (delyBD * r12) +
			       (-t3 * delyBD) + (delyAB * r12));
    dthetadr[2][1][2] = sc1 * ((-t1 * delzAB) + (delzBD * r12) +
			       (-t3 * delzBD) + (delzAB * r12));
    dthetadr[2][3][0] = sc1 * ((t3 * delxBD) - (delxAB * r12));
    dthetadr[2][3][1] = sc1 * ((t3 * delyBD) - (delyAB * r12));
    dthetadr[2][3][2] = sc1 * ((t3 * delzBD) - (delzAB * r12));

    // angleangle forces

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
	fabcd[i][j] = - 
	  ((aa_k1[type] * 
	    (dthABC*dthetadr[1][i][j] + dthCBD*dthetadr[0][i][j])) +
	   (aa_k2[type] * 
	    (dthABC*dthetadr[2][i][j] + dthABD*dthetadr[0][i][j])) +
	   (aa_k3[type] *
	    (dthABD*dthetadr[1][i][j] + dthCBD*dthetadr[2][i][j])));

    // apply force to each of 4 atoms

    if (NEWTON_BOND || i1 < nlocal) {
      f[i1][0] += fabcd[0][0];
      f[i1][1] += fabcd[0][1];
      f[i1][2] += fabcd[0][2];
    }

    if (NEWTON_BOND || i2 < nlocal) {
      f[i2][0] += fabcd[1][0];
      f[i2][1] += fabcd[1][1];
      f[i2][2] += fabcd[1][2];
    }

    if (NEWTON_BOND || i3 < nlocal) {
      f[i3][0] += fabcd[2][0];
      f[i3][1] += fabcd[2][1];
      f[i3][2] += fabcd[2][2];
    }

    if (NEWTON_BOND || i4 < nlocal) {
      f[i4][0] += fabcd[3][0];
      f[i4][1] += fabcd[3][1];
      f[i4][2] += fabcd[3][2];
    }

    if (EVFLAG)
      ev_tally_thr(this,i1,i2,i3,i4,nlocal,NEWTON_BOND,eimproper,
		   fabcd[0],fabcd[2],fabcd[3],delxAB,delyAB,delzAB,
		   delxBC,delyBC,delzBC,delxBD,delyBD,delzBD,thr);
  }
}
