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
#include "improper_umbrella_omp.h"
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

ImproperUmbrellaOMP::ImproperUmbrellaOMP(class LAMMPS *lmp)
  : ImproperUmbrella(lmp), ThrOMP(lmp,THR_IMPROPER)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void ImproperUmbrellaOMP::compute(int eflag, int vflag)
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

template <int EVFLAG, int EFLAG, int NEWTON_BOND>
void ImproperUmbrellaOMP::eval(int nfrom, int nto, ThrData * const thr)
{
  int i1,i2,i3,i4,n,type;
  double eimproper,f1[3],f2[3],f3[3],f4[3];
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z;
  double domega,c,a,s,projhfg,dhax,dhay,dhaz,dahx,dahy,dahz,cotphi;
  double ax,ay,az,ra2,rh2,ra,rh,rar,rhr,arx,ary,arz,hrx,hry,hrz;

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

    // 1st bond

    vb1x = x[i2][0] - x[i1][0];
    vb1y = x[i2][1] - x[i1][1];
    vb1z = x[i2][2] - x[i1][2];
    domain->minimum_image(vb1x,vb1y,vb1z);

    // 2nd bond

    vb2x = x[i3][0] - x[i1][0];
    vb2y = x[i3][1] - x[i1][1];
    vb2z = x[i3][2] - x[i1][2];
    domain->minimum_image(vb2x,vb2y,vb2z);

    // 3rd bond

    vb3x = x[i4][0] - x[i1][0];
    vb3y = x[i4][1] - x[i1][1];
    vb3z = x[i4][2] - x[i1][2];
    domain->minimum_image(vb3x,vb3y,vb3z);

    // c0 calculation
    // A = vb1 X vb2 is perpendicular to IJK plane

    ax = vb1y*vb2z-vb1z*vb2y;
    ay = vb1z*vb2x-vb1x*vb2z;
    az = vb1x*vb2y-vb1y*vb2x;
    ra2 = ax*ax+ay*ay+az*az;
    rh2 = vb3x*vb3x+vb3y*vb3y+vb3z*vb3z;
    ra = sqrt(ra2);
    rh = sqrt(rh2);
    if (ra < SMALL) ra = SMALL;
    if (rh < SMALL) rh = SMALL;

    rar = 1/ra;
    rhr = 1/rh;
    arx = ax*rar;
    ary = ay*rar;
    arz = az*rar;
    hrx = vb3x*rhr;
    hry = vb3y*rhr;
    hrz = vb3z*rhr;

    c = arx*hrx+ary*hry+arz*hrz;

    // error check

    if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) {
      int me = comm->me;

      if (screen) {
	char str[128];
	sprintf(str,"Improper problem: %d/%d " BIGINT_FORMAT " %d %d %d %d",
		me,thr->get_tid(),update->ntimestep,
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
    
    if (c > 1.0) s = 1.0;
    if (c < -1.0) s = -1.0;

    s = sqrt(1.0 - c*c);
    if (s < SMALL) s = SMALL;
    cotphi = c/s;

    projhfg = (vb3x*vb1x+vb3y*vb1y+vb3z*vb1z) / 
      sqrt(vb1x*vb1x+vb1y*vb1y+vb1z*vb1z); 
    projhfg += (vb3x*vb2x+vb3y*vb2y+vb3z*vb2z) / 
      sqrt(vb2x*vb2x+vb2y*vb2y+vb2z*vb2z);
    if (projhfg > 0.0) {
      s *= -1.0;
      cotphi *= -1.0;
    }
	
    //  force and energy
    // if w0 = 0: E = k * (1 - cos w)
    // if w0 != 0: E = 0.5 * C (cos w - cos w0)^2, C = k/(sin(w0)^2

    if (w0[type] == 0.0) {
      if (EFLAG) eimproper = kw[type] * (1.0-s);
      a = -kw[type];
    } else {
      domega = s - cos(w0[type]);
      a = 0.5 * C[type] * domega;
      if (EFLAG) eimproper = a * domega;
      a *= 2.0;
    }

    // dhax = diffrence between H and A in X direction, etc

    a = a*cotphi;
    dhax = hrx-c*arx;
    dhay = hry-c*ary;
    dhaz = hrz-c*arz;

    dahx = arx-c*hrx;
    dahy = ary-c*hry;
    dahz = arz-c*hrz;

    f2[0] = (dhay*vb1z - dhaz*vb1y)*rar;
    f2[1] = (dhaz*vb1x - dhax*vb1z)*rar;
    f2[2] = (dhax*vb1y - dhay*vb1x)*rar;

    f3[0] = (-dhay*vb2z + dhaz*vb2y)*rar;
    f3[1] = (-dhaz*vb2x + dhax*vb2z)*rar;
    f3[2] = (-dhax*vb2y + dhay*vb2x)*rar;

    f4[0] = dahx*rhr;
    f4[1] = dahy*rhr;
    f4[2] = dahz*rhr;

    f1[0] = -(f2[0] + f3[0] + f4[0]);
    f1[1] = -(f2[1] + f3[1] + f4[1]);
    f1[2] = -(f2[2] + f3[2] + f4[2]);

    // apply force to each of 4 atoms

    if (NEWTON_BOND || i1 < nlocal) {
      f[i1][0] += f1[0]*a;
      f[i1][1] += f1[1]*a;
      f[i1][2] += f1[2]*a;
    }

    if (NEWTON_BOND || i2 < nlocal) {
      f[i2][0] += f3[0]*a;
      f[i2][1] += f3[1]*a;
      f[i2][2] += f3[2]*a;
    }

    if (NEWTON_BOND || i3 < nlocal) {
      f[i3][0] += f2[0]*a;
      f[i3][1] += f2[1]*a;
      f[i3][2] += f2[2]*a;
    }

    if (NEWTON_BOND || i4 < nlocal) {
      f[i4][0] += f4[0]*a;
      f[i4][1] += f4[1]*a;
      f[i4][2] += f4[2]*a;
    }

    if (EVFLAG)
      ev_tally_thr(this,i1,i2,i3,i4,nlocal,NEWTON_BOND,eimproper,f1,f3,f4,
		   vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,thr);
  }
}
