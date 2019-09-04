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

#include <cmath>
#include "improper_umbrella_omp.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "timer.h"
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
  ev_init(eflag,vflag);

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
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, thr);

    if (inum > 0) {
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
    }
    thr->timer(Timer::BOND);
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

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const int5_t * _noalias const improperlist = (int5_t *) neighbor->improperlist[0];
  const int nlocal = atom->nlocal;

  for (n = nfrom; n < nto; n++) {
    i1 = improperlist[n].a;
    i2 = improperlist[n].b;
    i3 = improperlist[n].c;
    i4 = improperlist[n].d;
    type = improperlist[n].t;

    // 1st bond

    vb1x = x[i2].x - x[i1].x;
    vb1y = x[i2].y - x[i1].y;
    vb1z = x[i2].z - x[i1].z;

    // 2nd bond

    vb2x = x[i3].x - x[i1].x;
    vb2y = x[i3].y - x[i1].y;
    vb2z = x[i3].z - x[i1].z;

    // 3rd bond

    vb3x = x[i4].x - x[i1].x;
    vb3y = x[i4].y - x[i1].y;
    vb3z = x[i4].z - x[i1].z;

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
        sprintf(str,"Improper problem: %d/%d " BIGINT_FORMAT " "
                TAGINT_FORMAT " " TAGINT_FORMAT " "
                TAGINT_FORMAT " " TAGINT_FORMAT,
                me,thr->get_tid(),update->ntimestep,
                atom->tag[i1],atom->tag[i2],atom->tag[i3],atom->tag[i4]);
        error->warning(FLERR,str,0);
        fprintf(screen,"  1st atom: %d %g %g %g\n",
                me,x[i1].x,x[i1].y,x[i1].z);
        fprintf(screen,"  2nd atom: %d %g %g %g\n",
                me,x[i2].x,x[i2].y,x[i2].z);
        fprintf(screen,"  3rd atom: %d %g %g %g\n",
                me,x[i3].x,x[i3].y,x[i3].z);
        fprintf(screen,"  4th atom: %d %g %g %g\n",
                me,x[i4].x,x[i4].y,x[i4].z);
      }
    }

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

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

    f2[0] = (dhay*vb1z - dhaz*vb1y)*rar*a;
    f2[1] = (dhaz*vb1x - dhax*vb1z)*rar*a;
    f2[2] = (dhax*vb1y - dhay*vb1x)*rar*a;

    f3[0] = (-dhay*vb2z + dhaz*vb2y)*rar*a;
    f3[1] = (-dhaz*vb2x + dhax*vb2z)*rar*a;
    f3[2] = (-dhax*vb2y + dhay*vb2x)*rar*a;

    f4[0] = dahx*rhr*a;
    f4[1] = dahy*rhr*a;
    f4[2] = dahz*rhr*a;

    f1[0] = -(f2[0] + f3[0] + f4[0]);
    f1[1] = -(f2[1] + f3[1] + f4[1]);
    f1[2] = -(f2[2] + f3[2] + f4[2]);

    // apply force to each of 4 atoms

    if (NEWTON_BOND || i1 < nlocal) {
      f[i1].x += f1[0];
      f[i1].y += f1[1];
      f[i1].z += f1[2];
    }

    if (NEWTON_BOND || i2 < nlocal) {
      f[i2].x += f3[0];
      f[i2].y += f3[1];
      f[i2].z += f3[2];
    }

    if (NEWTON_BOND || i3 < nlocal) {
      f[i3].x += f2[0];
      f[i3].y += f2[1];
      f[i3].z += f2[2];
    }

    if (NEWTON_BOND || i4 < nlocal) {
      f[i4].x += f4[0];
      f[i4].y += f4[1];
      f[i4].z += f4[2];
    }

    if (EVFLAG) {

      // get correct 4-body geometry for virial tally

      vb1x = x[i1].x - x[i2].x;
      vb1y = x[i1].y - x[i2].y;
      vb1z = x[i1].z - x[i2].z;

      vb2x = x[i3].x - x[i2].x;
      vb2y = x[i3].y - x[i2].y;
      vb2z = x[i3].z - x[i2].z;

      vb3x = x[i4].x - x[i3].x;
      vb3y = x[i4].y - x[i3].y;
      vb3z = x[i4].z - x[i3].z;

      ev_tally_thr(this,i1,i2,i3,i4,nlocal,NEWTON_BOND,eimproper,f1,f2,f4,
                   vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,thr);
    }
  }
}
