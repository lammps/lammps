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

#include "omp_compat.h"
#include <cmath>
#include "dihedral_helix_omp.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "timer.h"
#include "force.h"
#include "update.h"
#include "math_const.h"
#include "error.h"

#include "suffix.h"
using namespace LAMMPS_NS;
using namespace MathConst;

#define TOLERANCE 0.05
#define SMALL     0.001
#define SMALLER   0.00001

/* ---------------------------------------------------------------------- */

DihedralHelixOMP::DihedralHelixOMP(class LAMMPS *lmp)
  : DihedralHelix(lmp), ThrOMP(lmp,THR_DIHEDRAL)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void DihedralHelixOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = neighbor->ndihedrallist;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, cvatom, thr);

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
void DihedralHelixOMP::eval(int nfrom, int nto, ThrData * const thr)
{

  int i1,i2,i3,i4,n,type;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
  double edihedral,f1[3],f2[3],f3[3],f4[3];
  double sb1,sb2,sb3,rb1,rb3,c0,b1mag2,b1mag,b2mag2;
  double b2mag,b3mag2,b3mag,ctmp,r12c1,c1mag,r12c2;
  double c2mag,sc1,sc2,s1,s12,c,pd,a,a11,a22;
  double a33,a12,a13,a23,sx2,sy2,sz2;
  double s2,cx,cy,cz,cmag,dx,phi,si,siinv,sin2;

  edihedral = 0.0;

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const int5_t * _noalias const dihedrallist = (int5_t *) neighbor->dihedrallist[0];
  const int nlocal = atom->nlocal;

  for (n = nfrom; n < nto; n++) {
    i1 = dihedrallist[n].a;
    i2 = dihedrallist[n].b;
    i3 = dihedrallist[n].c;
    i4 = dihedrallist[n].d;
    type = dihedrallist[n].t;

    // 1st bond

    vb1x = x[i1].x - x[i2].x;
    vb1y = x[i1].y - x[i2].y;
    vb1z = x[i1].z - x[i2].z;

    // 2nd bond

    vb2x = x[i3].x - x[i2].x;
    vb2y = x[i3].y - x[i2].y;
    vb2z = x[i3].z - x[i2].z;

    vb2xm = -vb2x;
    vb2ym = -vb2y;
    vb2zm = -vb2z;

    // 3rd bond

    vb3x = x[i4].x - x[i3].x;
    vb3y = x[i4].y - x[i3].y;
    vb3z = x[i4].z - x[i3].z;

    // c0 calculation

    sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z);
    sb2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z);
    sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z);

    rb1 = sqrt(sb1);
    rb3 = sqrt(sb3);

    c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3;

    // 1st and 2nd angle

    b1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z;
    b1mag = sqrt(b1mag2);
    b2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z;
    b2mag = sqrt(b2mag2);
    b3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z;
    b3mag = sqrt(b3mag2);

    ctmp = vb1x*vb2x + vb1y*vb2y + vb1z*vb2z;
    r12c1 = 1.0 / (b1mag*b2mag);
    c1mag = ctmp * r12c1;

    ctmp = vb2xm*vb3x + vb2ym*vb3y + vb2zm*vb3z;
    r12c2 = 1.0 / (b2mag*b3mag);
    c2mag = ctmp * r12c2;

    // cos and sin of 2 angles and final c

    sin2 = MAX(1.0 - c1mag*c1mag,0.0);
    sc1 = sqrt(sin2);
    if (sc1 < SMALL) sc1 = SMALL;
    sc1 = 1.0/sc1;

    sin2 = MAX(1.0 - c2mag*c2mag,0.0);
    sc2 = sqrt(sin2);
    if (sc2 < SMALL) sc2 = SMALL;
    sc2 = 1.0/sc2;

    s1 = sc1 * sc1;
    s2 = sc2 * sc2;
    s12 = sc1 * sc2;
    c = (c0 + c1mag*c2mag) * s12;

    cx = vb1y*vb2z - vb1z*vb2y;
    cy = vb1z*vb2x - vb1x*vb2z;
    cz = vb1x*vb2y - vb1y*vb2x;
    cmag = sqrt(cx*cx + cy*cy + cz*cz);
    dx = (cx*vb3x + cy*vb3y + cz*vb3z)/cmag/b3mag;

    // error check

    if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) {
      int me = comm->me;

      if (screen) {
        char str[128];
        sprintf(str,"Dihedral problem: %d/%d " BIGINT_FORMAT " "
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

    phi = acos(c);
    if (dx < 0.0) phi *= -1.0;
    si = sin(phi);
    if (fabs(si) < SMALLER) si = SMALLER;
    siinv = 1.0/si;

    pd = -aphi[type] + 3.0*bphi[type]*sin(3.0*phi)*siinv +
      cphi[type]*sin(phi + MY_PI4)*siinv;

    if (EFLAG) edihedral = aphi[type]*(1.0 - c) + bphi[type]*(1.0 + cos(3.0*phi)) +
                 cphi[type]*(1.0 + cos(phi + MY_PI4));

    a = pd;
    c = c * a;
    s12 = s12 * a;
    a11 = c*sb1*s1;
    a22 = -sb2 * (2.0*c0*s12 - c*(s1+s2));
    a33 = c*sb3*s2;
    a12 = -r12c1 * (c1mag*c*s1 + c2mag*s12);
    a13 = -rb1*rb3*s12;
    a23 = r12c2 * (c2mag*c*s2 + c1mag*s12);

    sx2  = a12*vb1x + a22*vb2x + a23*vb3x;
    sy2  = a12*vb1y + a22*vb2y + a23*vb3y;
    sz2  = a12*vb1z + a22*vb2z + a23*vb3z;

    f1[0] = a11*vb1x + a12*vb2x + a13*vb3x;
    f1[1] = a11*vb1y + a12*vb2y + a13*vb3y;
    f1[2] = a11*vb1z + a12*vb2z + a13*vb3z;

    f2[0] = -sx2 - f1[0];
    f2[1] = -sy2 - f1[1];
    f2[2] = -sz2 - f1[2];

    f4[0] = a13*vb1x + a23*vb2x + a33*vb3x;
    f4[1] = a13*vb1y + a23*vb2y + a33*vb3y;
    f4[2] = a13*vb1z + a23*vb2z + a33*vb3z;

    f3[0] = sx2 - f4[0];
    f3[1] = sy2 - f4[1];
    f3[2] = sz2 - f4[2];

    // apply force to each of 4 atoms

    if (NEWTON_BOND || i1 < nlocal) {
      f[i1].x += f1[0];
      f[i1].y += f1[1];
      f[i1].z += f1[2];
    }

    if (NEWTON_BOND || i2 < nlocal) {
      f[i2].x += f2[0];
      f[i2].y += f2[1];
      f[i2].z += f2[2];
    }

    if (NEWTON_BOND || i3 < nlocal) {
      f[i3].x += f3[0];
      f[i3].y += f3[1];
      f[i3].z += f3[2];
    }

    if (NEWTON_BOND || i4 < nlocal) {
      f[i4].x += f4[0];
      f[i4].y += f4[1];
      f[i4].z += f4[2];
    }

    if (EVFLAG)
      ev_tally_thr(this,i1,i2,i3,i4,nlocal,NEWTON_BOND,edihedral,f1,f3,f4,
                   vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,thr);
  }
}
