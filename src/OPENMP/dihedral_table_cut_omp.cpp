// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "dihedral_table_cut_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "neighbor.h"

#include <cmath>

#include "omp_compat.h"
#include "suffix.h"
using namespace LAMMPS_NS;
using namespace MathConst;

static constexpr double TOLERANCE = 0.05;
static constexpr double SMALL =     0.0000001;

/* ---------------------------------------------------------------------- */

DihedralTableCutOMP::DihedralTableCutOMP(class LAMMPS *lmp)
  : DihedralTableCut(lmp), ThrOMP(lmp,THR_DIHEDRAL)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void DihedralTableCutOMP::compute(int eflag, int vflag)
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
void DihedralTableCutOMP::eval(int nfrom, int nto, ThrData * const thr)
{
  int i1,i2,i3,i4,i,j,k,n,type;
  double edihedral;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
  double fphi,fpphi;
  double r1mag2,r1,r2mag2,r2,r3mag2,r3;
  double sb1,rb1,sb2,rb2,sb3,rb3,c0,r12c1;
  double r12c2,costh12,costh13,costh23,sc1,sc2,s1,s2,c;
  double sinphi,a11,a22,a33,a12,a13,a23,sx1,sx2;
  double sx12,sy1,sy2,sy12,sz1,sz2,sz12;
  double t1,t2,t3,t4;
  double da1,da2;
  double s12,sin2;
  double dcosphidr[4][3],dphidr[4][3],dthetadr[2][4][3];
  double fabcd[4][3];

  edihedral = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const int * const * const dihedrallist = neighbor->dihedrallist;
  const int nlocal = atom->nlocal;

  for (n = nfrom; n < nto; n++) {
    i1 = dihedrallist[n][0];
    i2 = dihedrallist[n][1];
    i3 = dihedrallist[n][2];
    i4 = dihedrallist[n][3];
    type = dihedrallist[n][4];

    // 1st bond

    vb1x = x[i1][0] - x[i2][0];
    vb1y = x[i1][1] - x[i2][1];
    vb1z = x[i1][2] - x[i2][2];

    // 2nd bond

    vb2x = x[i3][0] - x[i2][0];
    vb2y = x[i3][1] - x[i2][1];
    vb2z = x[i3][2] - x[i2][2];

    vb2xm = -vb2x;
    vb2ym = -vb2y;
    vb2zm = -vb2z;

    // 3rd bond

    vb3x = x[i4][0] - x[i3][0];
    vb3y = x[i4][1] - x[i3][1];
    vb3z = x[i4][2] - x[i3][2];

    // distances

    r1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z;
    r1 = sqrt(r1mag2);
    r2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z;
    r2 = sqrt(r2mag2);
    r3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z;
    r3 = sqrt(r3mag2);

    sb1 = 1.0/r1mag2;
    rb1 = 1.0/r1;
    sb2 = 1.0/r2mag2;
    rb2 = 1.0/r2;
    sb3 = 1.0/r3mag2;
    rb3 = 1.0/r3;

    c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3;

    // angles

    r12c1 = rb1*rb2;
    r12c2 = rb2*rb3;
    costh12 = (vb1x*vb2x + vb1y*vb2y + vb1z*vb2z) * r12c1;
    costh13 = c0;
    costh23 = (vb2xm*vb3x + vb2ym*vb3y + vb2zm*vb3z) * r12c2;

    // cos and sin of 2 angles and final c

    sin2 = MAX(1.0 - costh12*costh12,0.0);
    sc1 = sqrt(sin2);
    if (sc1 < SMALL) sc1 = SMALL;
    sc1 = 1.0/sc1;

    sin2 = MAX(1.0 - costh23*costh23,0.0);
    sc2 = sqrt(sin2);
    if (sc2 < SMALL) sc2 = SMALL;
    sc2 = 1.0/sc2;

    s1 = sc1 * sc1;
    s2 = sc2 * sc2;
    s12 = sc1 * sc2;
    c = (c0 + costh12*costh23) * s12;

    // error check

    if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE))
      problem(FLERR, i1, i2, i3, i4);

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    double phil = acos(c);

    sinphi = sqrt(1.0 - c*c);
    sinphi = MAX(sinphi,SMALL);

    // n123 = vb1 x vb2

    double n123x = vb1y*vb2z - vb1z*vb2y;
    double n123y = vb1z*vb2x - vb1x*vb2z;
    double n123z = vb1x*vb2y - vb1y*vb2x;
    double n123_dot_vb3 = n123x*vb3x + n123y*vb3y + n123z*vb3z;
    if (n123_dot_vb3 > 0.0) {
      phil = -phil;
      sinphi = -sinphi;
    }

    a11 = -c*sb1*s1;
    a22 = sb2 * (2.0*costh13*s12 - c*(s1+s2));
    a33 = -c*sb3*s2;
    a12 = r12c1 * (costh12*c*s1 + costh23*s12);
    a13 = rb1*rb3*s12;
    a23 = r12c2 * (-costh23*c*s2 - costh12*s12);

    sx1  = a11*vb1x + a12*vb2x + a13*vb3x;
    sx2  = a12*vb1x + a22*vb2x + a23*vb3x;
    sx12 = a13*vb1x + a23*vb2x + a33*vb3x;
    sy1  = a11*vb1y + a12*vb2y + a13*vb3y;
    sy2  = a12*vb1y + a22*vb2y + a23*vb3y;
    sy12 = a13*vb1y + a23*vb2y + a33*vb3y;
    sz1  = a11*vb1z + a12*vb2z + a13*vb3z;
    sz2  = a12*vb1z + a22*vb2z + a23*vb3z;
    sz12 = a13*vb1z + a23*vb2z + a33*vb3z;

    // set up d(cos(phi))/d(r) and dphi/dr arrays

    dcosphidr[0][0] = -sx1;
    dcosphidr[0][1] = -sy1;
    dcosphidr[0][2] = -sz1;
    dcosphidr[1][0] = sx2 + sx1;
    dcosphidr[1][1] = sy2 + sy1;
    dcosphidr[1][2] = sz2 + sz1;
    dcosphidr[2][0] = sx12 - sx2;
    dcosphidr[2][1] = sy12 - sy2;
    dcosphidr[2][2] = sz12 - sz2;
    dcosphidr[3][0] = -sx12;
    dcosphidr[3][1] = -sy12;
    dcosphidr[3][2] = -sz12;

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
        dphidr[i][j] = -dcosphidr[i][j] / sinphi;


    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
        fabcd[i][j] = 0;
    edihedral = 0;


    // set up d(theta)/d(r) array
    // dthetadr(i,j,k) = angle i, atom j, coordinate k

    for (i = 0; i < 2; i++)
      for (j = 0; j < 4; j++)
        for (k = 0; k < 3; k++)
          dthetadr[i][j][k] = 0.0;

    t1 = costh12 / r1mag2;
    t2 = costh23 / r2mag2;
    t3 = costh12 / r2mag2;
    t4 = costh23 / r3mag2;

    // angle12

    dthetadr[0][0][0] = sc1 * ((t1 * vb1x) - (vb2x * r12c1));
    dthetadr[0][0][1] = sc1 * ((t1 * vb1y) - (vb2y * r12c1));
    dthetadr[0][0][2] = sc1 * ((t1 * vb1z) - (vb2z * r12c1));

    dthetadr[0][1][0] = sc1 * ((-t1 * vb1x) + (vb2x * r12c1) +
                               (-t3 * vb2x) + (vb1x * r12c1));
    dthetadr[0][1][1] = sc1 * ((-t1 * vb1y) + (vb2y * r12c1) +
                               (-t3 * vb2y) + (vb1y * r12c1));
    dthetadr[0][1][2] = sc1 * ((-t1 * vb1z) + (vb2z * r12c1) +
                               (-t3 * vb2z) + (vb1z * r12c1));

    dthetadr[0][2][0] = sc1 * ((t3 * vb2x) - (vb1x * r12c1));
    dthetadr[0][2][1] = sc1 * ((t3 * vb2y) - (vb1y * r12c1));
    dthetadr[0][2][2] = sc1 * ((t3 * vb2z) - (vb1z * r12c1));

    // angle23

    dthetadr[1][1][0] = sc2 * ((t2 * vb2x) + (vb3x * r12c2));
    dthetadr[1][1][1] = sc2 * ((t2 * vb2y) + (vb3y * r12c2));
    dthetadr[1][1][2] = sc2 * ((t2 * vb2z) + (vb3z * r12c2));

    dthetadr[1][2][0] = sc2 * ((-t2 * vb2x) - (vb3x * r12c2) +
                               (t4 * vb3x) + (vb2x * r12c2));
    dthetadr[1][2][1] = sc2 * ((-t2 * vb2y) - (vb3y * r12c2) +
                               (t4 * vb3y) + (vb2y * r12c2));
    dthetadr[1][2][2] = sc2 * ((-t2 * vb2z) - (vb3z * r12c2) +
                               (t4 * vb3z) + (vb2z * r12c2));

    dthetadr[1][3][0] = -sc2 * ((t4 * vb3x) + (vb2x * r12c2));
    dthetadr[1][3][1] = -sc2 * ((t4 * vb3y) + (vb2y * r12c2));
    dthetadr[1][3][2] = -sc2 * ((t4 * vb3z) + (vb2z * r12c2));

    // angle/angle/torsion cutoff

    da1 = acos(costh12) - aat_theta0_1[type];
    da2 = acos(costh23) - aat_theta0_1[type];
    double dtheta = aat_theta0_2[type]-aat_theta0_1[type];

    fphi = 0.0;
    fpphi = 0.0;
    if (phil < 0) phil += MY_2PI;
    uf_lookup(type, phil, fphi, fpphi);

    double gt = aat_k[type];
    double gtt = aat_k[type];
    double gpt = 0;
    double gptt = 0;

    if (acos(costh12) > aat_theta0_1[type]) {
      gt *= 1-da1*da1/dtheta/dtheta;
      gpt = -aat_k[type]*2*da1/dtheta/dtheta;
    }

    if (acos(costh23) > aat_theta0_1[type]) {
      gtt *= 1-da2*da2/dtheta/dtheta;
      gptt = -aat_k[type]*2*da2/dtheta/dtheta;
    }

    if (EFLAG) edihedral = gt*gtt*fphi;

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
        fabcd[i][j] -= gt*gtt*fpphi*dphidr[i][j]
          - gt*gptt*fphi*dthetadr[1][i][j] + gpt*gtt*fphi*dthetadr[0][i][j];

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
      ev_tally_thr(this,i1,i2,i3,i4,nlocal,NEWTON_BOND,edihedral,
                   fabcd[0],fabcd[2],fabcd[3],
                   vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,thr);
  }
}
