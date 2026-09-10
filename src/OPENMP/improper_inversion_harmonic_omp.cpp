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
   Contributing authors: Hendrik Heenen (Technical University of Munich)
                         and Rochus Schmid (Ruhr-Universitaet Bochum)
                         (OMP threading by Axel Kohlmeyer, Temple U)
------------------------------------------------------------------------- */

#include "omp_compat.h"
#include "improper_inversion_harmonic_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"

#include <cmath>

#include "suffix.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ImproperInversionHarmonicOMP::ImproperInversionHarmonicOMP(class LAMMPS *lmp)
  : ImproperInversionHarmonic(lmp), ThrOMP(lmp,THR_IMPROPER)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void ImproperInversionHarmonicOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = neighbor->nimproperlist;

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
void ImproperInversionHarmonicOMP::eval(int nfrom, int nto, ThrData * const thr)
{
  int i1,i2,i3,i4,n,type;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z;
  double rrvb1,rrvb2,rrvb3,rr2vb1,rr2vb2,rr2vb3;

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const int5_t * _noalias const improperlist = (int5_t *) neighbor->improperlist[0];

  for (n = nfrom; n < nto; n++) {
    i1 = improperlist[n].a;
    i2 = improperlist[n].b;
    i3 = improperlist[n].c;
    i4 = improperlist[n].d;
    type = improperlist[n].t;

    // 1st bond - IJ

    vb1x = x[i2].x - x[i1].x;
    vb1y = x[i2].y - x[i1].y;
    vb1z = x[i2].z - x[i1].z;
    rrvb1 = 1.0/sqrt(vb1x*vb1x+vb1y*vb1y+vb1z*vb1z);
    rr2vb1 = rrvb1*rrvb1;

    // 2nd bond - IK

    vb2x = x[i3].x - x[i1].x;
    vb2y = x[i3].y - x[i1].y;
    vb2z = x[i3].z - x[i1].z;
    rrvb2 = 1.0/sqrt(vb2x*vb2x+vb2y*vb2y+vb2z*vb2z);
    rr2vb2 = rrvb2*rrvb2;

    // 3rd bond - IL

    vb3x = x[i4].x - x[i1].x;
    vb3y = x[i4].y - x[i1].y;
    vb3z = x[i4].z - x[i1].z;
    rrvb3 = 1.0/sqrt(vb3x*vb3x+vb3y*vb3y+vb3z*vb3z);
    rr2vb3 = rrvb3*rrvb3;

    // compute all three inversion angles
    invang<EVFLAG,EFLAG,NEWTON_BOND>(i1,i2,i3,i4, type,
           vb3x, vb3y, vb3z, rrvb3, rr2vb3,
           vb2x, vb2y, vb2z, rrvb2, rr2vb2,
           vb1x, vb1y, vb1z, rrvb1, rr2vb1, f, thr);
    invang<EVFLAG,EFLAG,NEWTON_BOND>(i1,i3,i4,i2, type,
           vb1x, vb1y, vb1z, rrvb1, rr2vb1,
           vb3x, vb3y, vb3z, rrvb3, rr2vb3,
           vb2x, vb2y, vb2z, rrvb2, rr2vb2, f, thr);
    invang<EVFLAG,EFLAG,NEWTON_BOND>(i1,i4,i2,i3, type,
           vb2x, vb2y, vb2z, rrvb2, rr2vb2,
           vb1x, vb1y, vb1z, rrvb1, rr2vb1,
           vb3x, vb3y, vb3z, rrvb3, rr2vb3, f, thr);
  }
}

/* ----------------------------------------------------------------------
   compute inversion angles + energy and forces
------------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_BOND>
void ImproperInversionHarmonicOMP::invang(const int &i1, const int &i2,
          const int &i3, const int &i4,
          const int &type,
          const double &vb1x, const double &vb1y, const double &vb1z,
          const double &rrvb1, const double &rr2vb1,
          const double &vb2x, const double &vb2y, const double &vb2z,
          const double &rrvb2, const double &rr2vb2,
          const double &vb3x, const double &vb3y, const double &vb3z,
          const double &rrvb3, const double &rr2vb3,
          dbl3_t * _noalias const f, ThrData * const thr)
{
  double eimproper,f1[3],f2[3],f3[3],f4[3];
  double omega,cosomega,domega,gomega,rjk,rjl;
  double upx,upy,upz,upn,rup,umx,umy,umz,umn,rum,wwr;
  double rucb,rudb,rvcb,rvdb,rupupn,rumumn;

  const int nlocal = atom->nlocal;

  eimproper = 0.0;

  // scalar products of IJ*IK and IJ*IL
  rjk=vb3x*vb2x+vb3y*vb2y+vb3z*vb2z;
  rjl=vb1x*vb3x+vb1y*vb3y+vb1z*vb3z;

  // unit-vector: IK+IL
  upx=vb2x*rrvb2+vb1x*rrvb1;
  upy=vb2y*rrvb2+vb1y*rrvb1;
  upz=vb2z*rrvb2+vb1z*rrvb1;
  upn=1.0/sqrt(upx*upx+upy*upy+upz*upz);
  upx=upx*upn;
  upy=upy*upn;
  upz=upz*upn;
  rup=vb3x*upx+vb3y*upy+vb3z*upz;

  // unit-vector: IK-IL
  umx=vb2x*rrvb2-vb1x*rrvb1;
  umy=vb2y*rrvb2-vb1y*rrvb1;
  umz=vb2z*rrvb2-vb1z*rrvb1;
  umn=1.0/sqrt(umx*umx+umy*umy+umz*umz);
  umx=umx*umn;
  umy=umy*umn;
  umz=umz*umn;
  rum=vb3x*umx+vb3y*umy+vb3z*umz;

  // angle theta
  wwr=sqrt(rup*rup+rum*rum);

  cosomega=wwr*rrvb3;

  if (cosomega > 1.0) cosomega = 1.0;

  omega=acos(cosomega);

  domega = acos(cosomega)-w0[type];
  if (EFLAG) eimproper = kw[type]*(domega*domega);

  // kw[type] is divided by 3 -> threefold contribution
  gomega=0.0;
  if (omega*omega > 1.0e-24)  gomega=2.0*kw[type]*(domega)/(sin(omega));

  // projection IK and IL on unit vectors and contribution on IK and IL
  rucb = rjk-rup*(vb2x*upx+vb2y*upy+vb2z*upz);
  rudb = rjl-rup*(vb1x*upx+vb1y*upy+vb1z*upz);
  rvcb = rjk-rum*(vb2x*umx+vb2y*umy+vb2z*umz);
  rvdb = rjl-rum*(vb1x*umx+vb1y*umy+vb1z*umz);

  rupupn = rup*upn;
  rumumn = rum*umn;

  // force contributions of angle
  f2[0]=gomega*(-cosomega*vb3x*rr2vb3+rrvb3*(rup*upx+rum*umx)/wwr);
  f2[1]=gomega*(-cosomega*vb3y*rr2vb3+rrvb3*(rup*upy+rum*umy)/wwr);
  f2[2]=gomega*(-cosomega*vb3z*rr2vb3+rrvb3*(rup*upz+rum*umz)/wwr);

  f3[0]=gomega*rrvb3*(rupupn*rrvb2*(vb3x-rup*upx-rucb*vb2x*rr2vb2) +
        rumumn*rrvb2*(vb3x-rum*umx-rvcb*vb2x*rr2vb2))/wwr;
  f3[1]=gomega*rrvb3*(rupupn*rrvb2*(vb3y-rup*upy-rucb*vb2y*rr2vb2) +
        rumumn*rrvb2*(vb3y-rum*umy-rvcb*vb2y*rr2vb2))/wwr;
  f3[2]=gomega*rrvb3*(rupupn*rrvb2*(vb3z-rup*upz-rucb*vb2z*rr2vb2) +
        rumumn*rrvb2*(vb3z-rum*umz-rvcb*vb2z*rr2vb2))/wwr;

  f4[0]=gomega*rrvb3*(rupupn*rrvb1*(vb3x-rup*upx-rudb*vb1x*rr2vb1) -
        rumumn*rrvb1*(vb3x-rum*umx-rvdb*vb1x*rr2vb1))/wwr;
  f4[1]=gomega*rrvb3*(rupupn*rrvb1*(vb3y-rup*upy-rudb*vb1y*rr2vb1) -
        rumumn*rrvb1*(vb3y-rum*umy-rvdb*vb1y*rr2vb1))/wwr;
  f4[2]=gomega*rrvb3*(rupupn*rrvb1*(vb3z-rup*upz-rudb*vb1z*rr2vb1) -
        rumumn*rrvb1*(vb3z-rum*umz-rvdb*vb1z*rr2vb1))/wwr;

  f1[0] = -(f2[0] + f3[0] + f4[0]);
  f1[1] = -(f2[1] + f3[1] + f4[1]);
  f1[2] = -(f2[2] + f3[2] + f4[2]);

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

  if (EVFLAG) {
    double rb3x, rb3y, rb3z;

    rb3x = vb1x - vb2x;
    rb3y = vb1y - vb2y;
    rb3z = vb1z - vb2z;

    ev_tally_thr(this,i1,i2,i3,i4,nlocal,NEWTON_BOND,eimproper,f2,f3,f4,
                 vb3x,vb3y,vb3z,vb2x,vb2y,vb2z,rb3x,rb3y,rb3z,thr);
  }
}
