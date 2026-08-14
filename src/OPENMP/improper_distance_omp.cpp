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

#include "improper_distance_omp.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "neighbor.h"

#include <cmath>

#include "omp_compat.h"
#include "suffix.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ImproperDistanceOMP::ImproperDistanceOMP(class LAMMPS *lmp)
  : ImproperDistance(lmp), ThrOMP(lmp,THR_IMPROPER)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void ImproperDistanceOMP::compute(int eflag, int vflag)
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
void ImproperDistanceOMP::eval(int nfrom, int nto, ThrData * const thr)
{
  int i1,i2,i3,i4,n,type;
  double xab, yab, zab; // bond 1-2
  double xac, yac, zac; // bond 1-3
  double xad, yad, zad; // bond 1-4
  double xbc, ybc, zbc; // bond 2-3
  double xbd, ybd, zbd; // bond 2-4
  double xna, yna, zna, rna; // normal
  double da;
  double eimproper,f1[3],f2[3],f3[3],f4[3];
  double domega,a;

  eimproper = 0.0;

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const int5_t * _noalias const improperlist = (int5_t *) neighbor->improperlist[0];
  const int nlocal = atom->nlocal;

  for (n = nfrom; n < nto; n++) {
    i1 = improperlist[n].a;
    i2 = improperlist[n].b;
    i3 = improperlist[n].c;
    i4 = improperlist[n].d;
    type = improperlist[n].t;

    // geometry of 4-body
    // 1 is the central atom

    // bond 1->2
    xab = x[i2].x - x[i1].x;
    yab = x[i2].y - x[i1].y;
    zab = x[i2].z - x[i1].z;
    domain->minimum_image(FLERR, xab,yab,zab);

    // bond 1->3
    xac = x[i3].x - x[i1].x;
    yac = x[i3].y - x[i1].y;
    zac = x[i3].z - x[i1].z;
    domain->minimum_image(FLERR, xac,yac,zac);

    // bond 1->4
    xad = x[i4].x - x[i1].x;
    yad = x[i4].y - x[i1].y;
    zad = x[i4].z - x[i1].z;
    domain->minimum_image(FLERR, xad,yad,zad);

    // bond 2-3
    xbc = x[i3].x - x[i2].x;
    ybc = x[i3].y - x[i2].y;
    zbc = x[i3].z - x[i2].z;
    domain->minimum_image(FLERR, xbc,ybc,zbc);

    // bond 2-4
    xbd = x[i4].x - x[i2].x;
    ybd = x[i4].y - x[i2].y;
    zbd = x[i4].z - x[i2].z;
    domain->minimum_image(FLERR, xbd,ybd,zbd);

    xna =   ybc*zbd - zbc*ybd;
    yna = -(xbc*zbd - zbc*xbd);
    zna =   xbc*ybd - ybc*xbd;
    rna = 1.0 / sqrt(xna*xna+yna*yna+zna*zna);
    xna *= rna;
    yna *= rna;
    zna *= rna;

    da = xna*xab + yna*yab + zna*zab;

    domega = k[type]*da*da + chi[type]*da*da*da*da;
    a =  2.0* (k[type]*da + 2.0*chi[type]*da*da*da);

    if (EFLAG) eimproper = domega;

    f1[0] = a*( xna);
    f1[1] = a*( yna);
    f1[2] = a*( zna);

    f2[0] = a*( -xna               -yab*(zbd-zbc)*rna +zab*(ybd-ybc)*rna -da*( -yna*(zbd-zbc) + zna*(ybd-ybc) )*rna);
    f2[1] = a*( +xab*(zbd-zbc)*rna -yna               +zab*(xbc-xbd)*rna -da*( +xna*(zbd-zbc) + zna*(xbc-xbd) )*rna);
    f2[2] = a*( -xab*(ybd-ybc)*rna -yab*(xbc-xbd)*rna -zna               -da*( +xna*(ybc-ybd) - yna*(xbc-xbd) )*rna);

    f3[0] = a*( (           yab*zbd -zab*ybd ) *rna +da*( -yna*zbd +zna*ybd )*rna);
    f3[1] = a*( ( -xab*zbd          +zab*xbd ) *rna +da*( +xna*zbd -zna*xbd )*rna);
    f3[2] = a*( ( +xab*ybd -yab*xbd          ) *rna +da*( -xna*ybd +yna*xbd )*rna);

    f4[0] = a*( (          -yab*zbc +zab*ybc ) *rna -da*( -yna*zbc +zna*ybc )*rna);
    f4[1] = a*( ( +xab*zbc          -zab*xbc ) *rna -da*( +xna*zbc -zna*xbc )*rna);
    f4[2] = a*( ( -xab*ybc +yab*xbc          ) *rna -da*( -xna*ybc +yna*xbc )*rna);

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
      ev_tally_thr(this,i1,i2,i3,i4,nlocal,NEWTON_BOND,eimproper,f2,f3,f4,
                   xab,yab,zab,xac,yac,zac,xad-xac,yad-yac,zad-zac,thr);
  }
}
