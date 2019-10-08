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
#include "improper_cossq_omp.h"
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

ImproperCossqOMP::ImproperCossqOMP(class LAMMPS *lmp)
  : ImproperCossq(lmp), ThrOMP(lmp,THR_IMPROPER)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void ImproperCossqOMP::compute(int eflag, int vflag)
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
void ImproperCossqOMP::eval(int nfrom, int nto, ThrData * const thr)
{
  int i1,i2,i3,i4,n,type;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z;
  double eimproper,f1[3],f2[3],f3[3],f4[3];
  double rjisq, rji, rlksq, rlk, cosphi, angfac;
  double cjiji, clkji, clklk, cfact1, cfact2, cfact3;

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

    /* separation vector between i1 and i2, (i2-i1) */
    vb1x = x[i2].x - x[i1].x;
    vb1y = x[i2].y - x[i1].y;
    vb1z = x[i2].z - x[i1].z;
    rjisq = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z ;
    rji = sqrt(rjisq);

    /* separation vector between i2 and i3 (i3-i2) */
    vb2x = x[i3].x - x[i2].x;
    vb2y = x[i3].y - x[i2].y;
    vb2z = x[i3].z - x[i2].z;

    /* separation vector between i3 and i4, (i4-i3) */
    vb3x = x[i4].x - x[i3].x;
    vb3y = x[i4].y - x[i3].y;
    vb3z = x[i4].z - x[i3].z;
    rlksq = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z ;
    rlk = sqrt(rlksq);

    cosphi = (vb3x*vb1x + vb3y*vb1y + vb3z*vb1z)/(rji * rlk);

    /* Check that cos(phi) is in the correct limits. */
    if (cosphi > 1.0 + TOLERANCE || cosphi < (-1.0 - TOLERANCE)) {
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


      /* Apply corrections to round-off errors. */
      if (cosphi > 1.0)  cosphi -= SMALL;
      if (cosphi < -1.0) cosphi += SMALL;

      /* Calculate the angle: */
      double torangle = acos(cosphi);
      cosphi = cos(torangle - chi[type]);

      if (EFLAG) eimproper = 0.5 * k[type] * cosphi * cosphi;

      /*
        printf("The tags: %d-%d-%d-%d, of type %d .\n",atom->tag[i1],atom->tag[i2],atom->tag[i3],atom->tag[i4],type);
        printf("The ji vector: %f, %f, %f.\nThe lk vector: %f, %f, %f.\n", vb1x,vb1y,vb1z,vb3x,vb3y,vb3z);
        printf("The cosine of the angle: %-1.16e.\n", cosphi);
        printf("The energy of the improper: %-1.16e with prefactor %-1.16e.\n", eimproper, 0.5*k[type]);
      */

      /* Work out forces. */
      angfac = - k[type] * cosphi;

      cjiji = rjisq;
      clklk = rlksq;
      /*CLKJI = RXLK * RXJI + RYLK * RYJI + RZLK * RZJI */
      clkji = vb3x*vb1x + vb3y*vb1y + vb3z*vb1z;

      /*CFACT1 = CLKLK * CJIJI
        CFACT1 = SQRT(CFACT1)
        CFACT1 = ANGFAC / CFACT1*/
      cfact1 = angfac / sqrt(clklk * cjiji);
      /*CFACT2 = CLKJI / CLKLK*/
      cfact2 = clkji / clklk;
      /*CFACT3 = CLKJI / CJIJI*/
      cfact3 = clkji / cjiji;

      /*FIX = -RXLK + CFACT3 * RXJI
        FIY = -RYLK + CFACT3 * RYJI
        FIZ = -RZLK + CFACT3 * RZJI*/
      f1[0] = - vb3x + cfact3 * vb1x;
      f1[1] = - vb3y + cfact3 * vb1y;
      f1[2] = - vb3z + cfact3 * vb1z;

      /*FJX = -FIX
        FJY = -FIY
        FJZ = -FIZ*/
      f2[0] = - f1[0];
      f2[1] = - f1[1];
      f2[2] = - f1[2];

      /*FKX = CFACT2 * RXLK - RXJI
        FKY = CFACT2 * RYLK - RYJI
        FKZ = CFACT2 * RZLK - RZJI*/
      f3[0] = cfact2 * vb3x - vb1x;
      f3[1] = cfact2 * vb3y - vb1y;
      f3[2] = cfact2 * vb3z - vb1z;

      /*FLX = -FKX
        FLY = -FKY
        FLZ = -FKZ*/
      f4[0] = - f3[0];
      f4[1] = - f3[1];
      f4[2] = - f3[2];

      /*FIX = FIX * CFACT1
        FIY = FIY * CFACT1
        FIZ = FIZ * CFACT1*/
      f1[0] *= cfact1;
      f1[1] *= cfact1;
      f1[2] *= cfact1;

      /*FJX = FJX * CFACT1
        FJY = FJY * CFACT1
        FJZ = FJZ * CFACT1*/
      f2[0] *= cfact1;
      f2[1] *= cfact1;
      f2[2] *= cfact1;

      /*FKX = FKX * CFACT1
        FKY = FKY * CFACT1
        FKZ = FKZ * CFACT1*/
      f3[0] *= cfact1;
      f3[1] *= cfact1;
      f3[2] *= cfact1;

      /*FLX = FLX * CFACT1
        FLY = FLY * CFACT1
        FLZ = FLZ * CFACT1*/
      f4[0] *= cfact1;
      f4[1] *= cfact1;
      f4[2] *= cfact1;

      /* Apply force to each of 4 atoms */
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
        ev_tally_thr(this,i1,i2,i3,i4,nlocal,NEWTON_BOND,eimproper,f1,f3,f4,
                     -vb1x,-vb1y,-vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,thr);
    }
  }
}
