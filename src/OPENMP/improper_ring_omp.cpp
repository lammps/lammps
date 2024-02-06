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

#include "improper_ring_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_special.h"
#include "neighbor.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace MathSpecial;

static constexpr double SMALL =     0.001;

/* ---------------------------------------------------------------------- */

ImproperRingOMP::ImproperRingOMP(class LAMMPS *lmp)
  : ImproperRing(lmp), ThrOMP(lmp,THR_IMPROPER)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void ImproperRingOMP::compute(int eflag, int vflag)
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
    thr->timer(Timer::BOND);
      reduce_thr(this, eflag, vflag, thr);
    }
  } // end of omp parallel region
}

template <int EVFLAG, int EFLAG, int NEWTON_BOND>
void ImproperRingOMP::eval(int nfrom, int nto, ThrData * const thr)
{
  /* Be careful!: "chi" is the equilibrium angle in radians. */
  int i1,i2,i3,i4,n,type;

  double eimproper;

  /* Compatibility variables. */
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z;
  double f1[3], f3[3], f4[3];

  /* Actual computation variables. */
  int at1[3], at2[3], at3[3], icomb;
  double   bvec1x[3], bvec1y[3], bvec1z[3],
    bvec2x[3], bvec2y[3], bvec2z[3],
    bvec1n[3], bvec2n[3], bend_angle[3];
  double   angle_summer, angfac, cfact1, cfact2, cfact3;
  double   cjiji, ckjji, ckjkj, fix, fiy, fiz, fjx, fjy, fjz, fkx, fky, fkz;

  eimproper = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const int * const * const improperlist = neighbor->improperlist;
  const int nlocal = atom->nlocal;

  /* A description of the potential can be found in
     Macromolecules 35, pp. 1463-1472 (2002). */
  for (n = nfrom; n < nto; n++) {
    /* Take the ids of the atoms contributing to the improper potential. */
    i1 = improperlist[n][0];   /* Atom "1" of Figure 1 from the above reference.*/
    i2 = improperlist[n][1];   /* Atom "2" ... */
    i3 = improperlist[n][2];   /* Atom "3" ... */
    i4 = improperlist[n][3];   /* Atom "9" ... */
    type = improperlist[n][4];

    /* Calculate the necessary variables for LAMMPS implementation.
       if (evflag) ev_tally(i1,i2,i3,i4,nlocal,newton_bond,eimproper,f1,f3,f4,
       vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z);
       Although, they are irrelevant to the calculation of the potential, we keep
       them for maximal compatibility. */
    vb1x = x[i1][0] - x[i2][0]; vb1y = x[i1][1] - x[i2][1]; vb1z = x[i1][2] - x[i2][2];

    vb2x = x[i3][0] - x[i2][0]; vb2y = x[i3][1] - x[i2][1]; vb2z = x[i3][2] - x[i2][2];

    vb3x = x[i4][0] - x[i3][0]; vb3y = x[i4][1] - x[i3][1]; vb3z = x[i4][2] - x[i3][2];


    /* Pass the atom tags to form the necessary combinations. */
    at1[0] = i1; at2[0] = i2; at3[0] = i4;  /* ids: 1-2-9 */
    at1[1] = i1; at2[1] = i2; at3[1] = i3;  /* ids: 1-2-3 */
    at1[2] = i4; at2[2] = i2; at3[2] = i3;  /* ids: 9-2-3 */


    /* Initialize the sum of the angles differences. */
    angle_summer = 0.0;
    /* Take a loop over the three angles, defined by each triad: */
    for (icomb = 0; icomb < 3; icomb ++) {

      /* Bond vector connecting the first and the second atom. */
      bvec1x[icomb] = x[at2[icomb]][0] - x[at1[icomb]][0];
      bvec1y[icomb] = x[at2[icomb]][1] - x[at1[icomb]][1];
      bvec1z[icomb] = x[at2[icomb]][2] - x[at1[icomb]][2];
      /* also calculate the norm of the vector: */
      bvec1n[icomb] = sqrt(  bvec1x[icomb]*bvec1x[icomb]
                             + bvec1y[icomb]*bvec1y[icomb]
                             + bvec1z[icomb]*bvec1z[icomb]);
      /* Bond vector connecting the second and the third atom. */
      bvec2x[icomb] = x[at3[icomb]][0] - x[at2[icomb]][0];
      bvec2y[icomb] = x[at3[icomb]][1] - x[at2[icomb]][1];
      bvec2z[icomb] = x[at3[icomb]][2] - x[at2[icomb]][2];
      /* also calculate the norm of the vector: */
      bvec2n[icomb] = sqrt(  bvec2x[icomb]*bvec2x[icomb]
                             + bvec2y[icomb]*bvec2y[icomb]
                             + bvec2z[icomb]*bvec2z[icomb]);

      /* Calculate the bending angle of the atom triad: */
      bend_angle[icomb] = (  bvec2x[icomb]*bvec1x[icomb]
                             + bvec2y[icomb]*bvec1y[icomb]
                             + bvec2z[icomb]*bvec1z[icomb]);
      bend_angle[icomb] /= (bvec1n[icomb] * bvec2n[icomb]);
      if (bend_angle[icomb] >  1.0) bend_angle[icomb] -= SMALL;
      if (bend_angle[icomb] < -1.0) bend_angle[icomb] += SMALL;

      /* Append the current angle to the sum of angle differences. */
      angle_summer += (bend_angle[icomb] - chi[type]);
    }
    if (EFLAG) eimproper = (1.0/6.0) *k[type] * powint(angle_summer,6);
    /*
      printf("The tags: %d-%d-%d-%d, of type %d .\n",atom->tag[i1],atom->tag[i2],atom->tag[i3],atom->tag[i4],type);
      // printf("The coordinates of the first: %f, %f, %f.\n", x[i1][0], x[i1][1], x[i1][2]);
      // printf("The coordinates of the second: %f, %f, %f.\n", x[i2][0], x[i2][1], x[i2][2]);
      // printf("The coordinates of the third: %f, %f, %f.\n", x[i3][0], x[i3][1], x[i3][2]);
      // printf("The coordinates of the fourth: %f, %f, %f.\n", x[i4][0], x[i4][1], x[i4][2]);
      printf("The angles are: %f / %f / %f equilibrium: %f.\n", bend_angle[0], bend_angle[1], bend_angle[2],chi[type]);
      printf("The energy of the improper: %f with prefactor %f.\n", eimproper,(1.0/6.0)*k[type]);
      printf("The sum of the angles: %f.\n", angle_summer);
    */

    /* Force calculation acting on all atoms.
       Calculate the derivatives of the potential. */
    angfac = k[type] * powint(angle_summer,5);

    f1[0] = 0.0; f1[1] = 0.0; f1[2] = 0.0;
    f3[0] = 0.0; f3[1] = 0.0; f3[2] = 0.0;
    f4[0] = 0.0; f4[1] = 0.0; f4[2] = 0.0;

    /* Take a loop over the three angles, defined by each triad: */
    for (icomb = 0; icomb < 3; icomb ++)
      {
        /* Calculate the squares of the distances. */
        cjiji = bvec1n[icomb] * bvec1n[icomb];  ckjkj = bvec2n[icomb] * bvec2n[icomb];

        ckjji =   bvec2x[icomb] * bvec1x[icomb]
          + bvec2y[icomb] * bvec1y[icomb]
          + bvec2z[icomb] * bvec1z[icomb] ;

        cfact1 = angfac / (sqrt(ckjkj * cjiji));
        cfact2 = ckjji / ckjkj;
        cfact3 = ckjji / cjiji;

        /* Calculate the force acted on the third atom of the angle. */
        fkx = cfact2 * bvec2x[icomb] - bvec1x[icomb];
        fky = cfact2 * bvec2y[icomb] - bvec1y[icomb];
        fkz = cfact2 * bvec2z[icomb] - bvec1z[icomb];

        /* Calculate the force acted on the first atom of the angle. */
        fix = bvec2x[icomb] - cfact3 * bvec1x[icomb];
        fiy = bvec2y[icomb] - cfact3 * bvec1y[icomb];
        fiz = bvec2z[icomb] - cfact3 * bvec1z[icomb];

        /* Finally, calculate the force acted on the middle atom of the angle.*/
        fjx = - fix - fkx;  fjy = - fiy - fky;  fjz = - fiz - fkz;

        /* Consider the appropriate scaling of the forces: */
        fix *= cfact1; fiy *= cfact1; fiz *= cfact1;
        fjx *= cfact1; fjy *= cfact1; fjz *= cfact1;
        fkx *= cfact1; fky *= cfact1; fkz *= cfact1;

        if      (at1[icomb] == i1)  {f1[0] += fix; f1[1] += fiy; f1[2] += fiz;}
        else if (at2[icomb] == i1)  {f1[0] += fjx; f1[1] += fjy; f1[2] += fjz;}
        else if (at3[icomb] == i1)  {f1[0] += fkx; f1[1] += fky; f1[2] += fkz;}

        if      (at1[icomb] == i3)  {f3[0] += fix; f3[1] += fiy; f3[2] += fiz;}
        else if (at2[icomb] == i3)  {f3[0] += fjx; f3[1] += fjy; f3[2] += fjz;}
        else if (at3[icomb] == i3)  {f3[0] += fkx; f3[1] += fky; f3[2] += fkz;}

        if      (at1[icomb] == i4)  {f4[0] += fix; f4[1] += fiy; f4[2] += fiz;}
        else if (at2[icomb] == i4)  {f4[0] += fjx; f4[1] += fjy; f4[2] += fjz;}
        else if (at3[icomb] == i4)  {f4[0] += fkx; f4[1] += fky; f4[2] += fkz;}


        /* Store the contribution to the global arrays: */
        /* Take the id of the atom from the at1[icomb] element, i1 = at1[icomb]. */
        if (NEWTON_BOND || at1[icomb] < nlocal) {
          f[at1[icomb]][0] += fix;
          f[at1[icomb]][1] += fiy;
          f[at1[icomb]][2] += fiz;
        }
        /* Take the id of the atom from the at2[icomb] element, i2 = at2[icomb]. */
        if (NEWTON_BOND || at2[icomb] < nlocal) {
          f[at2[icomb]][0] += fjx;
          f[at2[icomb]][1] += fjy;
          f[at2[icomb]][2] += fjz;
        }
        /* Take the id of the atom from the at3[icomb] element, i3 = at3[icomb]. */
        if (NEWTON_BOND || at3[icomb] < nlocal) {
          f[at3[icomb]][0] += fkx;
          f[at3[icomb]][1] += fky;
          f[at3[icomb]][2] += fkz;
        }

      }

    if (EVFLAG)
      ev_tally_thr(this,i1,i2,i3,i4,nlocal,NEWTON_BOND,eimproper,f1,f3,f4,
                   vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,thr);
  }
}
