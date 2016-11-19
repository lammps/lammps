/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include <math.h>
#include <string.h>
#include <stdint.h>
#include "pair_agni_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "math_special.h"
#include "math_const.h"

#include "suffix.h"
using namespace LAMMPS_NS;
using namespace MathSpecial;

/*
   Copyright (c) 2012,2013   Axel Kohlmeyer <akohlmey@gmail.com>
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
   * Neither the name of the <organization> nor the
     names of its contributors may be used to endorse or promote products
     derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/* faster versions of 2**x, e**x, and 10**x in single and double precision.
 *
 * Based on the Cephes math library 2.8
 */

/* internal definitions for the fastermath library */

/* IEEE 754 double precision floating point data manipulation */
typedef union
{
    double   f;
    uint64_t u;
    struct {int32_t  i0,i1;};
}  udi_t;
#define FM_DOUBLE_BIAS 1023
#define FM_DOUBLE_EMASK 2146435072
#define FM_DOUBLE_MBITS 20
#define FM_DOUBLE_MMASK 1048575
#define FM_DOUBLE_EZERO 1072693248

/* generate 2**num in floating point by bitshifting */
#define FM_DOUBLE_INIT_EXP(var,num)                 \
    var.i0 = 0;                                     \
    var.i1 = (((int) num) + FM_DOUBLE_BIAS) << 20

/* double precision constants */
#define FM_DOUBLE_LOG2OFE  1.4426950408889634074
#define FM_DOUBLE_LOGEOF2  6.9314718055994530942e-1
#define FM_DOUBLE_LOG2OF10 3.32192809488736234789
#define FM_DOUBLE_LOG10OF2 3.0102999566398119521e-1
#define FM_DOUBLE_LOG10OFE 4.3429448190325182765e-1
#define FM_DOUBLE_SQRT2    1.41421356237309504880
#define FM_DOUBLE_SQRTH    0.70710678118654752440

/* optimizer friendly implementation of exp2(x).
 *
 * strategy:
 *
 * split argument into an integer part and a fraction:
 * ipart = floor(x+0.5);
 * fpart = x - ipart;
 *
 * compute exp2(ipart) from setting the ieee754 exponent
 * compute exp2(fpart) using a pade' approximation for x in [-0.5;0.5[
 *
 * the result becomes: exp2(x) = exp2(ipart) * exp2(fpart)
 */

static const double fm_exp2_q[] = {
/*  1.00000000000000000000e0, */
    2.33184211722314911771e2,
    4.36821166879210612817e3
};
static const double fm_exp2_p[] = {
    2.30933477057345225087e-2,
    2.02020656693165307700e1,
    1.51390680115615096133e3
};

static double fm_exp2(double x)
{
    double   ipart, fpart, px, qx;
    udi_t    epart;

    ipart = floor(x+0.5);
    fpart = x - ipart;
    FM_DOUBLE_INIT_EXP(epart,ipart);

    x = fpart*fpart;

    px =        fm_exp2_p[0];
    px = px*x + fm_exp2_p[1];
    qx =    x + fm_exp2_q[0];
    px = px*x + fm_exp2_p[2];
    qx = qx*x + fm_exp2_q[1];

    px = px * fpart;

    x = 1.0 + 2.0*(px/(qx-px));
    return epart.f*x;
}

static double fm_exp(double x)
{
#if defined(__BYTE_ORDER__)
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    return fm_exp2(FM_DOUBLE_LOG2OFE * (x));
#endif
#endif
    return exp(x);
}

/* ---------------------------------------------------------------------- */

PairAGNIOMP::PairAGNIOMP(LAMMPS *lmp) :
  PairAGNI(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairAGNIOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else evflag = vflag_fdotr = 0;

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, thr);

    if (evflag) eval<1>(ifrom, ito, thr);
    else eval<0>(ifrom, ito, thr);

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int EVFLAG>
void PairAGNIOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,k,ii,jj,itype,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const int * _noalias const type = atom->type;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double fxtmp,fytmp,fztmp;
  double *Vx, *Vy, *Vz;

  // loop over full neighbor list of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    itype = map[type[i]];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    fxtmp = fytmp = fztmp = 0.0;

    const Param &iparam = params[elem2param[itype]];
    Vx = new double[iparam.numeta];
    Vy = new double[iparam.numeta];
    Vz = new double[iparam.numeta];
    memset(Vx,0,iparam.numeta*sizeof(double));
    memset(Vy,0,iparam.numeta*sizeof(double));
    memset(Vz,0,iparam.numeta*sizeof(double));

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;

      if ((rsq > 0.0) && (rsq < iparam.cutsq)) {
        const double r = sqrt(rsq);
        const double cF = 0.5*(cos((MathConst::MY_PI*r)/iparam.cut)+1.0);
        const double wX = cF*delx/r;
        const double wY = cF*dely/r;
        const double wZ = cF*delz/r;

        for (k = 0; k < iparam.numeta; ++k) {
          const double e = fm_exp(-(iparam.eta[k]*rsq));
          Vx[k] += wX*e;
          Vy[k] += wY*e;
          Vz[k] += wZ*e;
        }
      }
    }

    for (j = 0; j < iparam.numtrain; ++j) {
      double kx = 0.0;
      double ky = 0.0;
      double kz = 0.0;

      for(int k = 0; k < iparam.numeta; ++k) {
        const double xu = iparam.xU[k][j];
        kx += square(Vx[k] - xu);
        ky += square(Vy[k] - xu);
        kz += square(Vz[k] - xu);
      }
      const double e = -0.5/(square(iparam.sigma));
      fxtmp += iparam.alpha[j]*fm_exp(kx*e);
      fytmp += iparam.alpha[j]*fm_exp(ky*e);
      fztmp += iparam.alpha[j]*fm_exp(kz*e);
    }
    fxtmp += iparam.b;
    fytmp += iparam.b;
    fztmp += iparam.b;
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;

    if (EVFLAG) ev_tally_xyz_full_thr(this,i,0.0,0.0,
                                      fxtmp,fytmp,fztmp,
                                      delx,dely,delz,thr);
    delete [] Vx;
    delete [] Vy;
    delete [] Vz;
  }
}

/* ---------------------------------------------------------------------- */

double PairAGNIOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairAGNI::memory_usage();

  return bytes;
}
