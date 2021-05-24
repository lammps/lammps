// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_resquared_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_extra.h"
#include "neigh_list.h"
#include "suffix.h"

#include "omp_compat.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairRESquaredOMP::PairRESquaredOMP(LAMMPS *lmp) :
  PairRESquared(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairRESquaredOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    if (evflag) {
      if (eflag) {
        if (force->newton_pair) eval<1,1,1>(ifrom, ito, thr);
        else eval<1,1,0>(ifrom, ito, thr);
      } else {
        if (force->newton_pair) eval<1,0,1>(ifrom, ito, thr);
        else eval<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair) eval<0,0,1>(ifrom, ito, thr);
      else eval<0,0,0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairRESquaredOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double evdwl,one_eng,rsq,r2inv,r6inv,forcelj,factor_lj;
  double fforce[3],ttor[3],rtor[3],r12[3];
  int *ilist,*jlist,*numneigh,**firstneigh;
  RE2Vars wi,wj;

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];
  dbl3_t * _noalias const tor = (dbl3_t *) thr->get_torque()[0];
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * _noalias const special_lj = force->special_lj;

  double fxtmp,fytmp,fztmp,t1tmp,t2tmp,t3tmp;

  evdwl = 0.0;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    itype = type[i];
    fxtmp=fytmp=fztmp=t1tmp=t2tmp=t3tmp=0.0;

    // not a LJ sphere

    if (lshape[itype] != 0.0) precompute_i(i,wi);

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      // r12 = center to center vector

      r12[0] = x[j].x-x[i].x;
      r12[1] = x[j].y-x[i].y;
      r12[2] = x[j].z-x[i].z;
      rsq = MathExtra::dot3(r12,r12);
      jtype = type[j];

      // compute if less than cutoff

      if (rsq < cutsq[itype][jtype]) {
        fforce[0] = fforce[1] = fforce[2] = 0.0;

        switch (form[itype][jtype]) {

         case SPHERE_SPHERE:
          r2inv = 1.0/rsq;
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
          forcelj *= -r2inv;
          if (EFLAG) one_eng =
              r6inv*(r6inv*lj3[itype][jtype]-lj4[itype][jtype]) -
              offset[itype][jtype];
          fforce[0] = r12[0]*forcelj;
          fforce[1] = r12[1]*forcelj;
          fforce[2] = r12[2]*forcelj;
          break;

         case SPHERE_ELLIPSE:
          precompute_i(j,wj);
          if (NEWTON_PAIR || j < nlocal) {
            one_eng = resquared_lj(j,i,wj,r12,rsq,fforce,rtor,true);
            tor[j].x += rtor[0]*factor_lj;
            tor[j].y += rtor[1]*factor_lj;
            tor[j].z += rtor[2]*factor_lj;
          } else
            one_eng = resquared_lj(j,i,wj,r12,rsq,fforce,rtor,false);
          break;

         case ELLIPSE_SPHERE:
          one_eng = resquared_lj(i,j,wi,r12,rsq,fforce,ttor,true);
          t1tmp += ttor[0]*factor_lj;
          t2tmp += ttor[1]*factor_lj;
          t3tmp += ttor[2]*factor_lj;
          break;

         default:
          precompute_i(j,wj);
          one_eng = resquared_analytic(i,j,wi,wj,r12,rsq,fforce,ttor,rtor);
          t1tmp += ttor[0]*factor_lj;
          t2tmp += ttor[1]*factor_lj;
          t3tmp += ttor[2]*factor_lj;
          if (NEWTON_PAIR || j < nlocal) {
            tor[j].x += rtor[0]*factor_lj;
            tor[j].y += rtor[1]*factor_lj;
            tor[j].z += rtor[2]*factor_lj;
          }
         break;
        }

        fforce[0] *= factor_lj;
        fforce[1] *= factor_lj;
        fforce[2] *= factor_lj;
        fxtmp += fforce[0];
        fytmp += fforce[1];
        fztmp += fforce[2];

        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= fforce[0];
          f[j].y -= fforce[1];
          f[j].z -= fforce[2];
        }

        if (EFLAG) evdwl = factor_lj*one_eng;

        if (EVFLAG) ev_tally_xyz_thr(this,i,j,nlocal,NEWTON_PAIR,
                                     evdwl,0.0,fforce[0],fforce[1],fforce[2],
                                     -r12[0],-r12[1],-r12[2],thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
    tor[i].x += t1tmp;
    tor[i].y += t2tmp;
    tor[i].z += t3tmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairRESquaredOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairRESquared::memory_usage();

  return bytes;
}
