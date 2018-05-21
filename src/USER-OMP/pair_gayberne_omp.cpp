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

#include <cmath>
#include "pair_gayberne_omp.h"
#include "math_extra.h"
#include "atom.h"
#include "comm.h"
#include "atom_vec_ellipsoid.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"

#include "suffix.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairGayBerneOMP::PairGayBerneOMP(LAMMPS *lmp) :
  PairGayBerne(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairGayBerneOMP::compute(int eflag, int vflag)
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
void PairGayBerneOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double evdwl,one_eng,rsq,r2inv,r6inv,forcelj,factor_lj;
  double fforce[3],ttor[3],rtor[3],r12[3];
  double a1[3][3],b1[3][3],g1[3][3],a2[3][3],b2[3][3],g2[3][3],temp[3][3];
  int *ilist,*jlist,*numneigh,**firstneigh;
  double *iquat,*jquat;

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];
  dbl3_t * _noalias const tor = (dbl3_t *) thr->get_torque()[0];
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * _noalias const special_lj = force->special_lj;
  const int * _noalias const ellipsoid = atom->ellipsoid;

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;

  double fxtmp,fytmp,fztmp,t1tmp,t2tmp,t3tmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  evdwl = 0.0;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    itype = type[i];
    fxtmp=fytmp=fztmp=t1tmp=t2tmp=t3tmp=0.0;

    if (form[itype][itype] == ELLIPSE_ELLIPSE) {
      iquat = bonus[ellipsoid[i]].quat;
      MathExtra::quat_to_mat_trans(iquat,a1);
      MathExtra::diag_times3(well[itype],a1,temp);
      MathExtra::transpose_times3(a1,temp,b1);
      MathExtra::diag_times3(shape2[itype],a1,temp);
      MathExtra::transpose_times3(a1,temp,g1);
    }

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

      if (rsq < cutsq[itype][jtype]) {

        switch (form[itype][jtype]) {
        case SPHERE_SPHERE:
          r2inv = 1.0/rsq;
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
          forcelj *= -r2inv;
          if (EFLAG)
            one_eng = r6inv*(r6inv*lj3[itype][jtype]-lj4[itype][jtype]) -
              offset[itype][jtype];
          fforce[0] = r12[0]*forcelj;
          fforce[1] = r12[1]*forcelj;
          fforce[2] = r12[2]*forcelj;
          ttor[0] = ttor[1] = ttor[2] = 0.0;
          rtor[0] = rtor[1] = rtor[2] = 0.0;
          break;

        case SPHERE_ELLIPSE:
          jquat = bonus[ellipsoid[j]].quat;
          MathExtra::quat_to_mat_trans(jquat,a2);
          MathExtra::diag_times3(well[jtype],a2,temp);
          MathExtra::transpose_times3(a2,temp,b2);
          MathExtra::diag_times3(shape2[jtype],a2,temp);
          MathExtra::transpose_times3(a2,temp,g2);
          one_eng = gayberne_lj(j,i,a2,b2,g2,r12,rsq,fforce,rtor);
          ttor[0] = ttor[1] = ttor[2] = 0.0;
          break;

        case ELLIPSE_SPHERE:
          one_eng = gayberne_lj(i,j,a1,b1,g1,r12,rsq,fforce,ttor);
          rtor[0] = rtor[1] = rtor[2] = 0.0;
          break;

        default:
          jquat = bonus[ellipsoid[j]].quat;
          MathExtra::quat_to_mat_trans(jquat,a2);
          MathExtra::diag_times3(well[jtype],a2,temp);
          MathExtra::transpose_times3(a2,temp,b2);
          MathExtra::diag_times3(shape2[jtype],a2,temp);
          MathExtra::transpose_times3(a2,temp,g2);
          one_eng = gayberne_analytic(i,j,a1,a2,b1,b2,g1,g2,r12,rsq,
                                      fforce,ttor,rtor);
          break;
        }

        fforce[0] *= factor_lj;
        fforce[1] *= factor_lj;
        fforce[2] *= factor_lj;
        ttor[0] *= factor_lj;
        ttor[1] *= factor_lj;
        ttor[2] *= factor_lj;

        fxtmp += fforce[0];
        fytmp += fforce[1];
        fztmp += fforce[2];
        t1tmp += ttor[0];
        t2tmp += ttor[1];
        t3tmp += ttor[2];

        if (NEWTON_PAIR || j < nlocal) {
          rtor[0] *= factor_lj;
          rtor[1] *= factor_lj;
          rtor[2] *= factor_lj;
          f[j].x -= fforce[0];
          f[j].y -= fforce[1];
          f[j].z -= fforce[2];
          tor[j].x += rtor[0];
          tor[j].y += rtor[1];
          tor[j].z += rtor[2];
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

double PairGayBerneOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairGayBerne::memory_usage();

  return bytes;
}
