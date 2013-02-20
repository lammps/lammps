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

#include "math.h"
#include "pair_gauss_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"

#include "suffix.h"
using namespace LAMMPS_NS;

#define EPSILON 1.0e-10
/* ---------------------------------------------------------------------- */

PairGaussOMP::PairGaussOMP(LAMMPS *lmp) :
  PairGauss(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairGaussOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else evflag = vflag_fdotr = 0;

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;
  double occ = 0.0;

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(eflag,vflag) reduction(+:occ)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, thr);

    if (evflag) {
      if (eflag) {
        if (force->newton_pair) occ = eval<1,1,1>(ifrom, ito, thr);
        else occ = eval<1,1,0>(ifrom, ito, thr);
      } else {
        if (force->newton_pair) occ = eval<1,0,1>(ifrom, ito, thr);
        else occ = eval<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair) occ = eval<0,0,1>(ifrom, ito, thr);
      else occ = eval<0,0,0>(ifrom, ito, thr);
    }

    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region

  if (eflag_global) pvector[0] = occ;
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
double PairGaussOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,forcelj,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int occ = 0;

  evdwl = 0.0;

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * _noalias const special_lj = force->special_lj;
  double fxtmp,fytmp,fztmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp=fytmp=fztmp=0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      // define a Gaussian well to be occupied if
      // the site it interacts with is within the force maximum

      if (EFLAG)
        if (eflag_global && rsq < 0.5/b[itype][jtype]) occ++;

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        forcelj = - 2.0*a[itype][jtype]*b[itype][jtype] * rsq *
          exp(-b[itype][jtype]*rsq);
        fpair = factor_lj*forcelj*r2inv;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        if (EFLAG) {
          evdwl = -(a[itype][jtype]*exp(-b[itype][jtype]*rsq) -
                    offset[itype][jtype]);
          evdwl *= factor_lj;
        }

        if (EVFLAG) ev_tally_thr(this, i,j,nlocal,NEWTON_PAIR,
                                 evdwl,0.0,fpair,delx,dely,delz,thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
  return occ;
}

/* ---------------------------------------------------------------------- */

double PairGaussOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairGauss::memory_usage();

  return bytes;
}
