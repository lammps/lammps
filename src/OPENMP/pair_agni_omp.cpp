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

#include "pair_agni_omp.h"

#include "atom.h"
#include "comm.h"
#include "math_const.h"
#include "math_special.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>
#include <cstring>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace MathSpecial;

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

    if (atomic_feature_version == AGNI_VERSION_1) {
      if (evflag) eval<AGNI_VERSION_1,1>(ifrom, ito, thr);
      else eval<AGNI_VERSION_1,0>(ifrom, ito, thr);
    } else if (atomic_feature_version == AGNI_VERSION_2) {
      if (evflag) eval<AGNI_VERSION_2,1>(ifrom, ito, thr);
      else eval<AGNI_VERSION_2,0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int ATOMIC_FEATURE_VERSION, int EVFLAG>
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

    const Param &iparam = params[elem1param[itype]];
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
          double e = 0.0;

          if (ATOMIC_FEATURE_VERSION == AGNI_VERSION_1)
            e = fm_exp(-(iparam.eta[k]*rsq));
          else if (ATOMIC_FEATURE_VERSION == AGNI_VERSION_2)
            e = (1.0 / (square(iparam.eta[k]) * iparam.gwidth * sqrt(MathConst::MY_2PI)))
              * fm_exp(-(square(r - iparam.eta[k])) / (2.0 * square(iparam.gwidth)));

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

      for (int k = 0; k < iparam.numeta; ++k) {
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
