// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Chun-Yeol You (DGIST)
------------------------------------------------------------------------- */

#include "pair_mdpd_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "random_mars.h"
#include "suffix.h"
#include "update.h"

#include <cmath>

#include "omp_compat.h"

using namespace LAMMPS_NS;

static constexpr double EPSILON = 1.0e-10;

/* ---------------------------------------------------------------------- */

PairMDPDOMP::PairMDPDOMP(LAMMPS *lmp) :
  PairMDPD(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
  random_thr = nullptr;
  nthreads = 0;
}

/* ---------------------------------------------------------------------- */

PairMDPDOMP::~PairMDPDOMP()
{
  if (random_thr) {
    for (int i=1; i < nthreads; ++i)
      delete random_thr[i];

    delete[] random_thr;
    random_thr = nullptr;
  }
}

/* ----------------------------------------------------------------------
   the base class settings() re-creates the serial random number generator,
   so the pool of per-thread generators must be re-created on the next call
------------------------------------------------------------------------- */

void PairMDPDOMP::settings(int narg, char **arg)
{
  PairMDPD::settings(narg, arg);

  if (random_thr) {
    for (int i = 1; i < nthreads; ++i) delete random_thr[i];
    delete[] random_thr;
    random_thr = nullptr;
  }
  nthreads = 0;
}

/* ---------------------------------------------------------------------- */

void PairMDPDOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int inum = list->inum;

  // number of threads has changed. reallocate pool of pRNGs
  if (nthreads != comm->nthreads) {
    if (random_thr) {
      for (int i=1; i < nthreads; ++i)
        delete random_thr[i];

      delete[] random_thr;
    }

    nthreads = comm->nthreads;
    random_thr = new RanMars*[nthreads];
    for (int i=1; i < nthreads; ++i)
      random_thr[i] = nullptr;

    // to ensure full compatibility with the serial MDPD style
    // we use the serial random number generator instance for thread 0
    random_thr[0] = random;
  }

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    // generate a random number generator instance for
    // all threads != 0. make sure we use unique seeds.
    if ((tid > 0) && (random_thr[tid] == nullptr)) {
      // wrap around to stay within the valid seed range of RanMars
      const bigint tseed = (bigint) seed + comm->me + (bigint) comm->nprocs * tid;
      random_thr[tid] = new RanMars(Pair::lmp, (int) ((tseed - 1) % 900000000 + 1));
    }

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
void PairMDPDOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  double rsq,r,rinv,dot,wc,wc_r,wr,randnum,factor_dpd;
  double rhoi,rhoj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  const auto * _noalias const v = (dbl3_t *) atom->v[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const double * _noalias const rho = atom->rho;
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double *special_lj = force->special_lj;
  const double dtinvsqrt = 1.0/sqrt(update->dt);
  double fxtmp,fytmp,fztmp;
  RanMars &rng = *random_thr[thr->get_tid()];

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    vxtmp = v[i].x;
    vytmp = v[i].y;
    vztmp = v[i].z;
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    rhoi = rho[i];
    fxtmp=fytmp=fztmp=0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_dpd = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        if (r < EPSILON) continue;     // r can be 0.0 in MDPD systems
        rinv = 1.0/r;
        delvx = vxtmp - v[j].x;
        delvy = vytmp - v[j].y;
        delvz = vztmp - v[j].z;
        dot = delx*delvx + dely*delvy + delz*delvz;

        wc = 1.0 - r/cut[itype][jtype];
        wc_r = 1.0 - r/cut_r[itype][jtype];
        wc_r = MAX(wc_r,0.0);
        wr = wc;

        rhoj = rho[j];
        randnum = rng.gaussian();

        // conservative force = A_att * wc + B_rep*(rhoi+rhoj)*wc_r
        // drag force = -gamma * wr^2 * (delx dot delv) / r
        // random force = sigma * wr * rnd * dtinvsqrt;

        fpair = A_att[itype][jtype]*wc + B_rep[itype][jtype]*(rhoi+rhoj)*wc_r;
        fpair -= gamma[itype][jtype]*wr*wr*dot*rinv;
        fpair += sigma[itype][jtype]*wr*randnum*dtinvsqrt;
        fpair *= factor_dpd*rinv;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        if (EFLAG) {
          // unshifted eng of conservative term:
          // eng shifted to 0.0 at cutoff
          evdwl = 0.5*A_att[itype][jtype]*cut[itype][jtype] * wr*wr
                + 0.5*B_rep[itype][jtype]*cut_r[itype][jtype]*(rhoi+rhoj)*wc_r*wc_r;
          evdwl *= factor_dpd;
        }

        if (EVFLAG) ev_tally_thr(this,i,j,nlocal,NEWTON_PAIR,
                                 evdwl,0.0,fpair,delx,dely,delz,thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairMDPDOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairMDPD::memory_usage();
  bytes += (double)comm->nthreads * sizeof(RanMars*);
  bytes += (double)comm->nthreads * sizeof(RanMars);

  return bytes;
}
