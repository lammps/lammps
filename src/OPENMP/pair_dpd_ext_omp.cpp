// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_dpd_ext_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "update.h"
#include "random_mars.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;

#define EPSILON 1.0e-10

/* ---------------------------------------------------------------------- */

PairDPDExtOMP::PairDPDExtOMP(LAMMPS *lmp) :
  PairDPDExt(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
  random_thr = nullptr;
  nthreads = 0;
}

/* ---------------------------------------------------------------------- */

PairDPDExtOMP::~PairDPDExtOMP()
{
  if (random_thr) {
    for (int i=1; i < nthreads; ++i)
      delete random_thr[i];

    delete[] random_thr;
    random_thr = nullptr;
  }
}

/* ---------------------------------------------------------------------- */

void PairDPDExtOMP::compute(int eflag, int vflag)
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

    // to ensure full compatibility with the serial DPD style
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
    if ((tid > 0) && (random_thr[tid] == nullptr))
      random_thr[tid] = new RanMars(Pair::lmp, seed + comm->me
                                    + comm->nprocs*tid);

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
void PairDPDExtOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpairx,fpairy,fpairz,fpair;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  double rsq,r,rinv,dot,wd,wdPar,wdPerp,randnum,randnumx,randnumy,randnumz,factor_dpd,factor_sqrt;
  double P[3][3];
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  const auto * _noalias const v = (dbl3_t *) atom->v[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
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
    fxtmp=fytmp=fztmp=0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_dpd = special_lj[sbmask(j)];
      factor_sqrt = special_sqrt[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        if (r < EPSILON) continue;     // r can be 0.0 in DPD systems
        rinv = 1.0/r;
        delvx = vxtmp - v[j].x;
        delvy = vytmp - v[j].y;
        delvz = vztmp - v[j].z;
        dot = delx*delvx + dely*delvy + delz*delvz;

        P[0][0] = 1.0 - delx*delx*rinv*rinv;
        P[0][1] =     - delx*dely*rinv*rinv;
        P[0][2] =     - delx*delz*rinv*rinv;

        P[1][0] = P[0][1];
        P[1][1] = 1.0 - dely*dely*rinv*rinv;
        P[1][2] =     - dely*delz*rinv*rinv;

        P[2][0] = P[0][2];
        P[2][1] = P[1][2];
        P[2][2] = 1.0 - delz*delz*rinv*rinv;

        wd = 1.0 - r/cut[itype][jtype];
        wdPar = pow(wd,ws[itype][jtype]);
        wdPerp = pow(wd,wsT[itype][jtype]);

        randnum = rng.gaussian();
        randnumx = rng.gaussian();
        randnumy = rng.gaussian();
        randnumz = rng.gaussian();

        // conservative force
        fpair = a0[itype][jtype]*wd;

        // drag force - parallel
        fpair -= gamma[itype][jtype]*wdPar*wdPar*dot*rinv;
        fpair *= factor_dpd;

        // random force - parallel
        fpair += factor_sqrt*sigma[itype][jtype]*wdPar*randnum*dtinvsqrt;

        fpairx = fpair*rinv*delx;
        fpairy = fpair*rinv*dely;
        fpairz = fpair*rinv*delz;

        // drag force - perpendicular
        const double prefactor_g = factor_dpd * gammaT[itype][jtype]*wdPerp*wdPerp;
        fpairx -= prefactor_g * (P[0][0]*delvx + P[0][1]*delvy + P[0][2]*delvz);
        fpairy -= prefactor_g * (P[1][0]*delvx + P[1][1]*delvy + P[1][2]*delvz);
        fpairz -= prefactor_g * (P[2][0]*delvx + P[2][1]*delvy + P[2][2]*delvz);

        // random force - perpendicular
        const double prefactor_s = factor_sqrt * sigmaT[itype][jtype]*wdPerp * dtinvsqrt;
        fpairx += prefactor_s * (P[0][0]*randnumx + P[0][1]*randnumy + P[0][2]*randnumz);
        fpairy += prefactor_s * (P[1][0]*randnumx + P[1][1]*randnumy + P[1][2]*randnumz);
        fpairz += prefactor_s * (P[2][0]*randnumx + P[2][1]*randnumy + P[2][2]*randnumz);

        fxtmp += fpairx;
        fytmp += fpairy;
        fztmp += fpairz;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= fpairx;
          f[j].y -= fpairy;
          f[j].z -= fpairz;
        }

        if (EFLAG) {
          // unshifted eng of conservative term:
          // evdwl = -a0[itype][jtype]*r * (1.0-0.5*r/cut[itype][jtype]);
          // eng shifted to 0.0 at cutoff
          evdwl = 0.5*a0[itype][jtype]*cut[itype][jtype] * wd*wd;
          evdwl *= factor_dpd;
        }

        if (EVFLAG) ev_tally_xyz_thr(this, i,j,nlocal,NEWTON_PAIR,evdwl,0.0,
                                     fpairx,fpairy,fpairz,delx,dely,delz,thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairDPDExtOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairDPDExt::memory_usage();
  bytes += (double)comm->nthreads * sizeof(RanMars*);
  bytes += (double)comm->nthreads * sizeof(RanMars);

  return bytes;
}
