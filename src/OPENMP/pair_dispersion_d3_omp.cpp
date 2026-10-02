// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
      Ulrik Unneberg
      Marc L. Descoteaux
      Yizhong R. Hu
      William C. Witt
      Affiliation: Harvard University
------------------------------------------------------------------------- */

#include "pair_dispersion_d3_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_special.h"
#include "memory.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>
#include <cstring>

#include "omp_compat.h"

using namespace LAMMPS_NS;
using namespace LAMMPS_NS::DispersionD3;

/* ---------------------------------------------------------------------- */

PairDispersionD3OMP::PairDispersionD3OMP(LAMMPS *lmp) :
  PairDispersionD3(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void PairDispersionD3OMP::calc_coordination_number()
{
  const int nthreads = comm->nthreads;
  const int newton_pair = force->newton_pair;

  // cn and dc6 hold one copy per thread, each nall long.  The threads only
  // ever touch their own copy; data_reduce_thr() sums the copies into the
  // first nall elements at the end of the respective loops.

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->grow(cn, nthreads * nmax, "pair:cn");
    memory->grow(dc6, nthreads * nmax, "pair:dc6");
  }

  const int inum = list->inum;

  // Begin parallel region, the central atoms indexed by ii are assigned to different threads.
  #if defined(_OPENMP)
  #pragma omp parallel LMP_DEFAULT_NONE \
  firstprivate(inum,nthreads)
  #endif
  {
    int ifrom, ito, tid;

    // Set up the starting and ending indices for each thread
    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);

    // Calculate coordination number with the helper functions
    // The flags need to be constants for the template instantiation
    if (force->newton_pair) {
        eval_coordination<1>(ifrom,ito,thr);
    } else{
        eval_coordination<0>(ifrom,ito,thr);
    }

    thr->timer(Timer::PAIR);
  }

  // communicate coordination number
  communicationStage = 1;
  if (newton_pair) comm->reverse_comm(this);
  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

template <int NEWTON_PAIR>
void PairDispersionD3OMP::eval_coordination(int iifrom, int iito, ThrData * const thr)
{

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int tid = thr->get_tid();
  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  const int * const * const firstneigh = list->firstneigh;

  // this thread's private copy of the coordination number accumulator
  double * _noalias const thr_cn = cn + tid * nall;
  memset(thr_cn, 0, sizeof(double) * nall);

  for (int ii = iifrom; ii < iito; ii++) {

    int i = ilist[ii];
    int itype = type[i];
    const int * _noalias const jlist = firstneigh[i];
    int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {

      int j = jlist[jj];
      j &= NEIGHMASK;
      int jtype = type[j];

      double delrj[3];
      delrj[0] = x[i].x - x[j].x;
      delrj[1] = x[i].y - x[j].y;
      delrj[2] = x[i].z - x[j].z;

      double rsq = delrj[0] * delrj[0] + delrj[1] * delrj[1] + delrj[2] * delrj[2];

      // if the atoms are too far away don't consider the contribution
      if (rsq > cn_thr) continue;

      double rr = sqrt(rsq);
      double rcov_ij = (rcov[itype] + rcov[jtype]) * AUTOANG;
      double cn_ij = 1.0 / (1.0 + exp(-K1 * ((rcov_ij / rr) - 1.0)));

      // update coordination number on a thread-local array
      thr_cn[i] += cn_ij;
      if (NEWTON_PAIR || j < nlocal) { thr_cn[j] += cn_ij; }
    }
  }

  // sum the per thread copies into cn[0] ... cn[nall-1]

  sync_threads();
  data_reduce_thr(cn, nall, nthreads, 1, tid);
  sync_threads();
}

/* ---------------------------------------------------------------------- */

void PairDispersionD3OMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  // dampingCode is validated in PairDispersionD3::init_style(), so the loops
  // below need no default case (error->all() must not be called from inside a
  // parallel region)

  // First call coordination number calculation
  calc_coordination_number();

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

// Parallel direct force computation and some other quantities calculation.
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag) \
firstprivate(inum,nthreads,nall)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    // Call the helper eval function with the appropriate flags for the first phase of the computation
    // Again, the flags need to be constants for the template instantiation
    if (evflag) {
      if (eflag) {
        if (force->newton_pair) {
          eval_first_phase<1,1,1>(ifrom, ito, thr);
        }
        else {
          eval_first_phase<1,1,0>(ifrom, ito, thr);
        }
      } else {
        if (force->newton_pair) eval_first_phase<1,0,1>(ifrom, ito, thr);
        else eval_first_phase<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair) eval_first_phase<0,0,1>(ifrom, ito, thr);
      else eval_first_phase<0,0,0>(ifrom, ito, thr);
    }
    thr->timer(Timer::PAIR);
  } // end of omp parallel region

  // Both phases tally into the same ThrData, so ev_setup_thr() is called only
  // in the first region (it zeroes the per thread accumulators) and
  // reduce_thr() only at the end of the second one.  The per thread force
  // arrays therefore stay unreduced across the communication below, which is
  // safe because it only exchanges dc6.

  // Communication stage 2 for dc6 values in preparation for calculation of indirect forces in the second phase
  communicationStage = 2;
  if (force->newton_pair) {
    comm->reverse_comm(this);
  }

  comm->forward_comm(this);

  // Process the second phase with the combined dc6 values
  #if defined(_OPENMP)
  #pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag) \
  firstprivate(inum,nthreads,nall)
  #endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);

    // Call the helper eval function with the appropriate flags for the second phase of the computation
    // Again, the flags need to be constants for the template instantiation
    if (evflag) {
      if (eflag) {
        if (force->newton_pair) {
          eval_second_phase<1,1,1>(ifrom, ito, thr);
        }
        else {
          eval_second_phase<1,1,0>(ifrom, ito, thr);
        }
      } else {
        if (force->newton_pair) eval_second_phase<1,0,1>(ifrom, ito, thr);
        else eval_second_phase<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair) eval_second_phase<0,0,1>(ifrom, ito, thr);
      else eval_second_phase<0,0,0>(ifrom, ito, thr);
    }
    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } //end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairDispersionD3OMP::eval_first_phase(int iifrom, int iito, ThrData * const thr)
{
  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * _noalias const special_lj = force->special_lj;
  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  const int * const * const firstneigh = list->firstneigh;
  const int nall = nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int tid = thr->get_tid();
  double evdwl = 0.0;

  // this thread's private copy of the dE/dC6 accumulator
  double * _noalias const thr_dc6 = dc6 + tid * nall;
  memset(thr_dc6, 0, sizeof(double) * nall);

  // Loop over assigned atoms
  for (int ii = iifrom; ii < iito; ++ii) {
    int i = ilist[ii];

    double xtmp = x[i].x;
    double ytmp = x[i].y;
    double ztmp = x[i].z;
    int itype = type[i];
    int jnum = numneigh[i];
    const int * _noalias const jlist = firstneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      double factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      double delx = xtmp - x[j].x;
      double dely = ytmp - x[j].y;
      double delz = ztmp - x[j].z;

      double rsq = delx * delx + dely * dely + delz * delz;

      int jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {

        double r2inv = 1.0 / rsq;
        double r6inv = r2inv * r2inv * r2inv;
        double r8inv = r2inv * r2inv * r2inv * r2inv;
        double r10inv = r2inv * r2inv * r2inv * r2inv * r2inv;

        // get_dC6 writes {C6, dC6/dCN_i, dC6/dCN_j}
        double c6_res[3] = {};
        get_dC6(itype, jtype, cn[i], cn[j], c6_res);

        double C6 = c6_res[0];
        double C8 = 3.0 * C6 * r2r4[itype] * r2r4[jtype] * AUTOANG * AUTOANG;

        double alpha6 = alpha;
        double alpha8 = alpha + 2;

        double t6, t8, damp6, damp8, e6, e8;
        double tmp6, tmp8, fpair1, fpair2, fpair;
        t6 = t8 = e6 = e8 = evdwl = fpair = fpair1 = fpair2 = 0.0;

        // Damping code selection - now using the passed dampingCode parameter
        switch (dampingCode) {
          // Written to avoid using sqrt and pow()
          case 1: /* Original damping */
            {
              double ip6 = rs6 * r0ab[type[i]][type[j]];
              double ip8 = rs8 * r0ab[type[i]][type[j]];

              double half_alpha6 = 0.5 * alpha6;
              double half_alpha8 = 0.5 * alpha8;

              t6 = MathSpecial::powauto(ip6, alpha6) * MathSpecial::powauto(rsq, -half_alpha6);
              t8 = MathSpecial::powauto(ip8, alpha8) * MathSpecial::powauto(rsq, -half_alpha8);

              damp6 = 1.0 / (1.0 + 6.0 * t6);
              damp8 = 1.0 / (1.0 + 6.0 * t8);

              e6 = C6 * damp6 * r6inv;
              e8 = C8 * damp8 * r8inv;

              tmp6 = 6 * s6 * C6 * r8inv * damp6;
              tmp8 = 8 * s8 * C8 * r10inv * damp8;

              fpair1 = -tmp6 - tmp8;
              fpair2 = tmp6 * alpha6 * t6 * damp6 + (3.0 / 4.0) * tmp8 * alpha8 * t8 * damp8;

              fpair = fpair1 + fpair2;
              fpair *= factor_lj;
            } break;
          // Written to avoid pow
          case 2: {    // zerom

            double r = sqrt(rsq);
            double r0 = r0ab[type[i]][type[j]];

            t6 = MathSpecial::powauto((r / (rs6 * r0)) + rs8 * r0, -alpha6);
            damp6 = 1.0 / (1.0 + 6.0 * t6);
            t8 = MathSpecial::powauto((r / r0) + rs8 * r0, -alpha8);
            damp8 = 1.0 / (1.0 + 6.0 * t8);

            e6 = C6 * damp6 * r6inv;
            e8 = C8 * damp8 * r8inv;

            tmp6 = 6 * s6 * C6 * r8inv * damp6;
            tmp8 = 8 * s8 * C8 * r10inv * damp8;

            fpair1 = -tmp6 - tmp8;

            double fp26 = tmp6 * alpha6 * t6 * damp6 * r / (r + rs6 * rs8 * r0 * r0);
            double fp28 = tmp8 * alpha8 * t8 * damp8 * r / (r + rs8 * r0 * r0);

            fpair2 = fp26 + (3.0 / 4.0) * fp28;

            fpair = fpair1 + fpair2;
            fpair *= factor_lj;
          } break;

          case 3:      // bj
          case 4: {    // bjm, same functional form as bj, different parameters
            double r0 = sqrt(C8 / C6);

            double r4 = rsq * rsq;
            double r6 = rsq * rsq * rsq;
            double r8 = rsq * rsq * rsq * rsq;

            double d = a1 * r0 + a2;
            double d2 = d * d;
            double d4 = d2 * d2;

            t6 = r6 + MathSpecial::cube(d2);
            t8 = r8 + MathSpecial::square(d4);

            e6 = C6 / t6;
            e8 = C8 / t8;

            tmp6 = 6.0 * s6 * C6 * r4 / (t6 * t6);
            tmp8 = 8.0 * s8 * C8 * r6 / (t8 * t8);

            fpair = -(tmp6 + tmp8);
            fpair *= factor_lj;
          } break;
        }

        if (EFLAG) evdwl = -(s6 * e6 + s8 * e8) * factor_lj;

        double rest = (s6 * e6 + s8 * e8) / C6;

        // Update thread-local dc6
        double dc6_contrib_i = rest * c6_res[1];
        thr_dc6[i] += dc6_contrib_i;

        if (NEWTON_PAIR || j < nlocal) {
          double dc6_contrib_j = rest * c6_res[2];
          thr_dc6[j] += dc6_contrib_j;
        }

        f[i].x += delx * fpair;
        f[i].y += dely * fpair;
        f[i].z += delz * fpair;

        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx * fpair;
          f[j].y -= dely * fpair;
          f[j].z -= delz * fpair;
        }

        // Update energy and virial
        if (EVFLAG) ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR, evdwl, 0.0, fpair, delx, dely, delz, thr);
      }
    }
  }

  // sum the per thread copies into dc6[0] ... dc6[nall-1]

  sync_threads();
  data_reduce_thr(dc6, nall, nthreads, 1, tid);
  sync_threads();
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairDispersionD3OMP::eval_second_phase(int iifrom, int iito, ThrData * const thr)
{
  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * _noalias const special_lj = force->special_lj;
  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  const int * const * const firstneigh = list->firstneigh;

  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,factor_lj,dcn,rcovij,expterm,fpair,fxtmp,fytmp,fztmp,r;

  // Loop over assigned center atoms
  for (int ii = iifrom; ii < iito; ii++) {
    int i = ilist[ii];
    int itype = type[i];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;

    int jnum = numneigh[i];
    const int * _noalias const jlist = firstneigh[i];
    fxtmp=fytmp=fztmp=0.0;
    // Neighbor atom
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;

      rsq = delx * delx + dely * dely + delz * delz;
      int jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);

        if (rsq < cn_thr) {
          rcovij = (rcov[type[i]] + rcov[type[j]]) * AUTOANG;
          expterm = exp(-K1 * (rcovij / r - 1.0));
          dcn = -K1 * rcovij * expterm / (rsq * (expterm + 1.0) * (expterm + 1.0));

        } else {
          dcn = 0.0;
        }

        fpair = dcn * (dc6[i] + dc6[j]) / r;
        fpair *= factor_lj;

        fxtmp += delx * fpair;
        fytmp += dely * fpair;
        fztmp += delz * fpair;
        if (NEWTON_PAIR || j < nlocal) {
            f[j].x -= delx * fpair;
            f[j].y -= dely * fpair;
            f[j].z -= delz * fpair;
        }

        // Update virial (no energy contributions in this phase)
        if (EVFLAG) ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR, 0.0, 0.0, fpair, delx, dely, delz, thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairDispersionD3OMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairDispersionD3::memory_usage();
  // cn and dc6 hold comm->nthreads copies here, the base class counts one each
  bytes += (double) (comm->nthreads - 1) * nmax * 2 * sizeof(double);
  return bytes;
}
