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
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_special.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"
#include "update.h"   // needed for Update::ntimestep

#include <algorithm>
#include <cmath>
#include <cstdlib>  // for setenv
#include <cstring>
#include <unordered_map>
#include <vector>

#include "omp_compat.h"

using namespace LAMMPS_NS;

// global ad hoc parameters - copied from pair_dispersion_d3.cpp
static constexpr double K1 = 16.0;
static constexpr double K3 = -4.0;

static constexpr double AUTOANG = 0.52917725;    // atomic units (Bohr) to Angstrom
static constexpr double AUTOEV = 27.21140795;    // atomic units (Hartree) to eV

/* ---------------------------------------------------------------------- */

PairDispersionD3OMP::PairDispersionD3OMP(LAMMPS *lmp) :
  PairDispersionD3(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void PairDispersionD3OMP::calc_coordination_number()
{
  int nlocal = atom->nlocal;
  const int nthreads = comm->nthreads;
  int nall = nlocal + atom->nghost;

  int newton_pair = force->newton_pair;

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->grow(cn, nmax, "pair:cn");
    memory->grow(dc6, nmax, "pair:dc6");
  }

  // zero out coordination number
  memset(cn, 0, sizeof(double) * (newton_pair ? nall : nlocal));
  memset(dc6, 0, sizeof(double) * (newton_pair ? nall : nlocal));

  int inum = list->inum;

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
  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  const int * const * const firstneigh = list->firstneigh;

  // Thread-local cn array to avoid race conditions
  auto thr_cn = std::vector<double>(atom->nmax);  // Initialize to zero

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

  // Contribute thread-local cn to global cn
  for (int i = 0; i < atom->nmax; i++) {
    // No need for atomic if newton_pair is false
    #pragma omp atomic
    cn[i] += thr_cn[i];
  }
}

/* ---------------------------------------------------------------------- */

void PairDispersionD3OMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);
  // First call coordination number calculation
  calc_coordination_number();

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

  // Zero out dc6 values before OpenMP section
  memset(dc6, 0, sizeof(double) * nall);

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
    reduce_thr(this, eflag, vflag_either, thr);
  } //end of omp parallel region

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   Modified from serial code to avoid race conditions
------------------------------------------------------------------------- */

void PairDispersionD3OMP::get_dC6(int iat, int jat, double cni, double cnj, double c6_res[3])
{
  double c6_ref, cni_ref, cnj_ref;
  double c6mem, r_save, r;
  double expterm, term;
  double num, den, d_num_i, d_num_j, d_den_i, d_den_j;

  c6mem = -1.0e20, r_save = 1.0e20;
  num = 0;
  den = 0;
  d_num_i = 0;
  d_num_j = 0;
  d_den_i = 0;
  d_den_j = 0;

  for (int ci = 0; ci <= mxci[iat]; ci++) {
    for (int cj = 0; cj <= mxci[jat]; cj++) {

      c6_ref = c6ab[iat][jat][ci][cj][0];
      double autoang6 = AUTOANG * AUTOANG * AUTOANG;
      autoang6 = MathSpecial::square(autoang6);
      c6_ref *= AUTOEV * autoang6;

      if (c6_ref > 0) {
        cni_ref = c6ab[iat][jat][ci][cj][1];
        cnj_ref = c6ab[iat][jat][ci][cj][2];

        r = (cni - cni_ref) * (cni - cni_ref) + (cnj - cnj_ref) * (cnj - cnj_ref);

        if (r < r_save) {
          r_save = r;
          c6mem = c6_ref;
        }

        expterm = exp(static_cast<double>(K3) * static_cast<double>(r));

        num += c6_ref * expterm;
        den += expterm;

        expterm = expterm * 2.0 * K3;

        term = expterm * (cni - cni_ref);
        d_num_i += c6_ref * term;
        d_den_i += term;

        term = expterm * (cnj - cnj_ref);
        d_num_j += c6_ref * term;
        d_den_j += term;
      }
    }
  }

  if (den > 1.0E-99) {
    c6_res[0] = num / den;
    c6_res[1] = ((d_num_i * den) - (d_den_i * num)) / (den * den);
    c6_res[2] = ((d_num_j * den) - (d_den_j * num)) / (den * den);
  } else {
    c6_res[0] = c6mem;
    c6_res[1] = 0;
    c6_res[2] = 0;
  }
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
  double evdwl = 0.0;

  // Thread-local dc6 array to avoid race conditions
  auto thr_dc6 = std::vector<double>(atom->nmax);  // Initialize to zero.

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

        // Modified from original code to avoid race conditions
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

          case 3: {    // bj
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

          case 4: {    // bjm
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

          default: {
            error->all(FLERR, Error::NOLASTLINE, "Damping code {} unknown", dampingCode);
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

  // Contribute thread-local dc6 to global dc6
  for (int i = 0; i < atom->nmax; i++) {
    // Possibly no need for atomic if newton_pair is false
    #pragma omp atomic
    dc6[i] += thr_dc6[i];
  }
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

  double dc6tmp,xtmp,ytmp,ztmp,delx,dely,delz,rsq,factor_lj,dcn,rcovij,expterm,fpair,fxtmp,fytmp,fztmp,r;

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
  return bytes;
}
