/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_FUNCTOR_COUL_LONG_CS_H
#define LMP_FUNCTOR_COUL_LONG_CS_H

// Core-shell (adiabatic Drude) variant of the long-range Coulomb policy.  It is
// derived from CoulLong and changes only the real-space Coulomb kernel, exactly
// as src/CORESHELL/pair_*_coul_long_cs.cpp derive from their base styles:
//
//   * rsq is bumped by a tiny epsilon (rsq_epsilon, applied by the driver) so an
//     exactly overlapping core/shell pair does not divide by zero;
//   * the direct-erfc branch uses a higher-accuracy polynomial approximation (the
//     B-coefficients below, with their own EWALD_P) that stays accurate as
//     r -> 0; and
//   * for bonded pairs (factor_coul < 1) the singular 1/r is regularized by
//     shifting r -> r + EPS_EWALD and using r2inv = 1/(rsq + EPS_EWALD_SQR).
//
// Everything else -- g_ewald, the interpolation tables, the cutoff, restart, and
// single() (which, as in the original /cs styles, uses the inherited A-polynomial
// kernel) -- is inherited unchanged from CoulLong.

#include "functor_coul_long.h"    // CoulLong (base), ewald_const.h, kspace.h

namespace LAMMPS_NS::functor {

// constants for the core-shell erfc approximation (see CORESHELL package).  These
// intentionally differ from EwaldConst: the B-polynomial uses its own EWALD_P.
namespace CoreShellConst {
  inline constexpr double EWALD_F = 1.12837917;
  inline constexpr double EWALD_P = 9.95473818e-1;
  inline constexpr double B0 = -0.1335096380159268;
  inline constexpr double B1 = -2.57839507e-1;
  inline constexpr double B2 = -1.37203639e-1;
  inline constexpr double B3 = -8.88822059e-3;
  inline constexpr double B4 = -5.80844129e-3;
  inline constexpr double B5 = 1.14652755e-1;
  inline constexpr double EPS_EWALD = 1.0e-6;
  inline constexpr double EPS_EWALD_SQR = 1.0e-12;
}    // namespace CoreShellConst

struct CoulLongCS : CoulLong {
  // tiny offset added to rsq by the driver (the "EPSILON" of the /cs styles)
  static constexpr double rsq_epsilon = 1.0e-20;

  template <bool EFLAG, int CTABLE>
  PairContribution eval_coul(double rsq, double qi, double qj, double factor_coul) const
  {
    using namespace CoreShellConst;
    const double qiqj = qi * qj;
    double forcecoul, ecoul = 0.0, prefactor = 0.0, r2inv;

    if (!CTABLE || rsq <= tabinnersq) {
      const double r = sqrt(rsq);
      prefactor = qqrd2e * qiqj;
      if (factor_coul < 1.0) {
        // bonded core/shell: regularize the singular 1/r with EPS_EWALD
        const double grij = g_ewald * (r + EPS_EWALD);
        const double expm2 = exp(-grij * grij);
        const double t = 1.0 / (1.0 + EWALD_P * grij);
        const double u = 1.0 - t;
        const double erfc = t * (1.0 + u * (B0 + u * (B1 + u * (B2 + u * (B3 + u * (B4 + u * B5)))))) * expm2;
        prefactor /= (r + EPS_EWALD);
        forcecoul = prefactor * (erfc + EWALD_F * grij * expm2 - (1.0 - factor_coul));
        r2inv = 1.0 / (rsq + EPS_EWALD_SQR);
        if constexpr (EFLAG) ecoul = prefactor * (erfc - (1.0 - factor_coul));
      } else {
        const double grij = g_ewald * r;
        const double expm2 = exp(-grij * grij);
        const double t = 1.0 / (1.0 + EWALD_P * grij);
        const double u = 1.0 - t;
        const double erfc = t * (1.0 + u * (B0 + u * (B1 + u * (B2 + u * (B3 + u * (B4 + u * B5)))))) * expm2;
        prefactor /= r;
        forcecoul = prefactor * (erfc + EWALD_F * grij * expm2);
        r2inv = 1.0 / rsq;
        if constexpr (EFLAG) ecoul = prefactor * erfc;
      }
    } else {
      // tabulated branch: identical to CoulLong (uses the rsq already bumped by
      // rsq_epsilon, matching the original /cs)
      Pair::union_int_float_t rsq_lookup;
      rsq_lookup.f = rsq;
      int itable = rsq_lookup.i & ncoulmask;
      itable >>= ncoulshiftbits;
      const double fraction = ((double) rsq_lookup.f - rtable[itable]) * drtable[itable];
      double table = ftable[itable] + fraction * dftable[itable];
      forcecoul = qiqj * table;
      if (factor_coul < 1.0) {
        table = ctable[itable] + fraction * dctable[itable];
        prefactor = qiqj * table;
        forcecoul -= (1.0 - factor_coul) * prefactor;
      }
      r2inv = 1.0 / rsq;
      if constexpr (EFLAG) {
        table = etable[itable] + fraction * detable[itable];
        ecoul = qiqj * table;
        if (factor_coul < 1.0) ecoul -= (1.0 - factor_coul) * prefactor;
      }
    }

    PairContribution out;
    out.fpair = forcecoul * r2inv;
    out.energy = ecoul;
    return out;
  }
};

}    // namespace LAMMPS_NS::functor

#endif
