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

#ifndef LMP_FUNCTOR_COUL_SOFT_H
#define LMP_FUNCTOR_COUL_SOFT_H

// Soft-core (FEP) Coulomb policies for the FUNCTOR driver.  Reimplement the
// Coulomb half of src/FEP/pair_coul_cut_soft.cpp and pair_coul_long_soft.cpp.
//
// The soft-core coupling parameter lambda is per type pair and is shared with the
// van der Waals term, so it is owned by the *evaluator* and the two derived
// quantities the Coulomb kernel needs are read straight from the evaluator's
// Param (passed in by the driver): p.lam1 = lambda^nlambda (the prefactor) and
// p.lam2 = alphac * (1-lambda)^2 (the soft-core denominator offset).  Both soft
// evaluators (EvaluatorLJCutSoft, EvaluatorNoneSoft) expose those two members.
//
// The screened distance is denc = sqrt(p.lam2 + rsq); the force is
// qqrd2e*lam1*qi*qj/denc^3 and the energy qqrd2e*lam1*qi*qj/denc, both already in
// "force/r" form (no extra 1/r^2 factor), matching the driver's convention.  The
// long-range variant multiplies in the Ewald real-space correction, evaluated at
// the *true* distance r = sqrt(rsq).  Soft-core styles never use the Coulomb
// interpolation tables (has_table is false).

#include "functor_coul_policies.h"    // CoulCutoffBase, PairContribution

#include "ewald_const.h"
#include "kspace.h"

namespace LAMMPS_NS::functor {

// cutoff soft-core Coulomb (coul/cut/soft and .../coul/cut/soft)

struct CoulCutSoft : CoulCutoffBase {
  static constexpr bool has_coul = true;
  static constexpr bool has_table = false;
  static constexpr bool needs_charge = true;
  static constexpr bool needs_kspace = false;
  static constexpr bool per_pair_cutoff = true;

  void init_style(Pair *, LAMMPS *lmp) { qqrd2e = lmp->force->qqrd2e; }

  template <bool EFLAG, int /*CTABLE*/, class P>
  PairContribution eval_coul(double rsq, double qi, double qj, double factor_coul,
                             const P &p) const
  {
    const double denc = sqrt(p.lam2 + rsq);
    const double qiqj = qi * qj;
    const double forcecoul = qqrd2e * p.lam1 * qiqj / (denc * denc * denc);
    PairContribution out;
    out.fpair = factor_coul * forcecoul;
    if constexpr (EFLAG)
      out.energy = factor_coul * qqrd2e * p.lam1 * qiqj / denc;
    else
      out.energy = 0.0;
    return out;
  }

  template <class P>
  PairContribution single_coul(LAMMPS *lmp, int i, int j, double rsq, double factor_coul,
                               const P &p) const
  {
    const double *q = lmp->atom->q;
    const double denc = sqrt(p.lam2 + rsq);
    const double qiqj = q[i] * q[j];
    const double forcecoul = qqrd2e * p.lam1 * qiqj / (denc * denc * denc);
    PairContribution out;
    out.fpair = factor_coul * forcecoul;
    out.energy = factor_coul * qqrd2e * p.lam1 * qiqj / denc;
    return out;
  }
};

// long-range (Ewald/PPPM real space) soft-core Coulomb (coul/long/soft and
// .../coul/long/soft).  Global cutoff, no interpolation tables.

struct CoulLongSoft : CoulCutoffBase {
  static constexpr bool has_coul = true;
  static constexpr bool has_table = false;
  static constexpr bool needs_charge = true;
  static constexpr bool needs_kspace = true;
  static constexpr bool per_pair_cutoff = false;    // single global Coulomb cutoff

  double g_ewald = 0.0;

  void init_style(Pair *, LAMMPS *lmp)
  {
    qqrd2e = lmp->force->qqrd2e;
    if (lmp->force->kspace == nullptr) lmp->error->all(FLERR, "Pair style requires a KSpace style");
    g_ewald = lmp->force->kspace->g_ewald;
  }

  template <bool EFLAG, int /*CTABLE*/, class P>
  PairContribution eval_coul(double rsq, double qi, double qj, double factor_coul,
                             const P &p) const
  {
    using namespace EwaldConst;
    const double r = sqrt(rsq);
    const double grij = g_ewald * r;
    const double expm2 = exp(-grij * grij);
    const double t = 1.0 / (1.0 + EWALD_P * grij);
    const double erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;

    const double denc = sqrt(p.lam2 + rsq);
    const double qiqj = qi * qj;
    const double fprefactor = qqrd2e * p.lam1 * qiqj / (denc * denc * denc);
    double forcecoul = fprefactor * (erfc + EWALD_F * grij * expm2);
    if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * fprefactor;

    PairContribution out;
    out.fpair = forcecoul;
    if constexpr (EFLAG) {
      const double eprefactor = qqrd2e * p.lam1 * qiqj / denc;
      double ecoul = eprefactor * erfc;
      if (factor_coul < 1.0) ecoul -= (1.0 - factor_coul) * eprefactor;
      out.energy = ecoul;
    } else
      out.energy = 0.0;
    return out;
  }

  template <class P>
  PairContribution single_coul(LAMMPS *lmp, int i, int j, double rsq, double factor_coul,
                               const P &p) const
  {
    using namespace EwaldConst;
    const double *q = lmp->atom->q;
    const double r = sqrt(rsq);
    const double grij = g_ewald * r;
    const double expm2 = exp(-grij * grij);
    const double t = 1.0 / (1.0 + EWALD_P * grij);
    const double erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;

    const double denc = sqrt(p.lam2 + rsq);
    const double qiqj = q[i] * q[j];
    const double fprefactor = qqrd2e * p.lam1 * qiqj / (denc * denc * denc);
    double forcecoul = fprefactor * (erfc + EWALD_F * grij * expm2);
    if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * fprefactor;

    const double eprefactor = qqrd2e * p.lam1 * qiqj / denc;
    double phicoul = eprefactor * erfc;
    if (factor_coul < 1.0) phicoul -= (1.0 - factor_coul) * eprefactor;

    PairContribution out;
    out.fpair = forcecoul;
    out.energy = phicoul;
    return out;
  }
};

}    // namespace LAMMPS_NS::functor

#endif
