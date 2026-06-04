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

#ifndef LMP_FUNCTOR_COUL_LONG_H
#define LMP_FUNCTOR_COUL_LONG_H

// Long-range (Ewald/PPPM) real-space Coulomb policy for the FUNCTOR driver.
// Reimplements the Coulomb half of src/KSPACE/pair_lj_cut_coul_long.cpp.
//
// Only core headers are needed (the base KSpace class and g_ewald live in
// src/kspace.h, the Ewald constants in src/ewald_const.h, and init_tables() is a
// core Pair method), so this compiles without the KSPACE package; at run time it
// requires a kspace_style, which is only available when KSPACE is installed.

#include "functor_coul_policies.h"    // CoulCutoffBase, PairContribution

#include "ewald_const.h"
#include "kspace.h"

namespace LAMMPS_NS::functor {

// The Coulomb cutoff is global (the per-pair arrays from CoulCutoffBase are
// filled with the single global value for uniformity with the driver).  The
// real-space term uses the direct erfc polynomial for rsq <= tabinnersq and the
// bitmapped interpolation tables otherwise; the CTABLE template parameter
// selects whether tables are used at all (driver dispatch on ncoultablebits).

struct CoulLong : CoulCutoffBase {
  static constexpr bool has_coul = true;
  static constexpr bool has_table = true;
  static constexpr bool needs_charge = true;
  static constexpr bool needs_kspace = true;
  static constexpr bool per_pair_cutoff = false;    // single global Coulomb cutoff

  double g_ewald = 0.0;

  // Coulomb interpolation tables: pointers into the owning Pair, valid after
  // init_style() (which builds them via Pair::init_tables when ncoultablebits>0)
  int ncoultablebits = 0;
  int ncoulmask = 0, ncoulshiftbits = 0;
  double tabinnersq = 0.0;
  double *rtable = nullptr, *drtable = nullptr, *ftable = nullptr, *dftable = nullptr;
  double *ctable = nullptr, *dctable = nullptr, *etable = nullptr, *detable = nullptr;

  void init_style(Pair *p, LAMMPS *lmp)
  {
    qqrd2e = lmp->force->qqrd2e;
    if (lmp->force->kspace == nullptr) lmp->error->all(FLERR, "Pair style requires a KSpace style");
    g_ewald = lmp->force->kspace->g_ewald;

    ncoultablebits = p->ncoultablebits;
    if (ncoultablebits) {
      p->init_tables(cut_coul_global, nullptr);
      tabinnersq = p->tabinnersq;
      ncoulmask = p->ncoulmask;
      ncoulshiftbits = p->ncoulshiftbits;
      rtable = p->rtable;
      drtable = p->drtable;
      ftable = p->ftable;
      dftable = p->dftable;
      ctable = p->ctable;
      dctable = p->dctable;
      etable = p->etable;
      detable = p->detable;
    }
  }

  template <bool EFLAG, int CTABLE>
  PairContribution eval_coul(double rsq, double qi, double qj, double factor_coul) const
  {
    using namespace EwaldConst;
    const double r2inv = 1.0 / rsq;
    const double qiqj = qi * qj;
    double forcecoul, ecoul = 0.0, prefactor = 0.0;

    if (!CTABLE || rsq <= tabinnersq) {
      const double r = sqrt(rsq);
      const double grij = g_ewald * r;
      const double expm2 = exp(-grij * grij);
      const double t = 1.0 / (1.0 + EWALD_P * grij);
      const double erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
      prefactor = qqrd2e * qiqj / r;
      forcecoul = prefactor * (erfc + EWALD_F * grij * expm2);
      if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor;
      if constexpr (EFLAG) {
        ecoul = prefactor * erfc;
        if (factor_coul < 1.0) ecoul -= (1.0 - factor_coul) * prefactor;
      }
    } else {
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

  PairContribution single_coul(LAMMPS *lmp, int i, int j, double rsq, double factor_coul) const
  {
    using namespace EwaldConst;
    PairContribution out{0.0, 0.0};
    if (rsq >= cut_coul_global * cut_coul_global) return out;

    const double *q = lmp->atom->q;
    const double r2inv = 1.0 / rsq;
    const double qiqj = q[i] * q[j];
    double forcecoul, phicoul, prefactor = 0.0;

    if (!ncoultablebits || rsq <= tabinnersq) {
      const double r = sqrt(rsq);
      const double grij = g_ewald * r;
      const double expm2 = exp(-grij * grij);
      const double t = 1.0 / (1.0 + EWALD_P * grij);
      const double erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
      prefactor = qqrd2e * qiqj / r;
      forcecoul = prefactor * (erfc + EWALD_F * grij * expm2);
      if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor;
      phicoul = prefactor * erfc;
      if (factor_coul < 1.0) phicoul -= (1.0 - factor_coul) * prefactor;
    } else {
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
      table = etable[itable] + fraction * detable[itable];
      phicoul = qiqj * table;
      if (factor_coul < 1.0) phicoul -= (1.0 - factor_coul) * prefactor;
    }

    out.fpair = forcecoul * r2inv;
    out.energy = phicoul;
    return out;
  }
};

}    // namespace LAMMPS_NS::functor

#endif
