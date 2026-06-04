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

#ifndef LMP_FUNCTOR_COUL_POLICIES_H
#define LMP_FUNCTOR_COUL_POLICIES_H

// Compile-time Coulomb policies for the FUNCTOR pair driver (PairFunctor).
//
// A policy is selected as the second template argument of PairFunctor and
// controls whether and how electrostatics are added in the inner loop.  The
// driver inspects the constexpr flags below with "if constexpr", so the entire
// Coulomb code path is removed from the compiled kernel when has_coul is false.
//
// Contract (only the parts used by a given stage need to exist; unused methods
// live in discarded "if constexpr" branches and are never instantiated):
//   static constexpr bool has_coul;     // electrostatics contribute at all
//   static constexpr bool has_table;    // uses bitmapped Ewald tables (CTABLE)
//   static constexpr bool needs_charge; // driver must fetch atom->q
// A has_coul policy additionally provides:
//   void init_style(LAMMPS *);          // cache qqrd2e (and g_ewald for long)
//   template <bool EFLAG, int CTABLE>
//   PairContribution eval_coul(double rsq, double qi, double qj, double factor_coul) const;
//   PairContribution single_coul(LAMMPS *, int i, int j, double rsq, double factor_coul) const;
// The Coulomb cutoff itself (coul.cut_coulsq) is owned by the policy; the driver
// reads it for the outer-loop guard and to merge into the per-pair cutoff.

#include "functor_evaluator.h"    // PairContribution

#include "atom.h"
#include "comm.h"
#include "force.h"

#include <cmath>
#include <cstdio>

namespace LAMMPS_NS::functor {

// no electrostatics (pure van der Waals styles: lj/cut, morse, ...)

struct CoulNone {
  static constexpr bool has_coul = false;
  static constexpr bool has_table = false;
  static constexpr bool needs_charge = false;
};

// cutoff Coulomb with a single global cutoff (".../coul/cut" combined styles).
//   F_coul/r = factor_coul * qqrd2e * qi * qj / r^3 ,  E_coul = factor_coul * qqrd2e * qi * qj / r
// (reimplements the Coulomb half of src/pair_lj_cut_coul_cut.cpp).  A per-pair
// Coulomb cutoff and the bare coul/cut style are a planned extension.

struct CoulCut {
  static constexpr bool has_coul = true;
  static constexpr bool has_table = false;
  static constexpr bool needs_charge = true;

  double cut_coulsq = 0.0;    // global Coulomb cutoff squared
  double qqrd2e = 1.0;

  void init_style(LAMMPS *lmp) { qqrd2e = lmp->force->qqrd2e; }

  // persist the global Coulomb cutoff across restarts (qqrd2e is re-derived in
  // init_style, so it need not be stored)
  void write_restart(FILE *fp) const { fwrite(&cut_coulsq, sizeof(double), 1, fp); }
  void read_restart(FILE *fp, LAMMPS *lmp)
  {
    if (lmp->comm->me == 0)
      utils::sfread(FLERR, &cut_coulsq, sizeof(double), 1, fp, nullptr, lmp->error);
    MPI_Bcast(&cut_coulsq, 1, MPI_DOUBLE, 0, lmp->world);
  }

  template <bool EFLAG, int /*CTABLE*/>
  PairContribution eval_coul(double rsq, double qi, double qj, double factor_coul) const
  {
    const double r2inv = 1.0 / rsq;
    const double rinv = sqrt(r2inv);
    const double forcecoul = qqrd2e * qi * qj * rinv;
    PairContribution out;
    out.fpair = factor_coul * forcecoul * r2inv;
    if constexpr (EFLAG)
      out.energy = factor_coul * forcecoul;
    else
      out.energy = 0.0;
    return out;
  }

  PairContribution single_coul(LAMMPS *lmp, int i, int j, double rsq, double factor_coul) const
  {
    const double r2inv = 1.0 / rsq;
    const double rinv = sqrt(r2inv);
    const double *q = lmp->atom->q;
    const double forcecoul = qqrd2e * q[i] * q[j] * rinv;
    PairContribution out;
    out.fpair = factor_coul * forcecoul * r2inv;
    out.energy = factor_coul * forcecoul;
    return out;
  }
};

}    // namespace LAMMPS_NS::functor

#endif
