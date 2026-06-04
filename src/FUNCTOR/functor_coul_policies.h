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
//   static constexpr double rsq_epsilon;// tiny value added to rsq in the loop
//                                        // (0.0 for all but the core-shell "/cs"
//                                        // policies, which add 1e-20 to survive
//                                        // r -> 0 for overlapping core/shell)
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
#include <cstddef>
#include <cstdio>

namespace LAMMPS_NS::functor {

// no electrostatics (pure van der Waals styles: lj/cut, morse, ...)

struct CoulNone {
  static constexpr bool has_coul = false;
  static constexpr bool has_table = false;
  static constexpr bool needs_charge = false;
  static constexpr bool needs_kspace = false;
  static constexpr double rsq_epsilon = 0.0;
};

// Shared storage for Coulomb policies that carry a per-type-pair cutoff: the raw
// cutoff (for mixing + restart) and its square (hot loop), plus the global
// cutoff and qqrd2e.  The driver allocates/fills/frees these through the
// lifecycle hooks here; mixing of unset pairs is done by the driver (it has
// Pair::mix_distance).  Concrete policies add the kernel and init_style.

struct CoulCutoffBase {
  static constexpr double rsq_epsilon = 0.0;    // overridden by the "/cs" policies

  double *cut_coul = nullptr;      // per-pair raw cutoff,  flat [n*n]
  double *cut_coulsq = nullptr;    // per-pair squared cutoff, flat [n*n] (hot loop)
  int n = 0;                       // matrix stride (== driver's nparams)
  double cut_coul_global = 0.0;
  double qqrd2e = 1.0;

  void allocate(int n_)
  {
    n = n_;
    cut_coul = new double[(std::size_t) n * n];
    cut_coulsq = new double[(std::size_t) n * n];
  }
  void deallocate()
  {
    delete[] cut_coul;
    delete[] cut_coulsq;
    cut_coul = nullptr;
    cut_coulsq = nullptr;
  }

  // parse the optional per-pair Coulomb cutoff at arg[iarg] (else use the global)
  double parse_cut(int narg, char **arg, LAMMPS *lmp, int iarg) const
  {
    return (narg > iarg) ? utils::numeric(FLERR, arg[iarg], false, lmp) : cut_coul_global;
  }

  // restart: only the global cutoff is a "setting"; per-pair raw cutoffs are
  // written/read by the driver alongside the evaluator coefficients
  void write_restart_settings(FILE *fp) const { fwrite(&cut_coul_global, sizeof(double), 1, fp); }
  void read_restart_settings(FILE *fp, LAMMPS *lmp)
  {
    if (lmp->comm->me == 0)
      utils::sfread(FLERR, &cut_coul_global, sizeof(double), 1, fp, nullptr, lmp->error);
    MPI_Bcast(&cut_coul_global, 1, MPI_DOUBLE, 0, lmp->world);
  }
};

// cutoff Coulomb with a per-type-pair cutoff (bare coul/cut and ".../coul/cut"
// combined styles):
//   F_coul/r = factor_coul * qqrd2e * qi * qj / r^3 ,  E_coul = factor_coul * qqrd2e * qi * qj / r
// (reimplements the Coulomb half of src/pair_coul_cut.cpp / pair_lj_cut_coul_cut.cpp).

struct CoulCut : CoulCutoffBase {
  static constexpr bool has_coul = true;
  static constexpr bool has_table = false;
  static constexpr bool needs_charge = true;
  static constexpr bool needs_kspace = false;
  static constexpr bool per_pair_cutoff = true;    // accepts a per-pair pair_coeff cutoff

  void init_style(Pair *, LAMMPS *lmp) { qqrd2e = lmp->force->qqrd2e; }

  // the trailing per-pair Param is unused here; it carries soft-core (FEP)
  // coefficients that only the "/soft" Coulomb policies read (see functor_coul_soft.h)
  template <bool EFLAG, int /*CTABLE*/, class P>
  PairContribution eval_coul(double rsq, double qi, double qj, double factor_coul,
                             const P & /*p*/) const
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

  template <class P>
  PairContribution single_coul(LAMMPS *lmp, int i, int j, double rsq, double factor_coul,
                               const P & /*p*/) const
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
