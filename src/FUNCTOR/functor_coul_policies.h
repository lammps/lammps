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
// Stage 2/3 will add the cutoff- and long-range policies with their own
// settings()/init_style() members and an
//   template <bool EFLAG, int CTABLE>
//   PairContribution eval_coul(double rsq, double qi, double qj, double factor_coul) const;
// kernel (see functor_evaluator.h for PairContribution).

namespace LAMMPS_NS::functor {

// no electrostatics (pure van der Waals styles: lj/cut, morse, ...)

struct CoulNone {
  static constexpr bool has_coul = false;
  static constexpr bool has_table = false;
  static constexpr bool needs_charge = false;
};

}    // namespace LAMMPS_NS::functor

#endif
