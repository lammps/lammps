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

#ifndef LMP_EVALUATOR_NONE_H
#define LMP_EVALUATOR_NONE_H

#include "functor_evaluator.h"

namespace LAMMPS_NS::functor {

// "no van der Waals" evaluator: contributes nothing.  Combined with a Coulomb
// policy it yields a pure electrostatics style (e.g. coul/cut/functor).  Its
// cutoff is zero (the vdW term is never within range), so the overall pair
// cutoff is set entirely by the Coulomb policy.  has_vdw is false, which tells
// the driver that pair_style/pair_coeff carry no van der Waals cutoff and the
// leading cutoff argument belongs to the Coulomb policy.

struct EvaluatorNone {
  struct Coeff {
    double cut;    // always 0.0; present for driver uniformity
  };
  struct alignas(64) Param {
    double cutsq;    // always 0.0
  };

  static constexpr bool needs_charge = false;
  static constexpr bool has_mixing = true;
  static constexpr bool has_vdw = false;

  static Coeff parse(int /*narg*/, char ** /*arg*/, LAMMPS * /*lmp*/, double /*cut_global*/,
                     int &nconsumed)
  {
    nconsumed = 0;
    return Coeff{0.0};
  }

  static Coeff mix(const Coeff &, const Coeff &, Pair *) { return Coeff{0.0}; }

  static Param derive(const Coeff &, int /*offset_flag*/) { return Param{0.0}; }

  template <bool /*EFLAG*/>
  static PairContribution pair(double /*rsq*/, const Param &, double /*factor_lj*/)
  {
    return PairContribution{0.0, 0.0};
  }
};

}    // namespace LAMMPS_NS::functor

#endif
