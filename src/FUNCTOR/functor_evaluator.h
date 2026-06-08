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

#ifndef LMP_FUNCTOR_EVALUATOR_H
#define LMP_FUNCTOR_EVALUATOR_H

// Shared infrastructure for FUNCTOR "evaluator" functors.  Each potential lives
// in its own header (e.g. evaluator_lj_cut.h) that includes this file, so an
// optional package can add a pair style by shipping a single evaluator header
// plus a registration entry, without touching any shared file.
//
// An evaluator is a struct with no virtual functions; every method is static
// and inline so it inlines into the templated driver kernel.  Required members
// (see doc/src/Developer_write_pair_functor.rst):
//
//   struct Coeff;                           // raw pair_coeff input; needs "double cut;"
//   struct Param;                           // derived hot-loop data; needs "double cutsq;"
//   static constexpr bool needs_charge;     // does the vdW part itself read q?
//   static constexpr bool has_mixing;       // may off-diagonal pairs be mixed?
//   static constexpr const char *name();
//   static Coeff parse(int narg, char **arg, LAMMPS *lmp, double cut_global);
//   static Coeff mix(const Coeff &, const Coeff &, Pair *);   // iff has_mixing
//   static Param derive(const Coeff &, int offset_flag);
//   template <bool EFLAG>
//   static PairContribution pair(double rsq, const Param &, double factor_lj);
//
// The EFLAG template parameter lets the (relatively expensive) energy term be
// dropped from the force-only code path at compile time.

#include "error.h"
#include "lammps.h"
#include "pair.h"
#include "utils.h"

#include <cmath>

namespace LAMMPS_NS::functor {

// Result of one pairwise kernel: the scalar force divided by r, and the
// matching pair energy (zero on the force-only path).  Returned by both the
// evaluator's pair() and a Coulomb policy's eval_coul(), and consumed via
// structured bindings in the driver.

struct PairContribution {
  double fpair;
  double energy;
};

}    // namespace LAMMPS_NS::functor

#endif
