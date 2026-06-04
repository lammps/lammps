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

#ifndef LMP_FUNCTOR_EVALUATORS_H
#define LMP_FUNCTOR_EVALUATORS_H

// Per-potential "evaluator" functors for the FUNCTOR pair driver.
//
// An evaluator is a small struct (no virtual functions, all methods static and
// inline) describing the functional form of one pairwise-additive potential.
// It is the only thing an author writes to add a new pair style; the shared
// PairFunctor<EVAL,COUL> driver supplies the neighbor loop, parameter storage,
// mixing, restart I/O, and the OPT/OPENMP-style optimizations.
//
// Required members (see doc/src/Developer_write_pair_functor.rst):
//   struct Coeff;                       // raw pair_coeff input (must contain a
//                                       //   "double cut;" member)
//   struct Param;                       // derived hot-loop data (must contain a
//                                       //   "double cutsq;" member)
//   static constexpr bool needs_charge; // does the vdW part itself read q?
//   static constexpr bool has_mixing;   // can off-diagonal pairs be mixed?
//   static const char *name();
//   static Coeff parse(int narg, char **arg, LAMMPS *lmp, double cut_global);
//   static Coeff mix(const Coeff &, const Coeff &, Pair *);   // if has_mixing
//   static Param derive(const Coeff &, int offset_flag);
//   static void  pair(double rsq, const Param &, double factor_lj,
//                     double &fpair, double &evdwl, int eflag);

#include "error.h"
#include "lammps.h"
#include "pair.h"
#include "utils.h"

#include <cmath>

namespace LAMMPS_NS {
namespace functor {

  // ------------------------------------------------------------------
  // Lennard-Jones, cutoff (reimplements src/pair_lj_cut.cpp)
  // ------------------------------------------------------------------

  struct EvaluatorLJCut {
    struct Coeff {
      double epsilon, sigma, cut;
    };
    struct Param {
      double lj1, lj2, lj3, lj4, offset, cutsq;
    };

    static constexpr bool needs_charge = false;
    static constexpr bool has_mixing = true;

    static const char *name() { return "lj/cut/functor"; }

    static Coeff parse(int narg, char **arg, LAMMPS *lmp, double cut_global)
    {
      if (narg < 4 || narg > 5)
        lmp->error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
      Coeff c;
      c.epsilon = utils::numeric(FLERR, arg[2], false, lmp);
      c.sigma = utils::numeric(FLERR, arg[3], false, lmp);
      c.cut = (narg == 5) ? utils::numeric(FLERR, arg[4], false, lmp) : cut_global;
      return c;
    }

    static Coeff mix(const Coeff &a, const Coeff &b, Pair *p)
    {
      Coeff c;
      c.epsilon = p->mix_energy(a.epsilon, b.epsilon, a.sigma, b.sigma);
      c.sigma = p->mix_distance(a.sigma, b.sigma);
      c.cut = p->mix_distance(a.cut, b.cut);
      return c;
    }

    static Param derive(const Coeff &c, int offset_flag)
    {
      Param p;
      p.lj1 = 48.0 * c.epsilon * pow(c.sigma, 12.0);
      p.lj2 = 24.0 * c.epsilon * pow(c.sigma, 6.0);
      p.lj3 = 4.0 * c.epsilon * pow(c.sigma, 12.0);
      p.lj4 = 4.0 * c.epsilon * pow(c.sigma, 6.0);
      if (offset_flag && (c.cut > 0.0)) {
        double ratio = c.sigma / c.cut;
        p.offset = 4.0 * c.epsilon * (pow(ratio, 12.0) - pow(ratio, 6.0));
      } else
        p.offset = 0.0;
      p.cutsq = c.cut * c.cut;
      return p;
    }

    static inline void pair(double rsq, const Param &p, double factor_lj, double &fpair,
                            double &evdwl, int eflag)
    {
      const double r2inv = 1.0 / rsq;
      const double r6inv = r2inv * r2inv * r2inv;
      const double forcelj = r6inv * (p.lj1 * r6inv - p.lj2);
      fpair = factor_lj * forcelj * r2inv;
      if (eflag) evdwl = factor_lj * (r6inv * (p.lj3 * r6inv - p.lj4) - p.offset);
    }
  };

}    // namespace functor
}    // namespace LAMMPS_NS

#endif
