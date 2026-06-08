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

#ifndef LMP_EVALUATOR_NONE_SOFT_H
#define LMP_EVALUATOR_NONE_SOFT_H

#include "functor_evaluator.h"

namespace LAMMPS_NS::functor {

// Zero van der Waals evaluator for the pure soft-core Coulomb styles
// (coul/cut/soft, coul/long/soft).  Like EvaluatorNone it contributes no vdW
// term, but it still owns the per-type-pair free-energy coupling parameter
// lambda (a pair_coeff argument) and derives the two quantities the "/soft"
// Coulomb policies read from Param: lam1 = lambda^nlambda and
// lam2 = alphac*(1-lambda)^2.  nlambda and alphac are style-global (Global).

struct EvaluatorNoneSoft {
  struct Coeff {
    double lambda, cut;    // cut is the (zero) vdW cutoff, for driver uniformity
  };
  struct alignas(64) Param {
    double cutsq;    // always 0.0 (no vdW)
    double lam1;     // lambda^nlambda
    double lam2;     // alphac*(1-lambda)^2
  };
  struct Global {
    double nlambda, alphac;
  };

  static constexpr bool needs_charge = false;
  static constexpr bool has_mixing = true;
  static constexpr bool has_vdw = false;

  // pair_style coul/{cut,long}/soft nlambda alphac cutoff
  // (the trailing Coulomb cutoff is consumed by the policy, not here)
  static int settings(int narg, char **arg, LAMMPS *lmp, double &cut_global, Global &g)
  {
    if (narg < 2) lmp->error->all(FLERR, "Illegal pair_style command");
    g.nlambda = utils::numeric(FLERR, arg[0], false, lmp);
    g.alphac = utils::numeric(FLERR, arg[1], false, lmp);
    cut_global = 0.0;
    return 2;
  }

  // pair_coeff I J lambda
  static Coeff parse(int narg, char **arg, LAMMPS *lmp, double /*cut_global*/, int &nconsumed)
  {
    if (narg < 3)
      lmp->error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
    Coeff c;
    c.lambda = utils::numeric(FLERR, arg[2], false, lmp);
    c.cut = 0.0;
    nconsumed = 1;
    return c;
  }

  static Coeff mix(const Coeff &a, const Coeff &b, Pair *, LAMMPS *lmp)
  {
    if (a.lambda != b.lambda)
      lmp->error->all(FLERR, "Pair coul/.../soft different lambda values in mix");
    return Coeff{a.lambda, 0.0};
  }

  static Param derive(const Coeff &c, int /*offset_flag*/, const Global &g)
  {
    Param p;
    p.cutsq = 0.0;
    p.lam1 = pow(c.lambda, g.nlambda);
    p.lam2 = g.alphac * (1.0 - c.lambda) * (1.0 - c.lambda);
    return p;
  }

  template <bool /*EFLAG*/>
  static PairContribution pair(double /*rsq*/, const Param &, double /*factor_lj*/)
  {
    return PairContribution{0.0, 0.0};
  }
};

}    // namespace LAMMPS_NS::functor

#endif
