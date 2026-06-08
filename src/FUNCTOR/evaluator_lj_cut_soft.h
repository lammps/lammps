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

#ifndef LMP_EVALUATOR_LJ_CUT_SOFT_H
#define LMP_EVALUATOR_LJ_CUT_SOFT_H

#include "functor_evaluator.h"

namespace LAMMPS_NS::functor {

// Soft-core (FEP) Lennard-Jones (reimplements the van der Waals half of
// src/FEP/pair_lj_cut_soft.cpp).  The free-energy coupling parameter lambda is
// per type pair; nlambda and alphalj are style-global (Global), as is alphac (the
// soft-core parameter of the companion Coulomb policy, which is carried here so a
// single pair_style line sets it for both halves).  The derived per-pair Param
// also exposes lam1 (= lambda^nlambda, the shared prefactor) and lam2
// (= alphac*(1-lambda)^2) so the "/soft" Coulomb policies can read them.
//
//   denlj = alphalj*(1-lambda)^2 + (r/sigma)^6
//   F/r   = lambda^nlambda * epsilon * (48*r4sig6/denlj^3 - 24*r4sig6/denlj^2)
//   E     = lambda^nlambda * 4*epsilon * (1/denlj^2 - 1/denlj) - offset
// with r4sig6 = r^4/sigma^6, so that r^6/sigma^6 = rsq*r4sig6.

struct EvaluatorLJCutSoft {
  struct Coeff {
    double epsilon, sigma, lambda, cut;
  };
  struct alignas(64) Param {
    double lam1;         // lambda^nlambda  (prefactor for both LJ and Coulomb)
    double sig6;         // sigma^6
    double denlj_off;    // alphalj*(1-lambda)^2
    double epsilon;
    double offset;
    double cutsq;
    double lam2;         // alphac*(1-lambda)^2  (soft-core Coulomb offset)
  };
  struct Global {
    double nlambda, alphalj, alphac;
  };

  static constexpr bool needs_charge = false;
  static constexpr bool has_mixing = true;
  static constexpr bool has_vdw = true;

  // pair_style lj/cut/coul/{cut,long}/soft nlambda alphalj alphac cut_lj [cut_coul]
  static int settings(int narg, char **arg, LAMMPS *lmp, double &cut_global, Global &g)
  {
    if (narg < 4) lmp->error->all(FLERR, "Illegal pair_style command");
    g.nlambda = utils::numeric(FLERR, arg[0], false, lmp);
    g.alphalj = utils::numeric(FLERR, arg[1], false, lmp);
    g.alphac = utils::numeric(FLERR, arg[2], false, lmp);
    cut_global = utils::numeric(FLERR, arg[3], false, lmp);
    return 4;
  }

  // pair_coeff I J epsilon sigma lambda [cut_lj]
  static Coeff parse(int narg, char **arg, LAMMPS *lmp, double cut_global, int &nconsumed)
  {
    if (narg < 5)
      lmp->error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
    Coeff c;
    c.epsilon = utils::numeric(FLERR, arg[2], false, lmp);
    c.sigma = utils::numeric(FLERR, arg[3], false, lmp);
    c.lambda = utils::numeric(FLERR, arg[4], false, lmp);
    if (c.sigma <= 0.0)
      lmp->error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
    if (narg >= 6) {
      c.cut = utils::numeric(FLERR, arg[5], false, lmp);
      nconsumed = 4;
    } else {
      c.cut = cut_global;
      nconsumed = 3;
    }
    return c;
  }

  static Coeff mix(const Coeff &a, const Coeff &b, Pair *p, LAMMPS *lmp)
  {
    if (a.lambda != b.lambda)
      lmp->error->all(FLERR, "Pair lj/cut/.../soft different lambda values in mix");
    Coeff c;
    c.epsilon = p->mix_energy(a.epsilon, b.epsilon, a.sigma, b.sigma);
    c.sigma = p->mix_distance(a.sigma, b.sigma);
    c.lambda = a.lambda;
    c.cut = p->mix_distance(a.cut, b.cut);
    return c;
  }

  static Param derive(const Coeff &c, int offset_flag, const Global &g)
  {
    Param p;
    p.lam1 = pow(c.lambda, g.nlambda);
    p.sig6 = pow(c.sigma, 6.0);
    p.denlj_off = g.alphalj * (1.0 - c.lambda) * (1.0 - c.lambda);
    p.lam2 = g.alphac * (1.0 - c.lambda) * (1.0 - c.lambda);
    p.epsilon = c.epsilon;
    p.cutsq = c.cut * c.cut;
    if (offset_flag && (c.cut > 0.0)) {
      const double denlj = p.denlj_off + pow(c.cut / c.sigma, 6.0);
      p.offset = p.lam1 * 4.0 * c.epsilon * (1.0 / (denlj * denlj) - 1.0 / denlj);
    } else
      p.offset = 0.0;
    return p;
  }

  template <bool EFLAG>
  static PairContribution pair(double rsq, const Param &p, double factor_lj)
  {
    const double r4sig6 = rsq * rsq / p.sig6;
    const double denlj = p.denlj_off + rsq * r4sig6;
    const double forcelj = p.lam1 * p.epsilon *
        (48.0 * r4sig6 / (denlj * denlj * denlj) - 24.0 * r4sig6 / (denlj * denlj));
    PairContribution out;
    out.fpair = factor_lj * forcelj;
    if constexpr (EFLAG)
      out.energy = factor_lj * (p.lam1 * 4.0 * p.epsilon * (1.0 / (denlj * denlj) - 1.0 / denlj) -
                                p.offset);
    else
      out.energy = 0.0;
    return out;
  }
};

}    // namespace LAMMPS_NS::functor

#endif
