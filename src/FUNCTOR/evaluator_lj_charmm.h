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

#ifndef LMP_EVALUATOR_LJ_CHARMM_H
#define LMP_EVALUATOR_LJ_CHARMM_H

#include "functor_evaluator.h"

namespace LAMMPS_NS::functor {

// CHARMM Lennard-Jones with a force/energy switching function between an inner
// and outer cutoff (reimplements the LJ half of
// src/KSPACE/pair_lj_charmm_coul_long.cpp).  The inner/outer cutoffs are global
// (a single pair of values for all type pairs) and arrive through the Global
// struct.  The 1-4 parameters eps14/sigma14 are accepted on the pair_coeff line
// (and in data files) but, as in the reference, are not used by the pairwise
// kernel of this style, so they are parsed and ignored.

struct EvaluatorLJCharmm {
  struct Coeff {
    double epsilon, sigma, eps14, sigma14, cut;
  };
  struct alignas(64) Param {
    double lj1, lj2, lj3, lj4;
    double cut_lj_innersq, cutsq, denom_lj_inv;
  };
  struct Global {
    double cut_lj_innersq, cut_ljsq, denom_lj_inv;
  };

  static constexpr bool needs_charge = false;
  static constexpr bool has_mixing = true;
  static constexpr bool has_vdw = true;

  // pair_style lj/charmm/coul/long cut_lj_inner cut_lj [cut_coul]
  static int settings(int narg, char **arg, LAMMPS *lmp, double &cut_global, Global &g)
  {
    if (narg < 2) lmp->error->all(FLERR, "Illegal pair_style command");
    const double cut_lj_inner = utils::numeric(FLERR, arg[0], false, lmp);
    cut_global = utils::numeric(FLERR, arg[1], false, lmp);    // outer LJ cutoff
    g.cut_lj_innersq = cut_lj_inner * cut_lj_inner;
    g.cut_ljsq = cut_global * cut_global;
    const double denom = (g.cut_ljsq - g.cut_lj_innersq);
    g.denom_lj_inv = 1.0 / (denom * denom * denom);
    return 2;
  }

  // pair_coeff I J epsilon sigma [eps14 sigma14]
  static Coeff parse(int narg, char **arg, LAMMPS *lmp, double cut_global, int &nconsumed)
  {
    if (narg != 4 && narg != 6)
      lmp->error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
    Coeff c;
    c.epsilon = utils::numeric(FLERR, arg[2], false, lmp);
    c.sigma = utils::numeric(FLERR, arg[3], false, lmp);
    if (narg == 6) {
      c.eps14 = utils::numeric(FLERR, arg[4], false, lmp);
      c.sigma14 = utils::numeric(FLERR, arg[5], false, lmp);
      nconsumed = 4;
    } else {
      c.eps14 = c.epsilon;
      c.sigma14 = c.sigma;
      nconsumed = 2;
    }
    c.cut = cut_global;
    return c;
  }

  static Coeff mix(const Coeff &a, const Coeff &b, Pair *p)
  {
    Coeff c;
    c.epsilon = p->mix_energy(a.epsilon, b.epsilon, a.sigma, b.sigma);
    c.sigma = p->mix_distance(a.sigma, b.sigma);
    c.eps14 = p->mix_energy(a.eps14, b.eps14, a.sigma14, b.sigma14);
    c.sigma14 = p->mix_distance(a.sigma14, b.sigma14);
    c.cut = p->mix_distance(a.cut, b.cut);
    return c;
  }

  static Param derive(const Coeff &c, int /*offset_flag*/, const Global &g)
  {
    Param p;
    p.lj1 = 48.0 * c.epsilon * pow(c.sigma, 12.0);
    p.lj2 = 24.0 * c.epsilon * pow(c.sigma, 6.0);
    p.lj3 = 4.0 * c.epsilon * pow(c.sigma, 12.0);
    p.lj4 = 4.0 * c.epsilon * pow(c.sigma, 6.0);
    p.cut_lj_innersq = g.cut_lj_innersq;
    p.cutsq = g.cut_ljsq;
    p.denom_lj_inv = g.denom_lj_inv;
    return p;
  }

  template <bool EFLAG>
  static PairContribution pair(double rsq, const Param &p, double factor_lj)
  {
    const double r2inv = 1.0 / rsq;
    const double r6inv = r2inv * r2inv * r2inv;
    double forcelj = r6inv * (p.lj1 * r6inv - p.lj2);
    double evdwl = 0.0;
    if constexpr (EFLAG) evdwl = r6inv * (p.lj3 * r6inv - p.lj4);

    if (rsq > p.cut_lj_innersq) {
      const double d = p.cutsq - rsq;
      const double switch1 =
          d * d * (p.cutsq + 2.0 * rsq - 3.0 * p.cut_lj_innersq) * p.denom_lj_inv;
      const double philj = r6inv * (p.lj3 * r6inv - p.lj4);
      const double switch2 = 12.0 * rsq * d * (rsq - p.cut_lj_innersq) * p.denom_lj_inv;
      forcelj = forcelj * switch1 + philj * switch2;
      if constexpr (EFLAG) evdwl *= switch1;
    }

    PairContribution out;
    out.fpair = factor_lj * forcelj * r2inv;
    out.energy = factor_lj * evdwl;
    return out;
  }
};

}    // namespace LAMMPS_NS::functor

#endif
