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

#ifndef LMP_EVALUATOR_MORSE_H
#define LMP_EVALUATOR_MORSE_H

#include "functor_evaluator.h"

namespace LAMMPS_NS::functor {

// Morse potential (reimplements src/pair_morse.cpp).  Unlike Lennard-Jones it
// has no mixing rule, so has_mixing is false and every type pair must be set
// explicitly (the driver errors otherwise).

struct EvaluatorMorse {
  struct Coeff {
    double d0, alpha, r0, cut;
  };
  struct alignas(64) Param {
    double d0, alpha, r0, morse1, offset, cutsq;
  };

  struct Global {};

  static constexpr bool needs_charge = false;
  static constexpr bool has_mixing = false;
  static constexpr bool has_vdw = true;

  static int settings(int narg, char **arg, LAMMPS *lmp, double &cut_global, Global &)
  {
    if (narg < 1) lmp->error->all(FLERR, "Illegal pair_style command");
    cut_global = utils::numeric(FLERR, arg[0], false, lmp);
    return 1;
  }

  static Coeff parse(int narg, char **arg, LAMMPS *lmp, double cut_global, int &nconsumed)
  {
    if (narg < 5)
      lmp->error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
    Coeff c;
    c.d0 = utils::numeric(FLERR, arg[2], false, lmp);
    c.alpha = utils::numeric(FLERR, arg[3], false, lmp);
    c.r0 = utils::numeric(FLERR, arg[4], false, lmp);
    if (narg >= 6) {
      c.cut = utils::numeric(FLERR, arg[5], false, lmp);
      nconsumed = 4;
    } else {
      c.cut = cut_global;
      nconsumed = 3;
    }
    return c;
  }

  static Param derive(const Coeff &c, int offset_flag, const Global &)
  {
    Param p;
    p.d0 = c.d0;
    p.alpha = c.alpha;
    p.r0 = c.r0;
    p.morse1 = 2.0 * c.d0 * c.alpha;
    if (offset_flag) {
      const double alpha_dr = -c.alpha * (c.cut - c.r0);
      p.offset = c.d0 * (exp(2.0 * alpha_dr) - 2.0 * exp(alpha_dr));
    } else
      p.offset = 0.0;
    p.cutsq = c.cut * c.cut;
    return p;
  }

  template <bool EFLAG>
  static PairContribution pair(double rsq, const Param &p, double factor_lj)
  {
    const double r = sqrt(rsq);
    const double dr = r - p.r0;
    const double dexp = exp(-p.alpha * dr);
    PairContribution out;
    out.fpair = factor_lj * p.morse1 * (dexp * dexp - dexp) / r;
    if constexpr (EFLAG)
      out.energy = factor_lj * (p.d0 * (dexp * dexp - 2.0 * dexp) - p.offset);
    else
      out.energy = 0.0;
    return out;
  }
};

}    // namespace LAMMPS_NS::functor

#endif
