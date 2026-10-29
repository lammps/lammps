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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(sw,PairSW);
// clang-format on
#else

#ifndef LMP_PAIR_SW_H
#define LMP_PAIR_SW_H

#include "pair.h"

#include "math_const.h"

#include <cmath>

namespace LAMMPS_NS {

class PairSW : public Pair {
 public:
  PairSW(class LAMMPS *);
  ~PairSW() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void init_style() override;
  void settings(int, char **) override;
  double single(int, int, int, int, double, double, double, double &) override;

  static constexpr int NPARAMS_PER_LINE = 14;

  struct Param {
    double epsilon, sigma;
    double littlea, lambda, gamma, costheta;
    double biga, bigb;
    double powerp, powerq;
    double tol;
    double cut, cutsq;
    double sigma_gamma, lambda_epsilon, lambda_epsilon2;
    double c1, c2, c3, c4, c5, c6;
    int ielement, jelement, kelement;
  };

  struct ShortLeg {
    double del[3], r, rinvsq, gsrainvsq, expgsrainv;
    int j, jtype;
  };

 protected:
  double cutmax;              // max cutoff for all elements
  Param *params;              // parameter set for an I-J-K interaction
  int maxshort;               // size of short neighbor list array
  ShortLeg *shortleg;         // short neighbor list with per-neighbor three-body data
  int skip_threebody_flag;    // whether to run threebody loop
  int sw_mod_flag;            // 1 if the angular term is modified (sw/mod)
  double delta1, delta2;      // range of the modification (sw/mod)
  int params_mapped;          // whether parameters have been read and mapped to elements

  virtual void allocate();
  virtual void read_file(char *);
  virtual void setup_params();
  void twobody(Param *, double, double &, int, double &);

  inline void fill_shortleg(const Param *const param, const double rsq, const double dx,
                            const double dy, const double dz, const int j, const int jtype,
                            ShortLeg &s) const
  {
    s.del[0] = dx;
    s.del[1] = dy;
    s.del[2] = dz;
    s.j = j;
    s.jtype = jtype;
    const double r = sqrt(rsq);
    const double rainv = 1.0 / (r - param->cut);
    const double gsrainv = param->sigma_gamma * rainv;
    s.r = r;
    s.rinvsq = 1.0 / rsq;
    s.gsrainvsq = gsrainv * rainv / r;
    s.expgsrainv = exp(gsrainv);
  }

  template <int MODFLAG>
  inline void threebody_cached(const Param *const paramijk, const ShortLeg &s1, const ShortLeg &s2,
                               double *const fj, double *const fk, const int eflag,
                               double &eng) const
  {
    const double rinv12 = 1.0 / (s1.r * s2.r);
    const double cs =
        (s1.del[0] * s2.del[0] + s1.del[1] * s2.del[1] + s1.del[2] * s2.del[2]) * rinv12;
    double delcs = cs - paramijk->costheta;
    double dfactor = 1.0;
    if (MODFLAG) {
      const double adelcs = fabs(delcs);
      if (adelcs >= delta2) {
        delcs = 0.0;
        dfactor = 0.0;
      } else if (adelcs > delta1) {
        const double arg = MathConst::MY_PI * (adelcs - delta1) / (delta2 - delta1);
        const double factor = 0.5 + 0.5 * cos(arg);
        dfactor = factor - 0.5 * MathConst::MY_PI * adelcs * sin(arg) / (delta2 - delta1);
        delcs *= factor;
      }
    }
    const double delcssq = delcs * delcs;
    const double facexp = s1.expgsrainv * s2.expgsrainv;
    const double facrad = paramijk->lambda_epsilon * facexp * delcssq;
    const double frad1 = facrad * s1.gsrainvsq;
    const double frad2 = facrad * s2.gsrainvsq;
    double facang = paramijk->lambda_epsilon2 * facexp * delcs;
    if (MODFLAG) facang *= dfactor;
    const double facang12 = rinv12 * facang;
    const double csfacang = cs * facang;
    const double csfac1 = s1.rinvsq * csfacang;
    const double csfac2 = s2.rinvsq * csfacang;

    fj[0] = s1.del[0] * (frad1 + csfac1) - s2.del[0] * facang12;
    fj[1] = s1.del[1] * (frad1 + csfac1) - s2.del[1] * facang12;
    fj[2] = s1.del[2] * (frad1 + csfac1) - s2.del[2] * facang12;

    fk[0] = s2.del[0] * (frad2 + csfac2) - s1.del[0] * facang12;
    fk[1] = s2.del[1] * (frad2 + csfac2) - s1.del[1] * facang12;
    fk[2] = s2.del[2] * (frad2 + csfac2) - s1.del[2] * facang12;

    if (eflag) eng = facrad;
  }
};

}    // namespace LAMMPS_NS

#endif
#endif
