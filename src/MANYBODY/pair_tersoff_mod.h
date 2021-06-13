/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(tersoff/mod,PairTersoffMOD);
// clang-format on
#else

#ifndef LMP_PAIR_TERSOFF_MOD_H
#define LMP_PAIR_TERSOFF_MOD_H

#include "pair_tersoff.h"
#include <cmath>

namespace LAMMPS_NS {

class PairTersoffMOD : public PairTersoff {
 public:
  PairTersoffMOD(class LAMMPS *);
  ~PairTersoffMOD() {}

  static constexpr int NPARAMS_PER_LINE = 20;

 protected:
  virtual void read_file(char *);
  virtual void setup_params();
  double zeta(Param *, double, double, double *, double *);

  double ters_fc(double, Param *);
  double ters_fc_d(double, Param *);
  double ters_bij(double, Param *);
  double ters_bij_d(double, Param *);
  void ters_zetaterm_d(double, double *, double, double, double *, double, double, double *,
                       double *, double *, Param *);

  // inlined functions for efficiency
  // these replace but do not override versions in PairTersoff
  // since overriding virtual inlined functions is best avoided

  inline double ters_gijk_mod(const double costheta, const Param *const param) const
  {
    const double ters_c1 = param->c1;
    const double ters_c2 = param->c2;
    const double ters_c3 = param->c3;
    const double ters_c4 = param->c4;
    const double ters_c5 = param->c5;
    const double tmp_h = (param->h - costheta) * (param->h - costheta);

    return ters_c1 +
        (ters_c2 * tmp_h / (ters_c3 + tmp_h)) * (1.0 + ters_c4 * exp(-ters_c5 * tmp_h));
  }

  inline double ters_gijk_d_mod(const double costheta, const Param *const param) const
  {
    const double ters_c2 = param->c2;
    const double ters_c3 = param->c3;
    const double ters_c4 = param->c4;
    const double ters_c5 = param->c5;
    const double tmp_h = (param->h - costheta) * (param->h - costheta);
    const double g1 = (param->h - costheta) / (ters_c3 + tmp_h);
    const double g2 = exp(-ters_c5 * tmp_h);

    return -2.0 * ters_c2 * g1 *
        ((1 + ters_c4 * g2) * (1 + g1 * (costheta - param->h)) - tmp_h * ters_c4 * ters_c5 * g2);
  }
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot open Tersoff potential file %s

The specified potential file cannot be opened.  Check that the path
and name are correct.

E: Incorrect format in Tersoff potential file

Incorrect number of words per line in the potential file.

E: Illegal Tersoff parameter

One or more of the coefficients defined in the potential file is
invalid.

E: Potential file has duplicate entry

The potential file has more than one entry for the same element.

E: Potential file is missing an entry

The potential file does not have a needed entry.

*/
