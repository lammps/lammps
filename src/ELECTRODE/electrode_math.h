/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meißner (TUHH)
------------------------------------------------------------------------- */

#ifndef LMP_ELECTRODE_MATH_H
#define LMP_ELECTRODE_MATH_H

#include "math_const.h"

namespace LAMMPS_NS {

namespace ElectrodeMath {
  static constexpr double EWALD_P = 0.3275911;
  static constexpr double A1 = 0.254829592;
  static constexpr double A2 = -0.284496736;
  static constexpr double A3 = 1.421413741;
  static constexpr double A4 = -1.453152027;
  static constexpr double A5 = 1.061405429;
  static constexpr double ERFCMAX = 5.8;    // erfc(ERFCMAX) < machine epsilon(double)

  static double safe_erfc(double x)
  {
    if (x > ERFCMAX) return 0.;
    double expm2 = exp(-x * x);
    double t = 1.0 / (1.0 + EWALD_P * x);
    return t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
  }

  static double safe_derfcr(double x, double &erfc)
  {
    if (x > ERFCMAX) {
      erfc = 0.;
      return 0.;
    }
    double x2 = x * x;
    double expm2 = exp(-x2);
    double t = 1.0 / (1.0 + EWALD_P * x);
    erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
    return -erfc - 2 * expm2 * x / MathConst::MY_PIS;
  }
}    // namespace ElectrodeMath

}    // namespace LAMMPS_NS

#endif
