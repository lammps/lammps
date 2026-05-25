/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meissner (TUHH)
------------------------------------------------------------------------- */

#ifndef LMP_ELECTRODE_MATH_H
#define LMP_ELECTRODE_MATH_H

#include "math_const.h"
#include "math_special.h"

#include <cmath>

namespace LAMMPS_NS {

namespace ElectrodeMath {
  static constexpr double ERFCMAX = 5.8;    // erfc(ERFCMAX) < machine epsilon(double)

  inline double safe_erfc(double x)
  {
    if (x > ERFCMAX) return 0.0;
    return MathSpecial::my_erfcx(x) * MathSpecial::expmsq(x);
  }

  inline double safe_derfcr(double x, double &erfc)
  {
    if (x > ERFCMAX) {
      erfc = 0.0;
      return 0.0;
    }
    double expm2 = MathSpecial::expmsq(x);
    erfc = MathSpecial::my_erfcx(x) * expm2;
    return -erfc - 2.0 * expm2 * x / MathConst::MY_PIS;
  }
}    // namespace ElectrodeMath

}    // namespace LAMMPS_NS

#endif
