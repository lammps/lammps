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

// Common definition of taper function and its derivative for interlayer potentials

#ifndef LMP_INTERLAYER_TAPER_H
#define LMP_INTERLAYER_TAPER_H

namespace LAMMPS_NS {
namespace InterLayer {

  static constexpr double Tap_coeff[8] = {1.0, 0.0, 0.0, 0.0, -35.0, 84.0, -70.0, 20.0};

  /* ----Calculate the long-range cutoff term */
  static inline double calc_Tap(double r_ij, double Rcut)
  {
    double Tap, r;

    r = r_ij / Rcut;
    if (r >= 1.0) {
      Tap = 0.0;
    } else {
      Tap = Tap_coeff[7] * r + Tap_coeff[6];
      Tap = Tap * r + Tap_coeff[5];
      Tap = Tap * r + Tap_coeff[4];
      Tap = Tap * r + Tap_coeff[3];
      Tap = Tap * r + Tap_coeff[2];
      Tap = Tap * r + Tap_coeff[1];
      Tap = Tap * r + Tap_coeff[0];
    }

    return (Tap);
  }

  /* ----Calculate the derivatives of long-range cutoff term */
  static inline double calc_dTap(double r_ij, double Rcut)
  {
    double dTap, r;

    r = r_ij / Rcut;
    if (r >= 1.0) {
      dTap = 0.0;
    } else {
      dTap = 7.0 * Tap_coeff[7] * r + 6.0 * Tap_coeff[6];
      dTap = dTap * r + 5.0 * Tap_coeff[5];
      dTap = dTap * r + 4.0 * Tap_coeff[4];
      dTap = dTap * r + 3.0 * Tap_coeff[3];
      dTap = dTap * r + 2.0 * Tap_coeff[2];
      dTap = dTap * r + Tap_coeff[1];
      dTap = dTap / Rcut;
    }

    return (dTap);
  }
}    // namespace InterLayer
}    // namespace LAMMPS_NS
#endif
