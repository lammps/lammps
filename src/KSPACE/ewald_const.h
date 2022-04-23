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

#ifndef LMP_EWALD_CONST_H
#define LMP_EWALD_CONST_H

namespace LAMMPS_NS {
namespace EwaldConst {
  static constexpr double EWALD_F = 1.12837917;
  static constexpr double EWALD_P = 0.3275911;
  static constexpr double A1 = 0.254829592;
  static constexpr double A2 = -0.284496736;
  static constexpr double A3 = 1.421413741;
  static constexpr double A4 = -1.453152027;
  static constexpr double A5 = 1.061405429;
}    // namespace EwaldConst
}    // namespace LAMMPS_NS

#endif
