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

#ifndef LMP_PAIR_ZBL_CONST_H
#define LMP_PAIR_ZBL_CONST_H

namespace LAMMPS_NS {
namespace PairZBLConstants {

  // ZBL constants

  static constexpr double pzbl = 0.23;
  static constexpr double a0 = 0.46850;
  static constexpr double c1 = 0.02817;
  static constexpr double c2 = 0.28022;
  static constexpr double c3 = 0.50986;
  static constexpr double c4 = 0.18175;
  static constexpr double d1 = 0.20162;
  static constexpr double d2 = 0.40290;
  static constexpr double d3 = 0.94229;
  static constexpr double d4 = 3.19980;
}    // namespace PairZBLConstants
}    // namespace LAMMPS_NS
#endif
