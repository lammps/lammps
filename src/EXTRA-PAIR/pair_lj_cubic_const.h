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

#ifndef LMP_PAIR_LJ_CUBIC_CONST_H
#define LMP_PAIR_LJ_CUBIC_CONST_H

namespace LAMMPS_NS {
namespace PairLJCubicConstants {

  // LJ quantities scaled by epsilon and rmin = sigma*2^1/6

  static constexpr double RT6TWO = 1.1224620483093730;    // 2^1/6
  static constexpr double SS = 1.1086834179687215;        // inflection point (13/7)^1/6
  static constexpr double PHIS = -0.7869822485207097;     // energy at s
  static constexpr double DPHIDS = 2.6899008972047196;    // gradient at s
  static constexpr double A3 = 27.9335700460986445;       // cubic coefficient
  static constexpr double SM = 1.5475372709146737;        // cubic cutoff = s*67/48}
}    // namespace PairLJCubicConstants
}    // namespace LAMMPS_NS
#endif

// python script to compute the constants
//
// sixth = 1.0/6.0
// rmin = pow(2.0,sixth)
// rs   = pow(26.0/7.0,sixth)
// ss   = rs/rmin
// pow6 = pow(1.0/rs,6.0)
// phis = 4.0*pow6*(pow6-1.0)
// dpds = -24.0*pow6*(2.0*pow6-1.0)/rs*rmin
// a3   = 8.0*pow(dpds,3)/(9.0*phis*phis)
// sm   = 67.0/48.0*ss
// print("static constexpr double RT6TWO = %19.16f; // 2^1/6" % rmin)
// print("static constexpr double SS     = %19.16f; // inflection point (13/7)^1/6" % ss)
// print("static constexpr double PHIS   = %19.16f; // energy at s" % phis)
// print("static constexpr double DPHIDS = %19.16f; // gradient at s" % dpds)
// print("static constexpr double A3     = %19.16f; // cubic coefficient" % a3)
// print("static constexpr double SM     = %19.16f; // cubic cutoff = s*67/48" % sm)
