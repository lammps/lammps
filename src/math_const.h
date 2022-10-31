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

#ifndef LMP_MATH_CONST_H
#define LMP_MATH_CONST_H

namespace LAMMPS_NS {

namespace MathConst {
  static constexpr double THIRD = 1.0 / 3.0;
  static constexpr double TWOTHIRDS = 2.0 / 3.0;
  static constexpr double MY_PI = 3.14159265358979323846;           // pi
  static constexpr double MY_2PI = 6.28318530717958647692;          // 2pi
  static constexpr double MY_3PI = 9.42477796076937971538;          // 3pi
  static constexpr double MY_4PI = 12.56637061435917295384;         // 4pi
  static constexpr double MY_4PI3 = 4.18879020478639098461;         // 4/3pi
  static constexpr double MY_PI2 = 1.57079632679489661923;          // pi/2
  static constexpr double MY_PI4 = 0.78539816339744830962;          // pi/4
  static constexpr double MY_PIS = 1.77245385090551602729;          // sqrt(pi)
  static constexpr double MY_ISPI4 = 1.12837916709551257389;        // 1/sqrt(pi/4)
  static constexpr double MY_SQRT2 = 1.41421356237309504880;        // sqrt(2)
  static constexpr double MY_ISQRT2 = 0.707106781186547524401;      // 1/sqrt(2)
  static constexpr double MY_CUBEROOT2 = 1.25992104989487316476;    // 2^(1/3)
  static constexpr double MY_TWOBYSIXTH = 1.12246204830937298142;   // 2^(1/6)                                                                  //
  static constexpr double DEG2RAD = MY_PI / 180.0;                  // degree to radians
  static constexpr double RAD2DEG = 180.0 / MY_PI;                  // radians to degree
}    // namespace MathConst

}    // namespace LAMMPS_NS

#endif
