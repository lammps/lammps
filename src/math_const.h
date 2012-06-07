/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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
  static const double THIRD  = 1.0/3.0;
  static const double MY_PI  = 3.14159265358979323846; // pi
  static const double MY_2PI = 6.28318530717958647692; // 2pi
  static const double MY_3PI = 9.42477796076937971538; // 3pi
  static const double MY_4PI = 12.56637061435917295384; // 4pi
  static const double MY_PI2 = 1.57079632679489661923; // pi/2
  static const double MY_PI4 = 0.78539816339744830962; // pi/4
  static const double MY_PIS = 1.77245385090551602729; // sqrt(pi)
}

}

#endif
