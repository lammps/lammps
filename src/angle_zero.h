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

#ifdef ANGLE_CLASS
// clang-format off
AngleStyle(zero,AngleZero);
// clang-format on
#else

#ifndef LMP_ANGLE_ZERO_H
#define LMP_ANGLE_ZERO_H

#include "angle.h"

namespace LAMMPS_NS {

class AngleZero : public Angle {
 public:
  AngleZero(class LAMMPS *);
  virtual ~AngleZero();
  virtual void compute(int, int);
  virtual void coeff(int, char **);
  virtual void settings(int, char **);

  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);

  double single(int, int, int, int);

 protected:
  double *theta0;
  int coeffflag;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

UNDOCUMENTED

E: Incorrect args for angle coefficients

Self-explanatory.  Check the input script or data file.

*/
