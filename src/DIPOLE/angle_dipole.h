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
AngleStyle(dipole,AngleDipole);
// clang-format on
#else

#ifndef LMP_ANGLE_DIPOLE_H
#define LMP_ANGLE_DIPOLE_H

#include "angle.h"

namespace LAMMPS_NS {

class AngleDipole : public Angle {
 public:
  AngleDipole(class LAMMPS *);
  virtual ~AngleDipole();
  virtual void compute(int, int);
  virtual void init_style();
  virtual void coeff(int, char **);
  virtual double equilibrium_angle(int);
  virtual void write_restart(FILE *);
  virtual void read_restart(FILE *);
  virtual void write_data(FILE *);
  virtual double single(int, int, int, int);

 protected:
  double *k, *gamma0;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for angle coefficients

Self-explanatory.  Check the input script or data file.

*/
