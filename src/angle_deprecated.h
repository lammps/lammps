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

#ifdef ANGLE_CLASS

AngleStyle(DEPRECATED,AngleDeprecated)

#else

#ifndef LMP_ANGLE_DEPRECATED_H
#define LMP_ANGLE_DEPRECATED_H

#include "angle.h"

namespace LAMMPS_NS {

class AngleDeprecated : public Angle {
 public:
  AngleDeprecated(class LAMMPS *lmp) : Angle(lmp) {}
  virtual ~AngleDeprecated() {}

  virtual void compute(int, int) {}
  virtual void settings(int, char **);
  virtual void coeff(int, char **) {}
  virtual double equilibrium_angle(int) { return 0.0; }
  virtual void write_restart(FILE *) {}
  virtual void read_restart(FILE *) {}
  virtual double single(int, int, int, int) { return 0.0; }
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
