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

AngleStyle(quartic,AngleQuartic)

#else

#ifndef LMP_ANGLE_QUARTIC_H
#define LMP_ANGLE_QUARTIC_H

#include "angle.h"

namespace LAMMPS_NS {

class AngleQuartic : public Angle {
 public:
  AngleQuartic(class LAMMPS *);
  virtual ~AngleQuartic();
  virtual void compute(int, int);
  void coeff(int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);
  double single(int, int, int, int);

 protected:
  double *k2, *k3, *k4, *theta0;

  void allocate();
};

}

#endif
#endif
