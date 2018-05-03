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

AngleStyle(class2/p6,AngleClass2P6)

#else

#ifndef LMP_ANGLE_CLASS2_P6_H
#define LMP_ANGLE_CLASS2_P6_H

#include <cstdio>
#include "angle.h"

namespace LAMMPS_NS {

class AngleClass2P6 : public Angle {
 public:
  AngleClass2P6(class LAMMPS *);
  virtual ~AngleClass2P6();
  virtual void compute(int, int);
  void coeff(int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);
  double single(int, int, int, int);

 protected:
  double *theta0,*k2,*k3,*k4,*k5,*k6;
  double *bb_k,*bb_r1,*bb_r2;
  double *ba_k1,*ba_k2,*ba_r1,*ba_r2;
  int *setflag_a,*setflag_bb,*setflag_ba;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for angle coefficients

Self-explanatory.  Check the input script or data file.

*/
