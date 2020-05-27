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

AngleStyle(cosine/buck6d, AngleCosineBuck6d)

#else

#ifndef LMP_ANGLE_COSINE_BUCK6D_H
#define LMP_ANGLE_COSINE_BUCK6D_H

#include "angle.h"

namespace LAMMPS_NS {

class AngleCosineBuck6d : public Angle {
 public:
  AngleCosineBuck6d(class LAMMPS *);
  virtual ~AngleCosineBuck6d();
  virtual void compute(int, int);
  void coeff(int, char **);
  void init_style();
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);
  double single(int, int, int, int);

 protected:
  double *k,*th0;
  double *eps,*d0;
  double **buck6d1,**buck6d2,**buck6d3,**buck6d4,**cut_ljsq;
  double **c0,**c1,**c2,**c3,**c4,**c5,**rsmooth_sq,**offset;
  int *multiplicity;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for angle coefficients

Self-explanatory.  Check the input script or data file.

*/
