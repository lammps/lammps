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

#ifdef COMPUTE_CLASS

ComputeStyle(temp/sphere,ComputeTempSphere)

#else

#ifndef LMP_COMPUTE_TEMP_SPHERE_H
#define LMP_COMPUTE_TEMP_SPHERE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempSphere : public Compute {
 public:
  ComputeTempSphere(class LAMMPS *, int, char **);
  ~ComputeTempSphere();
  void init();
  void setup();
  double compute_scalar();
  void compute_vector();

  void remove_bias(int, double *);
  void restore_bias(int, double *);

 private:
  int mode;
  double tfactor;
  double *inertia;
  char *id_bias;
  Compute *tbias;     // ptr to additional bias compute

  void dof_compute();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute temp/sphere requires atom style sphere

Self-explanatory.

E: Could not find compute ID for temperature bias

Self-explanatory.

E: Bias compute does not calculate temperature

The specified compute must compute temperature.

E: Bias compute does not calculate a velocity bias

The specified compute must compute a bias for temperature.

E: Bias compute group does not match compute group

The specified compute must operate on the same group as the parent
compute.

*/
