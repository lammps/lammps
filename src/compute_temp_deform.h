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

ComputeStyle(temp/deform,ComputeTempDeform)

#else

#ifndef LMP_COMPUTE_TEMP_DEFORM_H
#define LMP_COMPUTE_TEMP_DEFORM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempDeform : public Compute {
 public:
  ComputeTempDeform(class LAMMPS *, int, char **);
  virtual ~ComputeTempDeform();
  void init();
  void setup();
  virtual double compute_scalar();
  virtual void compute_vector();

  void remove_bias(int, double *);
  void remove_bias_thr(int, double *, double *);
  void remove_bias_all();
  void restore_bias(int, double *);
  void restore_bias_thr(int, double *, double *);
  void restore_bias_all();
  double memory_usage();

 protected:
  double tfactor;

  virtual void dof_compute();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

W: Using compute temp/deform with inconsistent fix deform remap option

Fix nvt/sllod assumes deforming atoms have a velocity profile provided
by "remap v" or "remap none" as a fix deform option.

W: Using compute temp/deform with no fix deform defined

This is probably an error, since it makes little sense to use
compute temp/deform in this case.

E: Temperature compute degrees of freedom < 0

This should not happen if you are calculating the temperature
on a valid set of atoms.

*/
