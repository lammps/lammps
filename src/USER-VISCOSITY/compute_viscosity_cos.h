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

/* ----------------------------------------------------------------------
   Contributing author: Zheng GONG (ENS Lyon)
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(viscosity/cos,ComputeViscosityCos)

#else

#ifndef LMP_COMPUTE_VISCOSITY_COS_H
#define LMP_COMPUTE_VISCOSITY_COS_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeViscosityCos : public Compute {
 public:
  ComputeViscosityCos(class LAMMPS *, int, char **);
  virtual ~ComputeViscosityCos();
  void init() {}
  void setup();
  virtual double compute_scalar();
  virtual void compute_vector();

  void remove_bias(int, double *);
  void remove_bias_thr(int, double *, double *);
  void remove_bias_all();
  void restore_bias(int, double *);
  void restore_bias_thr(int, double *, double *);
  void restore_bias_all();

 protected:
  double tfactor;
  double V;

  void dof_compute();
  void calc_V();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Temperature compute degrees of freedom < 0

This should not happen if you are calculating the temperature
on a valid set of atoms.

*/
