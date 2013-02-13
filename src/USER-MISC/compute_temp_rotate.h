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

ComputeStyle(temp/rotate,ComputeTempRotate)

#else

#ifndef LMP_COMPUTE_TEMP_ROTATE_H
#define LMP_COMPUTE_TEMP_ROTATE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempRotate : public Compute {
 public:
  ComputeTempRotate(class LAMMPS *, int, char **);
  ~ComputeTempRotate();
  void init();
  void setup();
  double compute_scalar();
  void compute_vector();

  void remove_bias(int, double *);
  void remove_bias_thr(int, double *, double *);
  void remove_bias_all();
  void restore_bias(int, double *);
  void restore_bias_all();
  void restore_bias_thr(int, double *, double *);

  double memory_usage();

 private:
  int fix_dof;
  double tfactor,masstotal;
  double **vbiasall;  // stored velocity bias for all atoms
  int maxbias;        // size of vbiasall array

  void dof_compute(); //without virtual

};

}

#endif
#endif
