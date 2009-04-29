/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef COMPUTE_TEMP_PARTIAL_H
#define COMPUTE_TEMP_PARTIAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempPartial : public Compute {
 public:
  ComputeTempPartial(class LAMMPS *, int, char **);
  ~ComputeTempPartial();
  void init();
  double compute_scalar();
  void compute_vector();

  int dof_remove(int);
  void remove_bias(int, double *);
  void remove_bias_all();
  void restore_bias(int, double *);
  void restore_bias_all();
  double memory_usage();

 private:
  int xflag,yflag,zflag;
  int fix_dof;
  double tfactor;

  void dof_compute();
};

}

#endif
