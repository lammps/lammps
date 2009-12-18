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

#ifndef COMPUTE_ANGLE_LOCAL_H
#define COMPUTE_ANGLE_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeAngleLocal : public Compute {
 public:
  ComputeAngleLocal(class LAMMPS *, int, char **);
  ~ComputeAngleLocal();
  void init();
  void compute_local();
  double memory_usage();

 private:
  int nvalues,dflag,eflag;
  int *which;
  int ncount;

  int nmax;
  double *theta;
  double *energy;
  double **array;
  double *buf;

  int compute_angles(int);
  void reallocate(int);

  typedef void (ComputeAngleLocal::*FnPtrPack)(int);
  FnPtrPack *pack_choice;              // ptrs to pack functions

  void pack_theta(int);
  void pack_energy(int);
};

}

#endif
