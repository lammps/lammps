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

#ifndef COMPUTE_TEMP_ASPHERE_H
#define COMPUTE_TEMP_ASPHERE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempAsphere : public Compute {
 public:
  ComputeTempAsphere(class LAMMPS *, int, char **);
  ~ComputeTempAsphere();
  void init();
  double compute_scalar();
  void compute_vector();

  void remove_bias(int, double *);
  void remove_bias_all();
  void restore_bias(double *);
  void restore_bias_all();

 private:
  int fix_dof;
  double tfactor;
  double **inertia;

  Compute *tbias;     // ptr to additional bias compute

  void recount();
  void calculate_inertia();
};

}

#endif
