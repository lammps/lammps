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

#ifndef COMPUTE_EROTATE_ASPHERE_H
#define COMPUTE_EROTATE_ASPHERE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeERotateASphere : public Compute {
 public:
  ComputeERotateASphere(class LAMMPS *, int, char **);
  ~ComputeERotateASphere();
  void init();
  double compute_scalar();

 private:
  double pfactor;
  double **inertia;

  void calculate_inertia();
};

}

#endif
