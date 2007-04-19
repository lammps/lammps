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

#ifndef COMPUTE_VARIABLE_H
#define COMPUTE_VARIABLE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeVariable : public Compute {
 public:
  ComputeVariable(class LAMMPS *, int, char **);
  ~ComputeVariable();
  void init();
  double compute_scalar();

 private:
  char *varname;
};

}

#endif
