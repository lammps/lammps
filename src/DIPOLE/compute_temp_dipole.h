/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef COMPUTE_TEMP_DIPOLE_H
#define COMPUTE_TEMP_DIPOLE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempDipole : public Compute {
 public:
  ComputeTempDipole(class LAMMPS *, int, char **);
  ~ComputeTempDipole();
  void init();
  double compute_scalar();
  void compute_vector();

 private:
  int fix_dof;
  double tfactor;
  double *inertia;

  void recount();
};

}

#endif
