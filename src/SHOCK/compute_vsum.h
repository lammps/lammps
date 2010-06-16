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

// Compute the sum of the squared velocities, without usual mass terms
// for the kinetic energy.  This is used for the MSST (Larry Fried).

#ifdef COMPUTE_CLASS

ComputeStyle(vsum,ComputeVsum)

#else

#ifndef COMPUTE_VSUM_H
#define COMPUTE_VSUM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeVsum : public Compute {
 public:
  ComputeVsum(class LAMMPS *, int, char **);
  ~ComputeVsum();
  void init();
  double compute_scalar();

};

}

#endif
#endif
