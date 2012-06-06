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

#ifdef COMPUTE_CLASS

ComputeStyle(reduce/region,ComputeReduceRegion)

#else

#ifndef LMP_COMPUTE_REDUCE_REGION_H
#define LMP_COMPUTE_REDUCE_REGION_H

#include "compute_reduce.h"

namespace LAMMPS_NS {

class ComputeReduceRegion : public ComputeReduce {
 public:
  ComputeReduceRegion(class LAMMPS *, int, char **);
  ~ComputeReduceRegion() {}

 private:
  double compute_one(int, int);
  bigint count(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Fix used in compute reduce not computed at compatible time

Fixes generate their values on specific timesteps.  Compute reduce is
requesting a value on a non-allowed timestep.

*/
