/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ------------------------------------------------------------------------- */

/*------------------------------------------------------------------------
  Contributing Authors : Romain Vermorel (LFCR), Laurent Joly (ULyon)
  --------------------------------------------------------------------------*/

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(stress/mop,ComputeStressMop);
// clang-format on
#else

#ifndef LMP_COMPUTE_STRESS_MOP_H
#define LMP_COMPUTE_STRESS_MOP_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeStressMop : public Compute {
 public:
  ComputeStressMop(class LAMMPS *, int, char **);
  ~ComputeStressMop() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_vector() override;

 private:
  void compute_pairs();

  int me, nvalues, dir;
  int *which;

  double *values_local, *values_global;
  double pos, pos1, dt, nktv2p, ftm2v;
  double area;
  class NeighList *list;
};

}    // namespace LAMMPS_NS

#endif
#endif
