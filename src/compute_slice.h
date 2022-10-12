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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(slice,ComputeSlice);
// clang-format on
#else

#ifndef LMP_COMPUTE_SLICE_H
#define LMP_COMPUTE_SLICE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSlice : public Compute {
 public:
  ComputeSlice(class LAMMPS *, int, char **);
  ~ComputeSlice() override;
  void init() override;
  void compute_vector() override;
  void compute_array() override;

 private:
  int me;
  int nstart, nstop, nskip, nvalues;
  int *which, *argindex, *value2index;
  char **ids;

  void extract_one(int, double *, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
