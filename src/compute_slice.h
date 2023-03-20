/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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
  struct value_t {
    int which;
    int argindex;
    std::string id;
    union {
      class Compute *c;
      class Fix *f;
      int v;
    } val;
  };
  std::vector<value_t> values;
  int nstart, nstop, nskip;

  void extract_one(int, double *, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
