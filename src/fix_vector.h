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

#ifdef FIX_CLASS
// clang-format off
FixStyle(vector,FixVector);
// clang-format on
#else

#ifndef LMP_FIX_VECTOR_H
#define LMP_FIX_VECTOR_H

#include "fix.h"

namespace LAMMPS_NS {

class FixVector : public Fix {
 public:
  FixVector(class LAMMPS *, int, char **);
  ~FixVector() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void end_of_step() override;
  double compute_vector(int) override;
  double compute_array(int, int) override;

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

  bigint nextstep, initialstep;

  bigint ncount;       // # of values processed and stored into growing vector or array
  bigint ncountmax;    // max # of values vector/array can hold
  bigint nmaxval;      // maximum allowed number of values
  bigint nindex;       // start index of data, may wrap around

  double *vector;
  double **array;
};
}    // namespace LAMMPS_NS
#endif
#endif
