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
FixStyle(ave/atom,FixAveAtom);
// clang-format on
#else

#ifndef LMP_FIX_AVE_ATOM_H
#define LMP_FIX_AVE_ATOM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAveAtom : public Fix {
 public:
  FixAveAtom(class LAMMPS *, int, char **);
  ~FixAveAtom() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void end_of_step() override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

 private:
  struct value_t {
    int which;         // type of data: COMPUTE, FIX, VARIABLE
    int argindex;      // 1-based index if data is vector, else 0
    std::string id;    // compute/fix/variable ID
    union {
      class Compute *c;
      class Fix *f;
      int v;
    } val;
  };
  std::vector<value_t> values;

  int nrepeat, irepeat;
  bigint nvalid, nvalid_last;
  double **array;

  bigint nextvalid();
};

}    // namespace LAMMPS_NS

#endif
#endif
