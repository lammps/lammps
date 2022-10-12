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

#ifdef FIX_CLASS
// clang-format off
FixStyle(drude/transform/direct,FixDrudeTransform<false>);
FixStyle(drude/transform/inverse,FixDrudeTransform<true>);
// clang-format on
#else

#ifndef LMP_FIX_DRUDE_TRANSFORM_H
#define LMP_FIX_DRUDE_TRANSFORM_H

#include "fix.h"

namespace LAMMPS_NS {

template <bool inverse> class FixDrudeTransform : public Fix {
 public:
  FixDrudeTransform(class LAMMPS *, int, char **);
  ~FixDrudeTransform() override;
  int setmask() override;
  void init() override;
  void setup(int vflag) override;
  void reduced_to_real();
  void real_to_reduced();
  void initial_integrate(int vflag) override;
  void final_integrate() override;
  int pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc) override;
  void unpack_forward_comm(int n, int first, double *buf) override;

 protected:
  double *mcoeff;
  class FixDrude *fix_drude;
};

}    // namespace LAMMPS_NS

#endif
#endif
