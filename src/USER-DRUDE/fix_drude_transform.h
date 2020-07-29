/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(drude/transform/direct,FixDrudeTransform<false>)
FixStyle(drude/transform/inverse,FixDrudeTransform<true>)

#else

#ifndef LMP_FIX_DRUDE_TRANSFORM_H
#define LMP_FIX_DRUDE_TRANSFORM_H

#include "fix.h"

namespace LAMMPS_NS {

template <bool inverse>
class FixDrudeTransform : public Fix {
 public:
  FixDrudeTransform<inverse>(class LAMMPS *, int, char **);
  ~FixDrudeTransform<inverse>();
  int setmask();
  void init();
  void setup(int vflag);
  void reduced_to_real();
  void real_to_reduced();
  void initial_integrate(int vflag);
  void final_integrate();
  int pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc);
  void unpack_forward_comm(int n, int first, double *buf);
 protected:
  double * mcoeff;
  class FixDrude * fix_drude;
};

}

#endif
#endif

