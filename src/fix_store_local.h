/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(STORE/LOCAL,FixStoreLocal);
// clang-format on
#else

#ifndef LMP_FIX_STORE_LOCAL_H
#define LMP_FIX_STORE_LOCAL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixStoreLocal : public Fix {
 public:
  FixStoreLocal(class LAMMPS *, int, char **);
  ~FixStoreLocal() override;
  int setmask() override;
  void post_force(int) override;
  double memory_usage() override;
  void add_data(double *, int, int);
  int nvalues;

 private:
  int nmax;

  double *vector;
  double **array;

  int ncount;
  int nreset;

  void reallocate(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
