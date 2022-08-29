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
FixStyle(STORE/GLOBAL,FixStoreGlobal);
// clang-format on
#else

#ifndef LMP_FIX_STORE_GLOBAL_H
#define LMP_FIX_STORE_GLOBAL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixStoreGlobal : public Fix {
 public:
  int nrow, ncol;     // copy of n1,n2 for array for classes to access
  double *vstore;     // vector storage
  double **astore;    // array storage

  FixStoreGlobal(class LAMMPS *, int, char **);
  ~FixStoreGlobal() override;
  int setmask() override;
  void reset_global(int, int);

  void write_restart(FILE *) override;
  void restart(char *) override;

  double memory_usage() override;

 private:
  int vecflag;      // 1 if ncol=1 or nvalues=1
  int arrayflag;    // 1 if ncol > 1

  int n1, n2;      // size of 3d dims of data struct
  double *rbuf;    // restart buffer for GLOBAL vec/array/tensor
};

}    // namespace LAMMPS_NS

#endif
#endif
