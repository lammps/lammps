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
FixStyle(pair,FixPair);
// clang-format on
#else

#ifndef LMP_FIX_PAIR_H
#define LMP_FIX_PAIR_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPair : public Fix {
 public:
  FixPair(class LAMMPS *, int, char **);
  ~FixPair() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void setup_pre_force(int) override;
  void pre_force(int) override;
  void min_pre_force(int) override;
  void post_force(int) override;
  void min_post_force(int) override;

  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

  double memory_usage() override;

 private:
  int nevery, nfield, ncols;
  bigint lasttime;
  char *pairname;
  char **fieldname, **triggername;
  int *trigger;
  int **triggerptr;

  class Pair *pstyle;
  double *vector;
  double **array;

  void query_pstyle(LAMMPS *lmp);
};

}    // namespace LAMMPS_NS

#endif
#endif
