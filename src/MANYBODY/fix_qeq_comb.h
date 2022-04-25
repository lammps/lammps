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
FixStyle(qeq/comb,FixQEQComb);
// clang-format on
#else

#ifndef LMP_FIX_QEQ_COMB_H
#define LMP_FIX_QEQ_COMB_H

#include "fix.h"

namespace LAMMPS_NS {

class FixQEQComb : public Fix {
 public:
  FixQEQComb(class LAMMPS *, int, char **);
  ~FixQEQComb() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  double memory_usage() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

  void min_post_force(int) override;

 protected:
  int me, firstflag;
  double precision;
  int ilevel_respa;
  bigint ngroup;
  FILE *fp;

  class PairComb *comb;
  class PairComb3 *comb3;
  int nmax;
  double *qf, *q1, *q2;
};

}    // namespace LAMMPS_NS

#endif
#endif
