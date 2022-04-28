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
FixStyle(store/force,FixStoreForce);
// clang-format on
#else

#ifndef LMP_FIX_STORE_FORCE_H
#define LMP_FIX_STORE_FORCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixStoreForce : public Fix {
 public:
  FixStoreForce(class LAMMPS *, int, char **);
  ~FixStoreForce() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double memory_usage() override;

 private:
  int nlevels_respa;
  int nmax;
  double **foriginal;    // stored force on atoms
};

}    // namespace LAMMPS_NS

#endif
#endif
