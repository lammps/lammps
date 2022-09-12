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
FixStyle(latte,FixLatte);
// clang-format on
#else

#ifndef LMP_FIX_LATTE_H
#define LMP_FIX_LATTE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLatte : public Fix {
 public:
  FixLatte(class LAMMPS *, int, char **);
  ~FixLatte() override;
  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void setup(int) override;
  void min_setup(int) override;
  void setup_pre_reverse(int, int) override;
  void initial_integrate(int) override;
  void pre_reverse(int, int) override;
  void post_force(int) override;
  void min_post_force(int) override;
  void final_integrate() override;
  void reset_dt() override;
  double compute_scalar() override;
  double memory_usage() override;

 protected:
  int coulomb, pbcflag, pe_peratom, virial_global, virial_peratom, neighflag;
  int exclude, excludebit;
  int eflag_caller;
  char *id_pe,*id_exclude;
  int *exclusion_group_ptr;
  int setupflag, newsystem;
  bigint natoms_last;

  int flags_latte[6];

  int nmax;
  double *qpotential;
  double **flatte;
  double latte_energy;

  class NeighList *list;
  class Compute *c_pe;

  void latte_wrapper_all();
  void latte_wrapper_exclude();
};

}    // namespace LAMMPS_NS

#endif
#endif
