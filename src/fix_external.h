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
FixStyle(external,FixExternal);
// clang-format on
#else

#ifndef LMP_FIX_EXTERNAL_H
#define LMP_FIX_EXTERNAL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixExternal : public Fix {
 public:
  double **fexternal;

  FixExternal(class LAMMPS *, int, char **);
  ~FixExternal() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void setup_pre_reverse(int, int) override;
  void min_setup(int) override;
  void pre_reverse(int, int) override;
  void post_force(int) override;
  void min_post_force(int) override;
  double compute_scalar() override;
  double compute_vector(int) override;

  void set_energy_global(double);
  void set_virial_global(double *);
  void set_energy_peratom(double *);
  void set_virial_peratom(double **);
  void set_vector_length(int);
  void set_vector(int, double);

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

  typedef void (*FnPtr)(void *, bigint, int, tagint *, double **, double **);
  void set_callback(FnPtr, void *);

  void *extract(const char *, int &) override;

 private:
  int mode, ncall, napply, eflag_caller;
  FnPtr callback;
  void *ptr_caller;
  double user_energy;
  double user_virial[6];
  double *caller_vector;
};

}    // namespace LAMMPS_NS

#endif
#endif
