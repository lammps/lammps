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
  ~FixExternal();
  int setmask();
  void init();
  void setup(int);
  void setup_pre_reverse(int, int);
  void min_setup(int);
  void pre_reverse(int, int);
  void post_force(int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);

  void set_energy_global(double);
  void set_virial_global(double *);
  void set_energy_peratom(double *);
  void set_virial_peratom(double **);
  void set_vector_length(int);
  void set_vector(int, double);

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

  typedef void (*FnPtr)(void *, bigint, int, tagint *, double **, double **);
  void set_callback(FnPtr, void *);

  void *extract(const char *, int &);

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

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix external callback function not set

This must be done by an external program in order to use this fix.

E: Invalid set_vector index in fix external

UNDOCUMENTED

*/
