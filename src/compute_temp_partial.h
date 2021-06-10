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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(temp/partial,ComputeTempPartial);
// clang-format on
#else

#ifndef LMP_COMPUTE_TEMP_PARTIAL_H
#define LMP_COMPUTE_TEMP_PARTIAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempPartial : public Compute {
 public:
  ComputeTempPartial(class LAMMPS *, int, char **);
  virtual ~ComputeTempPartial();
  void init() {}
  void setup();
  double compute_scalar();
  void compute_vector();

  int dof_remove(int);
  void remove_bias(int, double *);
  void remove_bias_thr(int, double *, double *);
  void remove_bias_all();
  void reapply_bias_all();
  void restore_bias(int, double *);
  void restore_bias_thr(int, double *, double *);
  void restore_bias_all();
  double memory_usage();

 protected:
  int xflag, yflag, zflag;
  double tfactor;

  void dof_compute();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute temp/partial cannot use vz for 2d systemx

Self-explanatory.

E: Temperature compute degrees of freedom < 0

This should not happen if you are calculating the temperature
on a valid set of atoms.

*/
