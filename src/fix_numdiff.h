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

FixStyle(numdiff,FixNumDiff)

#else

#ifndef LMP_FIX_NUMDIFF_H
#define LMP_FIX_NUMDIFF_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNumDiff : public Fix {
 public:
  FixNumDiff(class LAMMPS *, int, char **);
  ~FixNumDiff();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double memory_usage();

private:
  double delta;
  int maxatom;
  int ilevel_respa;

  int pair_compute_flag;            // 0 if pair->compute is skipped
  int kspace_compute_flag;          // 0 if kspace->compute is skipped

  char *id_pe;
  class Compute *pe;

  double **numdiff_forces;          // finite diff forces
  double **temp_x;                  // original coords
  double **temp_f;                  // original forces

  void calculate_forces();
  void displace_atoms(int, int, int);
  void restore_atoms(int, int);
  double update_energy();
  void force_clear(double **);
  void reallocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute ID for fix numdiff does not exist

Self-explanatory.

*/
