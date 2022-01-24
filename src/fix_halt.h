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
FixStyle(halt,FixHalt);
// clang-format on
#else

#ifndef LMP_FIX_HALT_H
#define LMP_FIX_HALT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixHalt : public Fix {
 public:
  FixHalt(class LAMMPS *, int, char **);
  ~FixHalt();
  int setmask();
  void init();
  void end_of_step();
  void min_post_force(int);
  void post_run();

 private:
  int attribute, operation, eflag, msgflag, ivar;
  bigint nextstep, thisstep;
  double value, tratio;
  char *idvar;
  char *dlimit_path;

  double bondmax();
  double tlimit();
  double diskfree();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not find fix halt variable name

Self-explanatory.

E: Fix halt variable is not equal-style variable

Self-explanatory.

E: Invalid fix halt attribute

Self-explanatory.

E: Invalid fix halt operator

Self-explanatory.

E: Disk limit not supported by OS or illegal path

Self-explanatory.

W: Fix halt condition for fix-id %s met on step %ld with value %g

Self explanatory.

*/
