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
FixStyle(enforce2d,FixEnforce2D);
// clang-format on
#else

#ifndef LMP_FIX_ENFORCE2D_H
#define LMP_FIX_ENFORCE2D_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEnforce2D : public Fix {
 public:
  FixEnforce2D(class LAMMPS *, int, char **);
  ~FixEnforce2D();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);

 protected:
  int nfixlist;
  class Fix **flist;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use fix enforce2d with 3d simulation

Self-explanatory.

E: Fix enforce2d must be defined after fix %s

UNDOCUMENTED

*/
