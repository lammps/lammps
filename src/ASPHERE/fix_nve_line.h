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
FixStyle(nve/line,FixNVELine);
// clang-format on
#else

#ifndef LMP_FIX_NVE_LINE_H
#define LMP_FIX_NVE_LINE_H

#include "fix_nve.h"

namespace LAMMPS_NS {

class FixNVELine : public FixNVE {
 public:
  FixNVELine(class LAMMPS *, int, char **);
  ~FixNVELine() {}
  int setmask();
  void init();
  void initial_integrate(int);
  void final_integrate();

 private:
  double MINUSPI, TWOPI;
  class AtomVecLine *avec;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix nve/line requires atom style line

Self-explanatory.

E: Fix nve/line can only be used for 2d simulations

Self-explanatory.

E: Fix nve/line requires line particles

Self-explanatory.

*/
