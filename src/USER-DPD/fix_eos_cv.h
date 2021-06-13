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
FixStyle(eos/cv,FixEOScv);
// clang-format on
#else

#ifndef LMP_FIX_EOS_CV_H
#define LMP_FIX_EOS_CV_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEOScv : public Fix {
 public:
  FixEOScv(class LAMMPS *, int, char **);
  virtual ~FixEOScv() {}
  int setmask();
  virtual void init();
  virtual void post_integrate();
  virtual void end_of_step();

 protected:
  double cvEOS;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E:  FixEOScv requires atom_style with internal temperature and energies (e.g. dpd)

Self-explanatory.

E: EOS cv must be > 0.0

The constant volume heat capacity must be larger than zero.

E: Internal temperature < zero

Self-explanatory.  EOS may not be valid under current simulation conditions.

*/
