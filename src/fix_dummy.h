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

FixStyle(DUMMY,FixDummy)

#else

#ifndef LMP_FIX_DUMMY_H
#define LMP_FIX_DUMMY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDummy : public Fix {
 public:
  FixDummy(class LAMMPS *, int, char **);
  virtual ~FixDummy() {}
  int setmask();

 protected:
  int initial_integrate_flag,final_integrate_flag;
  int pre_exchange_flag,pre_neighbor_flag;
  int pre_force_flag,post_force_flag;
  int end_of_step_flag;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
