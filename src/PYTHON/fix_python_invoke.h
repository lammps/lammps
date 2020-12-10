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

FixStyle(python/invoke,FixPythonInvoke)
FixStyle(python,FixPythonInvoke)

#else

#ifndef LMP_FIX_PYTHON_INVOKE_H
#define LMP_FIX_PYTHON_INVOKE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPythonInvoke : public Fix {
 public:
  FixPythonInvoke(class LAMMPS *, int, char **);
  virtual ~FixPythonInvoke() {}
  int setmask();
  virtual void end_of_step();
  virtual void post_force(int);

 private:
  void * pFunc;
  int selected_callback;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Unsupported callback name for fix python/invoke

UNDOCUMENTED

E: Could not initialize embedded Python

UNDOCUMENTED

E: Could not find Python function

UNDOCUMENTED

*/
