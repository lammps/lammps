/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Pair zero is a dummy pair interaction useful for requiring a
   force cutoff distance in the absence of pair-interactions or
   with hybrid/overlay if a larger force cutoff distance is required.

   This can be used in conjunction with bond/create to create bonds
   that are longer than the cutoff of a given force field, or to
   calculate radial distribution functions for models without
   pair interactions.

------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(python/move,FixPythonMove);
// clang-format on
#else

#ifndef LMP_FIX_PYTHON_MOVE_H
#define LMP_FIX_PYTHON_MOVE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPythonMove : public Fix {
 public:
  FixPythonMove(LAMMPS *lmp, int narg, char **arg);
  virtual ~FixPythonMove();

  int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual void initial_integrate_respa(int, int, int);
  virtual void final_integrate_respa(int, int);
  virtual void reset_dt();

 protected:
  void *py_move;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Fix python/integrate requires fully qualified class name

UNDOCUMENTED

E: Loading python integrator module failure

UNDOCUMENTED

E: Could not find integrator class in module'

UNDOCUMENTED

E: Could not instantiate instance of integrator class'

UNDOCUMENTED

E: Could not find 'init' method'

UNDOCUMENTED

E: Could not find 'initial_integrate' method'

UNDOCUMENTED

E: Could not find 'final_integrate' method'

UNDOCUMENTED

E: Could not find 'initial_integrate_respa' method'

UNDOCUMENTED

E: Could not find 'final_integrate_respa' method'

UNDOCUMENTED

E: Could not find 'reset_dt' method'

UNDOCUMENTED

U: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
