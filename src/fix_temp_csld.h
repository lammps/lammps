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

FixStyle(temp/csld,FixTempCSLD)

#else

#ifndef LMP_FIX_TEMP_CSLD_H
#define LMP_FIX_TEMP_CSLD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTempCSLD : public Fix {
 public:
  FixTempCSLD(class LAMMPS *, int, char **);
  ~FixTempCSLD();
  int setmask();
  void init();
  void end_of_step();
  int modify_param(int, char **);
  void reset_target(double);
  virtual double compute_scalar();
  virtual void *extract(const char *, int &);

 private:
  double t_start,t_stop,t_period,t_target;
  double **vhold;
  double energy;
  int nmax,which;
  int tstyle,tvar;
  char *tstr;

  char *id_temp;
  class Compute *temperature;
  int tflag;

  class RanMars *random;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix temp/csld is not compatible with fix rattle or fix shake

These two commands cannot currently be used together with fix temp/csld.

E: Variable name for fix temp/csld does not exist

Self-explanatory.

E: Variable for fix temp/csld is invalid style

Only equal-style variables can be used.

E: Temperature ID for fix temp/csld does not exist

Self-explanatory.

W: Cannot thermostat atoms in rigid bodies

UNDOCUMENTED

E: Fix temp/csld variable returned negative temperature

Self-explanatory.

E: Could not find fix_modify temperature ID

The compute ID for computing temperature does not exist.

E: Fix_modify temperature ID does not compute temperature

The compute ID assigned to the fix must compute temperature.

W: Group for fix_modify temp != fix group

The fix_modify command is specifying a temperature computation that
computes a temperature on a different group of atoms than the fix
itself operates on.  This is probably not what you want to do.

*/
