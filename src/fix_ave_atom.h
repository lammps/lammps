/* ----------------------------------------------------------------------
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

FixStyle(ave/atom,FixAveAtom)

#else

#ifndef LMP_FIX_AVE_ATOM_H
#define LMP_FIX_AVE_ATOM_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixAveAtom : public Fix {
 public:
  FixAveAtom(class LAMMPS *, int, char **);
  ~FixAveAtom();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  void reset_timestep(bigint);

 private:
  int nvalues;
  int nrepeat,irepeat;
  bigint nvalid;
  int *which,*argindex,*value2index;
  char **ids;
  double **array;

  bigint nextvalid();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute ID for fix ave/atom does not exist

Self-explanatory.

E: Fix ave/atom compute does not calculate per-atom values

A compute used by fix ave/atom must generate per-atom values.

E: Fix ave/atom compute does not calculate a per-atom vector

A compute used by fix ave/atom must generate per-atom values.

E: Fix ave/atom compute does not calculate a per-atom array

Self-explanatory.

E: Fix ave/atom compute array is accessed out-of-range

Self-explanatory.

E: Fix ID for fix ave/atom does not exist

Self-explanatory.

E: Fix ave/atom fix does not calculate per-atom values

A fix used by fix ave/atom must generate per-atom values.

E: Fix ave/atom fix does not calculate a per-atom vector

A fix used by fix ave/atom must generate per-atom values.

E: Fix ave/atom fix does not calculate a per-atom array

Self-explanatory.

E: Fix ave/atom fix array is accessed out-of-range

Self-explanatory.

E: Fix for fix ave/atom not computed at compatible time

Fixes generate their values on specific timesteps.  Fix ave/atom is
requesting a value on a non-allowed timestep.

E: Variable name for fix ave/atom does not exist

Self-explanatory.

E: Fix ave/atom variable is not atom-style variable

A variable used by fix ave/atom must generate per-atom values.

E: Fix ave/atom missed timestep

You cannot reset the timestep to a value beyond where the fix
expects to next perform averaging.

*/
