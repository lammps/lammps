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

FixStyle(reax/bonds,FixReaxBonds)

#else

#ifndef LMP_FIX_REAX_BONDS_H
#define LMP_FIX_REAX_BONDS_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixReaxBonds : public Fix {
 public:
  FixReaxBonds(class LAMMPS *, int, char **);
  ~FixReaxBonds();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();

 private:
  int me;
  int nfreq;
  FILE *fp;

  void OutputReaxBonds(bigint, FILE*);
  int nint(const double&);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot open fix reax/bonds file %s

The output file for the fix reax/bonds command cannot be opened.
Check that the path and name are correct.

E: Cannot use fix reax/bonds without pair_style reax

Self-explantory.

E: Fix reax/bonds numbonds > nsbmax_most

The limit of the number of bonds expected by the ReaxFF force field
was exceeded.

*/
