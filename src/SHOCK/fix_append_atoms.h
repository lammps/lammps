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

FixStyle(append_atoms,FixAppendAtoms)

#else

#ifndef FIX_APPEND_ATOMS_H
#define FIX_APPEND_ATOMS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAppendAtoms : public Fix {
 public:
  FixAppendAtoms(class LAMMPS *, int, char **);
  ~FixAppendAtoms();
  int setmask();
  void setup(int);
  void pre_exchange();
  void initial_integrate(int);
  void post_force(int);

 private:
  int get_spatial();
  int spatflag, xloflag, xhiflag, yloflag, yhiflag, zloflag, zhiflag;
  int ranflag, tempflag, xseed, tseed;
  double ranx, rany, ranz, t_target, t_period, t_extent;
  class RanMars *randomx;
  class RanMars *randomt;
  int scaleflag, freq;
  int *basistype, nbasis;
  int advance, advance_sum;
  double size, spatlead;
  char *spatialid;
  double tfactor;
  double *gfactor1,*gfactor2;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Only zhi currently implemented for append_atoms

UNDOCUMENTED

E: Append boundary must be shrink/minimum

UNDOCUMENTED

E: Only zhi currently implemented for append_atom

UNDOCUMENTED

E: Bad fix ID in fix append_atoms command

UNDOCUMENTED

E: Cannot use append_atoms in periodic dimension

UNDOCUMENTED

E: Cannot append atoms to a triclinic box

UNDOCUMENTED

E: Use of fix append_atoms with undefined lattice

UNDOCUMENTED

E: Fix ID for fix ave/spatial does not exist

Self-explanatory.

E: Must define lattice to append_atoms

UNDOCUMENTED

U: must define lattice to append_atoms

UNDOCUMENTED

*/
