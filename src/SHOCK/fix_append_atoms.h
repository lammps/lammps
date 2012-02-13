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

E: Only zhi currently implemented for fix append/atoms

Self-explanatory.

E: Append boundary must be shrink/minimum

The boundary style of the face where atoms are added
must be of type m (shrink/minimum). 

E: Bad fix ID in fix append/atoms command

The value of the fix_id for keyword spatial must start with the suffix
f_.

E: Cannot use append/atoms in periodic dimension

The boundary style of the face where atoms are added can not be of
type p (periodic).

E: Cannot append atoms to a triclinic box

The simulation box must be defined with edges alligned with the
Cartesian axes.

E: Use of fix append/atoms with undefined lattice

A lattice must be defined before using this fix.

E: Fix ID for fix ave/spatial does not exist

Self-explanatory.

E: Must define lattice to append/atoms

A lattice must be defined before using this fix.

*/
