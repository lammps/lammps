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

FixStyle(restrain,FixRestrain)

#else

#ifndef LMP_FIX_RESTRAIN_H
#define LMP_FIX_RESTRAIN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRestrain : public Fix {
 public:
  FixRestrain(class LAMMPS *, int, char **);
  ~FixRestrain();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();

 private:
  int nlevels_respa;
  int nrestrain,maxrestrain;
  int *rstyle;
  int **ids;
  double *kstart,*kstop,*target;
  double *cos_target,*sin_target;
  double energy,energy_all;

  void restrain_bond(int);
  void restrain_angle(int);
  void restrain_dihedral(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix restrain requires an atom map, see atom_modify

Self-explanatory.

E: Restrain atoms %d %d missing on proc %d at step %ld

The 2 atoms in a restrain bond specified by the fix restrain
command are not all accessible to a processor.  This probably means an
atom has moved too far.

E: Restrain atoms %d %d %d missing on proc %d at step %ld

The 3 atoms in a restrain angle specified by the fix restrain
command are not all accessible to a processor.  This probably means an
atom has moved too far.

E: Restrain atoms %d %d %d %d missing on proc %d at step %ld

The 4 atoms in a restrain dihedral specified by the fix restrain
command are not all accessible to a processor.  This probably means an
atom has moved too far.

W: Restrain problem: %d %ld    

UNDOCUMENTED

U: Restrain problem: %d %ld %d %d %d %d

Conformation of the 4 listed dihedral atoms is extreme; you may want
to check your simulation geometry.

*/
