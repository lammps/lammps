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

FixStyle(bond/swap,FixBondSwap)

#else

#ifndef LMP_FIX_BONDSWAP_H
#define LMP_FIX_BONDSWAP_H

#include "fix.h"
#include "pair.h"

namespace LAMMPS_NS {

class FixBondSwap : public Fix {
 public:
  FixBondSwap(class LAMMPS *, int, char **);
  ~FixBondSwap();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void post_integrate();
  int modify_param(int, char **);
  double compute_vector(int);
  double memory_usage();

 private:
  double fraction,cutsq;
  int nmax,tflag;
  int *alist;
  int naccept,foursome;
  int angleflag;
  char *id_temp;
  int *type;
  double **x;

  class NeighList *list;
  class Compute *temperature;
  class RanMars *random;

  double dist_rsq(int, int);
  double pair_eng(int, int);
  double bond_eng(int, int, int);
  double angle_eng(int, int, int, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use fix bond/swap with non-molecular systems

Only systems with bonds that can be changed can be used.  Atom_style
template does not qualify.

E: Must use atom style with molecule IDs with fix bond/swap

Self-explanatory.

E: Temperature ID for fix bond/swap does not exist

Self-explanatory.

E: Fix bond/swap requires pair and bond styles

Self-explanatory.

E: Pair style does not support fix bond/swap

The pair style does not have a single() function, so it can
not be invoked by fix bond/swap.

W: Fix bond/swap will ignore defined angles

See the doc page for fix bond/swap for more info on this
restriction.

E: Fix bond/swap cannot use dihedral or improper styles

These styles cannot be defined when using this fix.

E: Fix bond/swap requires special_bonds = 0,1,1

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
