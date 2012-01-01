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

FixStyle(bond/create,FixBondCreate)

#else

#ifndef LMP_FIX_BOND_CREATE_H
#define LMP_FIX_BOND_CREATE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBondCreate : public Fix {
 public:
  FixBondCreate(class LAMMPS *, int, char **);
  ~FixBondCreate();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void post_integrate();
  void post_integrate_respa(int, int);

  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  void grow_arrays(int);
  void copy_arrays(int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  double compute_vector(int);
  double memory_usage();

 private:
  int me;
  int iatomtype,jatomtype;
  int btype,seed;
  int imaxbond,jmaxbond;
  int inewtype,jnewtype;
  double cutsq,fraction;

  int createcount,createcounttotal;   // bond formation stats

  int nmax;
  int *bondcount;        // count of created bonds this atom is part of
  int *partner;          // ID of preferred atom for this atom to bond to
  double *distsq;        // distance to preferred bond partner
  double *probability;   // random # to use in decision to form bond

  class RanMars *random;
  class NeighList *list;
  int countflag,commflag;
  int nlevels_respa;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid atom type in fix bond/create command

Self-explanatory.

E: Invalid bond type in fix bond/create command

Self-explanatory.

E: Cannot use fix bond/create with non-molecular systems

Self-explanatory.

E: Inconsistent iparam/jparam values in fix bond/create command

If itype and jtype are the same, then their maxbond and newtype 
settings must also be the same.

E: Fix bond/create cutoff is longer than pairwise cutoff

This is not allowed because bond creation is done using the
pairwise neighbor list.

E: Fix bond/create requires special_bonds lj = 0,1,1

Self-explanatory.

E: Fix bond/create requires special_bonds coul = 0,1,1

Self-explanatory.

W: Created bonds will not create angles, dihedrals, or impropers

See the doc page for fix bond/create for more info on this
restriction.

E: Could not count initial bonds in fix bond/create

Could not find one of the atoms in a bond on this processor.

E: New bond exceeded bonds per atom in fix bond/create

See the read_data command for info on setting the "extra bond per
atom" header value to allow for additional bonds to be formed.

E: New bond exceeded special list size in fix bond/create

See the special_bonds extra command for info on how to leave space in
the special bonds list to allow for additional bonds to be formed.

*/
