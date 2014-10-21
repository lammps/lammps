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

  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  void grow_arrays(int);
  void copy_arrays(int, int, int);
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
  int atype,dtype,itype;
  int angleflag,dihedralflag,improperflag;
  int overflow;
  tagint lastcheck;

  int *bondcount;
  int createcount,createcounttotal;
  int nmax;
  tagint *partner,*finalpartner;
  double *distsq,*probability;

  int ncreate,maxcreate;
  tagint **created;

  tagint *copy;

  class RanMars *random;
  class NeighList *list;
  
  int countflag,commflag;
  int nlevels_respa;
  int nangles,ndihedrals,nimpropers;

  void check_ghosts();
  void update_topology();
  void rebuild_special(int);
  void create_angles(int);
  void create_dihedrals(int);
  void create_impropers(int);
  int dedup(int, int, tagint *);

  // DEBUG

  void print_bb();
  void print_copy(const char *, tagint, int, int, int, int *);
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

Only systems with bonds that can be changed can be used.  Atom_style
template does not qualify.

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
