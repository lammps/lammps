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

/* ----------------------------------------------------------------------
   Contributing Author: Jacob Gissinger (jacob.gissinger@colorado.edu)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(bond/react,FixBondReact)

#else

#ifndef LMP_FIX_BOND_REACT_H
#define LMP_FIX_BOND_REACT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBondReact : public Fix {
 public:
  FixBondReact(class LAMMPS *, int, char **);
  ~FixBondReact();
  int setmask();
  void post_constructor();
  void init();
  void init_list(int, class NeighList *);
  void post_integrate();
  void post_integrate_respa(int, int);

  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double compute_vector(int);
  double memory_usage();

 private:
  int me,nprocs;
  int nreacts;
  int *nevery;
  FILE *fp;
  int *iatomtype,*jatomtype;
  int *seed;
  double *cutsq,*fraction;
  tagint lastcheck;
  int stabilization_flag;
  int *stabilize_steps_flag;
  int status;
  int *groupbits;

  int rxnID; // integer ID for identifying current bond/react
  int *reaction_count;
  int *reaction_count_total;
  int nmax; // max num local atoms
  int max_natoms; // max natoms in a molecule template
  tagint *partner,*finalpartner;
  double *distsq,*probability;
  int *ncreate;
  int maxcreate;
  int allncreate;
  tagint ***created;
  int *local_ncreate;

  class Molecule *onemol; // pre-reacted molecule template
  class Molecule *twomol; // post-reacted molecule template
  class FixNVELimit *fix1; // used to relax reaction sites
  class FixPropertyAtom *fix2; // used to indicate relaxing atoms
  class FixPropertyAtom *fix3; // used to indicate system-wide thermostat
  class FixPropertyAtom *fix4; // indicates to which 'react' atom belongs
  class RanMars **random;
  class NeighList *list;

  int *reacted_mol,*unreacted_mol;
  int *limit_duration; // indicates how long to relax
  char *nve_limit_xmax; // indicates max distance allowed to move when relaxing
  char *id_fix1; // id of internally created fix nve/limit
  char *id_fix2; // id of internally created fix per-atom property (recently reacted)
  char *id_fix3; // id of internally created fix per-atom property (system-wide thermostat)
  char *id_fix4; // id of internally created fix per-atom property (ID of 'react' argument)
  char *master_group; // group containing relaxing atoms from all fix rxns
  char *exclude_group; // group for system-wide thermostat

  int countflag,commflag;
  int nlevels_respa;

  void superimpose_algorithm(); // main function of the superimpose algorithm

  int *ibonding,*jbonding;
  int nedge,nequivalent; // number of edge, equivalent atoms in mapping file
  int attempted_rxn; // there was an attempt!
  int ghostcheck_flag; // idicates whether a reaction instances contains a nonlocal atom
  int this_rxn_count; // num of local reaction occurrences
  int *local_rxn_count;
  int *ghostly_rxn_count;
  int avail_guesses; // num of restore points available
  int *guess_branch; // used when there is more than two choices when guessing
  int **restore_pt; // contains info about restore points
  tagint **restore; // contaings info about restore points
  int *pioneer_count; // counts pioneers

  int **edge; // atoms in molecule templates with incorrect valences
  int ***equivalences; // relation between pre- and post-reacted templates
  int ***reverse_equiv; // re-ordered equivalences
  int **landlocked_atoms; // all atoms at least three bonds away from edge atoms

  int pion,neigh,trace; // important indices for various loops. required for restore points
  int lcl_inst; // reaction instance
  tagint **glove; // 1st colmn: pre-reacted template, 2nd colmn: global IDs
  // for all mega_gloves and global_mega_glove: first row is the ID of bond/react
  tagint **local_mega_glove; // consolidation local of reaction instances
  tagint **ghostly_mega_glove; // consolidation nonlocal of reaction instances
  tagint **global_mega_glove; // consolidation (inter-processor) of gloves containing nonlocal atoms
  int local_num_mega; // num of local reaction instances
  int ghostly_num_mega; // num of ghostly reaction instances
  int global_megasize; // num of reaction instances in global_mega_glove
  int *pioneers; // during Superimpose Algorithm, atoms which have been assigned, but whose first neighbors haven't
  int glove_counter; // used to determine when to terminate Superimpose Algorithm

  void read(int);
  void EdgeIDs(char *,int);
  void Equivalences(char *,int);

  void make_a_guess ();
  void neighbor_loop();
  void check_a_neighbor();
  void crosscheck_the_neighbor();
  void inner_crosscheck_loop();
  void ring_check();

  void open(char *);
  void readline(char *);
  void parse_keyword(int, char *, char *);
  void skip_lines(int, char *);
  int parse(char *, char **, int);

  void find_landlocked_atoms(int);
  void glove_ghostcheck();
  void ghost_glovecast();
  void update_everything();
  void unlimit_bond();
  void limit_bond(int);
  void dedup_mega_gloves(int); //dedup global mega_glove

  // DEBUG (currently obsolete)

  void print_bb();

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid exclude group name

Exclude group name should not previously be defined.

E: Cannot use fix bond/react with non-molecular systems

Only systems with bonds that can be changed can be used.  Atom_style
template does not qualify.

E: Fix bond/react cutoff is longer than pairwise cutoff

This is not allowed because bond creation is done using the
pairwise neighbor list.

E: Molecule template ID for fix bond/react does not exist

A valid molecule template must have been created with the molecule command.

E: Superimpose file errors:

Please ensure superimpose file is properly formatted.

E: Atom affected by reaction too close to template edge

This means an atom which changes type during the reaction is too close
to an 'edge' atom defined in the superimpose file. This could cause incorrect
assignment of bonds, angle, etc. Generally, this means you must include
more atoms in your templates, such that there are at least two atoms
between each atom involved in the reaction and an edge atom.

E: Fix bond/react needs ghost atoms from farther away

This is because a processor needs to superimpose the entire unreacted
molecule template onto simulation atoms it can 'see.' The comm_modify cutoff
command can be used to extend the communication range.

E: Excessive iteration of superimpose algorithm

You may have discovered a bug! But first, please double check that your
molecule template atom types, bond types, etc. are consistent with your simulation,
and that all atoms affected by a reaction are sufficently separated from edge atoms.
If this issue persists, please contact the developer.

*/
