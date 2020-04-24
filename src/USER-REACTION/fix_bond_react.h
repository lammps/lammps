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

  enum {MAXLINE=256};

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
  int newton_bond;
  int nreacts;
  int *nevery;
  FILE *fp;
  int *iatomtype,*jatomtype;
  int *seed;
  double **cutsq,*fraction;
  int *max_rxn,*nlocalskips,*nghostlyskips;
  tagint lastcheck;
  int stabilization_flag;
  int custom_exclude_flag;
  int *stabilize_steps_flag;
  int *update_edges_flag;
  int nconstraints;
  int narrhenius;
  double **constraints;
  int status;
  int *groupbits;

  int rxnID; // integer ID for identifying current bond/react
  char **rxn_name; // name of reaction
  int *reaction_count;
  int *reaction_count_total;
  int nmax; // max num local atoms
  int max_natoms; // max natoms in a molecule template
  tagint *partner,*finalpartner;
  double **distsq,*probability;
  int *ncreate;
  int maxcreate;
  int allncreate;
  tagint ***created;

  class Molecule *onemol; // pre-reacted molecule template
  class Molecule *twomol; // post-reacted molecule template
  Fix *fix1;              // nve/limit used to relax reaction sites
  Fix *fix2;              // properties/atom used to indicate 1) relaxing atoms
                          //                                  2) to which 'react' atom belongs
  Fix *fix3;              // property/atom used for system-wide thermostat
  class RanMars **random; // random number for 'prob' keyword
  class RanMars **rrhandom; // random number for Arrhenius constraint
  class NeighList *list;

  int *reacted_mol,*unreacted_mol;
  int *limit_duration; // indicates how long to relax
  char *nve_limit_xmax; // indicates max distance allowed to move when relaxing
  char *id_fix1; // id of internally created fix nve/limit
  char *id_fix2; // id of internally created fix per-atom properties
  char *id_fix3; // id of internally created 'stabilization group' per-atom property fix
  char *statted_id; // name of 'stabilization group' per-atom property
  char *master_group; // group containing relaxing atoms from all fix rxns
  char *exclude_group; // group for system-wide thermostat

  int countflag,commflag;
  int nlevels_respa;

  void superimpose_algorithm(); // main function of the superimpose algorithm

  int *ibonding,*jbonding;
  int *closeneigh; // indicates if bonding atoms of a rxn are 1-2, 1-3, or 1-4 neighbors
  int nedge,nequivalent,ncustom,ndelete,nchiral,nconstr; // # edge, equivalent, custom atoms in mapping file
  int attempted_rxn; // there was an attempt!
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
  int **custom_edges; // atoms in molecule templates with incorrect valences
  int **delete_atoms; // atoms in pre-reacted templates to delete
  int ***chiral_atoms; // pre-react chiral atoms. 1) flag 2) orientation 3-4) ordered atom types

  int **nxspecial,**onemol_nxspecial,**twomol_nxspecial; // full number of 1-4 neighbors
  tagint **xspecial,**onemol_xspecial,**twomol_xspecial; // full 1-4 neighbor list

  int pion,neigh,trace; // important indices for various loops. required for restore points
  int lcl_inst; // reaction instance
  tagint **glove; // 1st colmn: pre-reacted template, 2nd colmn: global IDs
  // for all mega_gloves and global_mega_glove: first row is the ID of bond/react
  tagint **local_mega_glove; // consolidation local of reaction instances
  tagint **ghostly_mega_glove; // consolidation nonlocal of reaction instances
  tagint **global_mega_glove; // consolidation (inter-processor) of gloves containing nonlocal atoms
  int *localsendlist; // indicates ghosts of other procs
  int local_num_mega; // num of local reaction instances
  int ghostly_num_mega; // num of ghostly reaction instances
  int global_megasize; // num of reaction instances in global_mega_glove
  int *pioneers; // during Superimpose Algorithm, atoms which have been assigned, but whose first neighbors haven't
  int glove_counter; // used to determine when to terminate Superimpose Algorithm

  void read(int);
  void EdgeIDs(char *,int);
  void Equivalences(char *,int);
  void CustomEdges(char *,int);
  void DeleteAtoms(char *,int);
  void ChiralCenters(char *,int);
  void Constraints(char *,int);

  void make_a_guess ();
  void neighbor_loop();
  void check_a_neighbor();
  void crosscheck_the_neighbor();
  void inner_crosscheck_loop();
  void ring_check();
  int check_constraints();
  double get_temperature();
  int get_chirality(double[12]); // get handedness given an ordered set of coordinates

  void open(char *);
  void readline(char *);
  void parse_keyword(int, char *, char *);
  void skip_lines(int, char *);
  int parse(char *, char **, int);

  void far_partner();
  void close_partner();
  void get_molxspecials();
  void find_landlocked_atoms(int);
  void glove_ghostcheck();
  void ghost_glovecast();
  void update_everything();
  void unlimit_bond();
  void limit_bond(int);
  void dedup_mega_gloves(int); //dedup global mega_glove
  virtual void write_restart(FILE *);
  virtual void restart(char *buf);

  struct Set {
    int nreacts;
    char rxn_name[MAXLINE];
    int reaction_count_total;
  };
  Set *set;

  // DEBUG

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

E: Bond/react: Cannot use fix bond/react with non-molecular systems

Only systems with bonds that can be changed can be used. Atom_style
template does not qualify.

E: Bond/react: Rmax cutoff is longer than pairwise cutoff

This is not allowed because bond creation is done using the pairwise
neighbor list.

E: Bond/react: Molecule template ID for fix bond/react does not exist

A valid molecule template must have been created with the molecule
command.

E: Bond/react: Reaction templates must contain the same number of atoms

There should be a one-to-one correspondence between atoms in the
pre-reacted and post-reacted templates, as specified by the map file.

E: Bond/react: Unknown section in map file

Please ensure reaction map files are properly formatted.

E: Bond/react: Invalid template atom ID in map file

Atom IDs in molecule templates range from 1 to the number of atoms in the template.

E or W: Bond/react: Atom affected by reaction %s too close to template edge

This means an atom which changes type or connectivity during the
reaction is too close to an 'edge' atom defined in the map file. This
could cause incorrect assignment of bonds, angle, etc. Generally, this
means you must include more atoms in your templates, such that there
are at least two atoms between each atom involved in the reaction and
an edge atom.

E: Bond/react: Fix bond/react needs ghost atoms from farther away

This is because a processor needs to map the entire unreacted molecule
template onto simulation atoms it knows about. The comm_modify cutoff
command can be used to extend the communication range.

E: Bond/react: A deleted atom cannot remain bonded to an atom that is not deleted

Self-explanatory.

E: Bond/react: First neighbors of chiral atoms must be of mutually different types

Self-explanatory.

E: Bond/react: Chiral atoms must have exactly four first neighbors

Self-explanatory.

E: Bond/react: Molecule template 'Coords' section required for chiralIDs keyword

The coordinates of atoms in the pre-reacted template are used to determine chirality.

E: Bond/react special bond generation overflow

The number of special bonds per-atom created by a reaction exceeds the
system setting. See the read_data or create_box command for how to
specify this value.

E: Bond/react topology/atom exceed system topology/atom

The number of bonds, angles etc per-atom created by a reaction exceeds
the system setting. See the read_data or create_box command for how to
specify this value.

*/
