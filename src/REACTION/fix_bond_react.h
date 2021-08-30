/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing Author: Jacob Gissinger (jacob.r.gissinger@gmail.com)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(bond/react,FixBondReact);
// clang-format on
#else

#ifndef LMP_FIX_BOND_REACT_H
#define LMP_FIX_BOND_REACT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBondReact : public Fix {
 public:
  enum { MAXLINE = 256 };    // max length of line read from files
  enum { MAXCONIDS = 4 };    // max # of IDs used by any constraint
  enum { MAXCONPAR = 5 };    // max # of constraint parameters

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
  int me, nprocs;
  int newton_bond;
  int nreacts;
  int *nevery;
  FILE *fp;
  int *iatomtype, *jatomtype;
  int *seed;
  double **cutsq, *fraction;
  int *max_rxn, *nlocalskips, *nghostlyskips;
  tagint lastcheck;
  int stabilization_flag;
  int reset_mol_ids_flag;
  int custom_exclude_flag;
  int *stabilize_steps_flag;
  int *custom_charges_fragid;
  int *create_atoms_flag;
  int *modify_create_fragid;
  double *overlapsq;
  int *molecule_keyword;
  int maxnconstraints;
  int *nconstraints;
  char **constraintstr;
  int nrxnfunction;
  std::vector<std::string> rxnfunclist;
  int narrhenius;
  int **var_flag, **var_id;    // for keyword values with variable inputs
  int status;
  int *groupbits;

  int rxnID;          // integer ID for identifying current bond/react
  char **rxn_name;    // name of reaction
  int *reaction_count;
  int *reaction_count_total;
  int nmax;          // max num local atoms
  int max_natoms;    // max natoms in a molecule template
  tagint *partner, *finalpartner;
  double **distsq;
  int *nattempt;
  int maxattempt;
  int allnattempt;
  tagint ***attempt;

  class Molecule *onemol;      // pre-reacted molecule template
  class Molecule *twomol;      // post-reacted molecule template
  Fix *fix1;                   // nve/limit used to relax reaction sites
  Fix *fix2;                   // properties/atom used to indicate 1) relaxing atoms
                               //                                  2) to which 'react' atom belongs
  Fix *fix3;                   // property/atom used for system-wide thermostat
  class RanMars **random;      // random number for 'prob' keyword
  class RanMars **rrhandom;    // random number for Arrhenius constraint
  class NeighList *list;
  class ResetMolIDs *reset_mol_ids;    // class for resetting mol IDs

  int *reacted_mol, *unreacted_mol;
  int *limit_duration;     // indicates how long to relax
  char *nve_limit_xmax;    // indicates max distance allowed to move when relaxing
  char *id_fix1;           // id of internally created fix nve/limit
  char *id_fix2;           // id of internally created fix per-atom properties
  char *id_fix3;           // id of internally created 'stabilization group' per-atom property fix
  char *statted_id;        // name of 'stabilization group' per-atom property
  char *master_group;      // group containing relaxing atoms from all fix rxns
  char *exclude_group;     // group for system-wide thermostat

  int countflag, commflag;
  int nlevels_respa;

  void superimpose_algorithm();    // main function of the superimpose algorithm

  int *ibonding, *jbonding;
  int *closeneigh;    // indicates if bonding atoms of a rxn are 1-2, 1-3, or 1-4 neighbors
  int nedge, nequivalent, ndelete, ncreate, nchiral;    // # edge, equivalent atoms in mapping file
  int attempted_rxn;                                    // there was an attempt!
  int *local_rxn_count;
  int *ghostly_rxn_count;
  int avail_guesses;     // num of restore points available
  int *guess_branch;     // used when there is more than two choices when guessing
  int **restore_pt;      // contains info about restore points
  tagint **restore;      // contaings info about restore points
  int *pioneer_count;    // counts pioneers

  int **edge;                // atoms in molecule templates with incorrect valences
  int ***equivalences;       // relation between pre- and post-reacted templates
  int ***reverse_equiv;      // re-ordered equivalences
  int **landlocked_atoms;    // all atoms at least three bonds away from edge atoms
  int **custom_charges;      // atoms whose charge should be updated
  int **delete_atoms;        // atoms in pre-reacted templates to delete
  int **create_atoms;        // atoms in post-reacted templates to create
  int ***chiral_atoms;    // pre-react chiral atoms. 1) flag 2) orientation 3-4) ordered atom types

  int **nxspecial, **onemol_nxspecial, **twomol_nxspecial;    // full number of 1-4 neighbors
  tagint **xspecial, **onemol_xspecial, **twomol_xspecial;    // full 1-4 neighbor list

  int pion, neigh, trace;    // important indices for various loops. required for restore points
  int lcl_inst;              // reaction instance
  tagint **glove;            // 1st colmn: pre-reacted template, 2nd colmn: global IDs
  // for all mega_gloves and global_mega_glove: first row is the ID of bond/react
  tagint **local_mega_glove;      // consolidation local of reaction instances
  tagint **ghostly_mega_glove;    // consolidation nonlocal of reaction instances
  tagint **global_mega_glove;    // consolidation (inter-processor) of gloves containing nonlocal atoms
  int *localsendlist;        // indicates ghosts of other procs
  int local_num_mega;        // num of local reaction instances
  int ghostly_num_mega;      // num of ghostly reaction instances
  int global_megasize;       // num of reaction instances in global_mega_glove
  int *pioneers;    // during Superimpose Algorithm, atoms which have been assigned, but whose first neighbors haven't
  int glove_counter;    // used to determine when to terminate Superimpose Algorithm

  void read(int);
  void EdgeIDs(char *, int);
  void Equivalences(char *, int);
  void DeleteAtoms(char *, int);
  void CreateAtoms(char *, int);
  void CustomCharges(int, int);
  void ChiralCenters(char *, int);
  void ReadConstraints(char *, int);
  void readID(char *, int, int, int);

  void make_a_guess();
  void neighbor_loop();
  void check_a_neighbor();
  void crosscheck_the_neighbor();
  void inner_crosscheck_loop();
  int ring_check();
  int check_constraints();
  void get_IDcoords(int, int, double *);
  double get_temperature(tagint **, int, int);
  void customvarnames();   // get per-atom variables names used by custom constraint
  void get_customvars();   // evaluate local values for variables names used by custom constraint
  double custom_constraint(std::string);   // evaulate expression for custom constraint
  double rxnfunction(std::string, std::string, std::string);   // eval rxn_sum and rxn_ave
  int get_chirality(double[12]);    // get handedness given an ordered set of coordinates

  void open(char *);
  void readline(char *);
  void parse_keyword(int, char *, char *);

  void far_partner();
  void close_partner();
  void get_molxspecials();
  void find_landlocked_atoms(int);
  void glove_ghostcheck();
  void ghost_glovecast();
  void update_everything();
  int insert_atoms(tagint **, int);
  void unlimit_bond();
  void limit_bond(int);
  void dedup_mega_gloves(int);    //dedup global mega_glove
  virtual void write_restart(FILE *);
  virtual void restart(char *buf);

  struct Set {
    int nreacts;
    char rxn_name[MAXLINE];
    int reaction_count_total;
  };
  Set *set;

  struct Constraint {
    int type;
    int id[MAXCONIDS];
    int idtype[MAXCONIDS];
    double par[MAXCONPAR];
    std::string str;
  };
  int ncustomvars;
  std::vector<std::string> customvarstrs;
  int nvvec;
  double **vvec;    // per-atom vector to store variable constraint atom-style variable values
  std::vector<std::vector<Constraint>> constraints;

  // DEBUG

  void print_bb();
};

}    // namespace LAMMPS_NS

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
        Bond/react: Atom type affected by reaction %s too close to template edge
        Bond/react: Bond type affected by reaction %s too close to template edge

This means an atom (or bond) that changes type or connectivity during the
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

E: Bond/react: Variable name does not exist

Self-explanatory.

E: Bond/react: Variable is not equal-style

Self-explanatory.

E: Bond/react: Molecule fragment does not exist

Self-explanatory.

*/
