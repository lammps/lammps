/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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

#include <map>
#include <set>

namespace LAMMPS_NS {

class FixBondReact : public Fix {
 public:
  enum { MAXLINE = 256 };    // max length of line read from files
  enum { MAXCONIDS = 4 };    // max # of IDs used by any constraint
  enum { MAXCONPAR = 5 };    // max # of constraint parameters

  FixBondReact(class LAMMPS *, int, char **);
  ~FixBondReact() override;
  int setmask() override;
  void post_constructor() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void post_integrate() override;
  void post_integrate_respa(int, int) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double compute_vector(int) override;
  double memory_usage() override;

 private:
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
  int **rate_limit;
  int **store_rxn_count;
  int *stabilize_steps_flag;
  int *custom_charges_fragid;
  int *rescale_charges_flag;
  int *create_atoms_flag;
  int *modify_create_fragid;
  double *overlapsq;
  int *molecule_keyword;
  int maxnconstraints;
  int *nconstraints;
  char **constraintstr;
  int nrxnfunction;
  std::vector<std::string> rxnfunclist;     // lists current special rxn function
  std::vector<int> peratomflag; // 1 if special rxn function uses per-atom variable (vs. per-bond)
  int atoms2bondflag;           // 1 if atoms2bond map has been populated on this timestep
  int narrhenius;
  int **var_flag, **var_id;     // for keyword values with variable inputs
  int status;
  int *groupbits;

  int rxnID;          // integer ID for identifying current bond/react
  char **rxn_name;    // name of reaction
  int *reaction_count;
  int *reaction_count_total;
  int nmax;          // max num local atoms
  int max_natoms;    // max natoms in a molecule template
  int max_rate_limit_steps;    // max rate limit interval
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
  class ResetAtomsMol *reset_mol_ids;    // class for resetting mol IDs

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
  tagint **global_mega_glove;     // consolidation (inter-processor) of gloves
                                  // containing nonlocal atoms

  int *localsendlist;      // indicates ghosts of other procs
  int local_num_mega;      // num of local reaction instances
  int ghostly_num_mega;    // num of ghostly reaction instances
  int global_megasize;     // num of reaction instances in global_mega_glove
  int *pioneers;           // during Superimpose Algorithm, atoms which have been assigned,
                           // but whose first neighbors haven't
  int glove_counter;       // used to determine when to terminate Superimpose Algorithm

  void read_variable_keyword(const char *, int, int);
  void read_map_file(int);
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
  void customvarnames();    // get per-atom variables names used by custom constraint
  void get_customvars();    // evaluate local values for variables names used by custom constraint
  double custom_constraint(const std::string &);    // evaulate expression for custom constraint
  double rxnfunction(const std::string &, const std::string &,
                     const std::string &);    // eval rxn_sum and rxn_ave
  void get_atoms2bond(int);
  int get_chirality(double[12]);              // get handedness given an ordered set of coordinates

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
  void unlimit_bond(); // removes atoms from stabilization, and other post-reaction every-step operations
  void dedup_mega_gloves(int);    //dedup global mega_glove
  void write_restart(FILE *) override;
  void restart(char *buf) override;

  // store restart data
  struct Set {
    int nreacts;
    char rxn_name[MAXLINE];
    int reaction_count_total;
    int max_rate_limit_steps;
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
  double **vvec;    // per-atom vector to store custom constraint atom-style variable values
  class Compute *cperbond;    // pointer to 'compute bond/local' used by custom constraint ('rxnbond' function)
  std::map<std::set<tagint>, int> atoms2bond;    // maps atom pair to index of local bond array
  std::vector<std::vector<Constraint>> constraints;

  // DEBUG

  void print_bb();
};

}    // namespace LAMMPS_NS

#endif
#endif
