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
   Contributing Author: Jacob Gissinger (jgissing@stevens.edu)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(bond/react,FixBondReact);
// clang-format on
#else

#ifndef LMP_FIX_BOND_REACT_H
#define LMP_FIX_BOND_REACT_H

#include "fix.h"

#include <array>
#include <deque>
#include <map>
#include <set>

namespace LAMMPS_NS {

class FixBondReact : public Fix {
 public:
  FixBondReact(class LAMMPS *, int, char **);
  ~FixBondReact() override;
  int setmask() override;
  void post_constructor() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void post_integrate() override;
  void post_integrate_respa(int, int) override;
  void post_force(int) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double compute_vector(int) override;
  double memory_usage() override;

 private:
  static constexpr double BIG = 1.0e20;
  static constexpr int MAXGUESS = 20;                      // max # of guesses allowed by superimpose algorithm
  static constexpr int MAXCONARGS = 14;                    // max # of arguments for any type of constraint + rxnID
  static constexpr int MAXLINE = 1024;                     // max length of line read from files
  static constexpr int MAXNAME = 256;                      // max character length of react-ID
  enum class Status { ACCEPT, REJECT, PROCEED,
                      CONTINUE, GUESSFAIL, RESTORE };      // values for superimpose algorithm status
  enum class Reset_Mol_IDs { YES, NO, MOLMAP };            // values for reset_mol_ids keyword
  enum class Molecule_Keys { OFF, INTER, INTRA };          // values for molecule_keyword
  enum class Dedup_Modes { LOCAL, GLOBAL };                // flag for one-proc vs shared reaction sites

  int newton_bond;
  FILE *fp;
  tagint lastcheck;
  int stabilization_flag;
  Reset_Mol_IDs molid_mode;
  int custom_exclude_flag;
  int rescale_charges_anyflag;                             // indicates if any reactions do charge rescaling
  int nrxnfunction;
  std::vector<std::string> rxnfunclist;                    // lists current special rxn function
  std::vector<int> peratomflag;                            // 1 if special rxn function uses per-atom variable (vs. per-bond)
  int atoms2bondflag;                                      // 1 if atoms2bond map has been populated on this timestep
  Status status;

  struct Reaction {
    int ID;
    class Molecule *reactant;                              // pre-reacted molecule template
    class Molecule *product;                               // post-reacted molecule template
    std::string name, constraintstr;
    std::string mapfilename;
    int nevery, groupbits;
    int iatomtype, jatomtype;
    int ibonding, jbonding;
    int closeneigh;                                        // indicates if bonding atoms of a rxn are 1-2, 1-3, or 1-4 neighbors
    double rminsq, rmaxsq;
    double fraction;
    double mol_total_charge;                               // sum of charges of post-reaction atoms whose charges are updated
    int reaction_count, reaction_count_total;
    int local_rxn_count, ghostly_rxn_count;
    int nlocalkeep, nghostlykeep;
    int seed, limit_duration;
    int stabilize_steps_flag;
    int custom_charges_fragid;
    int rescale_charges_flag;                              // if nonzero, indicates number of atoms whose charges are updated
    int create_atoms_flag, modify_create_fragid;
    double overlapsq;
    Molecule_Keys molecule_keyword;
    int v_nevery, v_rmin, v_rmax, v_prob;                  // ID of variable, -1 if static
    int nnewmolids;                                        // number of unique new molids needed for each reaction
    std::vector<std::array<tagint, 2>> attempts;           // stores sim atom IDs of initiator atoms

    struct ReactionAtomFlags {
      int edge;                                            // true if atom in molecule template has incorrect valency
      int landlocked;                                      // true if atom is at least three bonds away from edge atoms
      int recharged;                                       // true if atom whose charge should be updated
      int deleted;                                         // true if atom in pre-reacted template to delete
      int created;                                         // true if atom in post-reacted template to create
      int newmolid;                                        // for molmap option: mol IDs in post, but not in pre, re-indexed from 1
      std::array<int, 6> chiral;                           // pre-react chiral atoms. 1) flag 2) orientation 3-4) ordered atom types
      std::array<int, 2> amap;                             // atom map: clmn 1 = product atom IDs, clmn 2: reactant atom IDs
      std::array<int, 2> ramap;                            // reverse amap
    };
    std::vector<ReactionAtomFlags> atoms;

    struct Constraint {
      int ID;
      enum class Type { DISTANCE, ANGLE, DIHEDRAL, ARRHENIUS, RMSD, CUSTOM } type;
      struct Distance { double rminsq, rmaxsq; } distance;
      struct Angle { double amin, amax; } angle;
      struct Dihedral { double amin, amax, amin2, amax2; } dihedral;
      struct RMSD { double rmsdmax; } rmsd;
      struct Arrhenius { double A, n, E_a, seed; class RanMars *rrhandom; } arrhenius;
      struct Custom { std::string str; } custom;
      enum class IDType { ATOM, FRAG };
      static constexpr int MAXCONIDS = 4;
      std::array<int, MAXCONIDS> ids;
      std::array<IDType, MAXCONIDS> idtypes{};
      bool satisfied;
    };
    std::vector<Constraint> constraints;
  };
  std::vector<Reaction> rxns;

  int ncustomvars;
  std::vector<std::string> customvarstrs;
  int nvvec;
  double **vvec;                                           // per-atom vector to store custom constraint atom-style variable values
  class Compute *cperbond;                                 // pointer to 'compute bond/local' used by custom constraint ('rxnbond' function)
  std::map<std::set<tagint>, int> atoms2bond;              // maps atom pair to index of local bond array

  int nmax;                                                // max num local atoms
  int max_natoms;                                          // max natoms in a molecule template
  tagint *partner, *finalpartner;
  double **distsq;
  int allnattempt;
  unsigned shuffle_seed;                                   // user-provided value for the 'shuffle_seed' common keyword

  Fix *fix1;                                               // nve/limit used to relax reaction sites
  Fix *fix2;                                               // properties/atom used to indicate 1) relaxing atoms 2) to which 'react' atom belongs
  Fix *fix3;                                               // property/atom used for system-wide thermostat
  class RanMars **random;                                  // random number for 'prob' keyword
  class NeighList *list;
  class ResetAtomsMol *reset_mol_ids;                      // class for resetting mol IDs

  std::string nve_limit_xmax;                              // indicates max distance allowed to move when relaxing
  std::string id_fix1;                                     // id of internally created fix nve/limit
  std::string id_fix2;                                     // id of internally created fix per-atom properties
  std::string id_fix3;                                     // id of internally created 'stabilization group' per-atom property fix
  std::string statted_id;                                  // name of 'stabilization group' per-atom property
  std::string master_group;                                // group containing relaxing atoms from all fix rxns
  std::string exclude_group;                               // group for system-wide thermostat

  Reaction *rxnptr;                                        // for reverse_comm
  int countflag, commflag;
  int nlevels_respa;

  struct Superimpose {
    int avail_guesses;                                     // num of restore points available
    std::vector<int> guess_branch;                         // used when there is more than two choices when guessing
    struct StatePoint {
      int pion, neigh, trace, glove_counter;
      std::vector<tagint> glove, pioneer_count, pioneers;
    } sp;
  };
  std::vector<Superimpose::StatePoint> restore_pts;

  int **nxspecial;                                         // full number of 1-4 neighbors
  tagint **xspecial;                                       // full 1-4 neighbor list

  int cuff;                                                // extra space in mega_gloves: default = 1, w/ rescale_charges_flag = 2
  std::vector<std::vector<double>> my_mega_glove;          // local + ghostly reaction instances. for all mega_gloves: first row = rxnID.
  std::vector<std::vector<double>> local_mega_glove;       // consolidation of local reaction instances
  std::vector<std::vector<double>> ghostly_mega_glove;     // consolidation of nonlocal reaction instances
  double **global_mega_glove;                              // consolidation (inter-processor) of gloves containing nonlocal atoms

  int *localsendlist;                                      // indicates ghosts of other procs
  int my_num_mega;                                         // local + ghostly reaction instances (on this proc)
  int local_num_mega;                                      // num of local reaction instances
  int ghostly_num_mega;                                    // num of ghostly reaction instances
  int global_megasize;                                     // num of reaction instances in global_mega_glove

  void validate_variable_keyword(const char *, int);
  void read_map_file(Reaction &);
  void EdgeIDs(char *, Reaction &, int);
  void Equivalences(char *, Reaction &, int);
  void DeleteAtoms(char *, Reaction &, int);
  void CreateAtoms(char *, Reaction &, int);
  void CustomCharges(int, Reaction &);
  void ChiralCenters(char *, Reaction &, int);
  void ReadConstraints(char *, Reaction &);
  void readID(char *, Reaction::Constraint &, Reaction &, int);

  void superimpose_algorithm();
  void make_a_guess(Superimpose &, Reaction &);
  void neighbor_loop(Superimpose &, Reaction &);
  void check_a_neighbor(Superimpose &, Reaction &);
  void crosscheck_the_neighbor(Superimpose &, Reaction &);
  void inner_crosscheck_loop(Superimpose &, Reaction &);
  int ring_check(Reaction &, std::vector<tagint> &);
  int check_constraints(Reaction &, std::vector<tagint> &);
  void get_IDcoords(Reaction::Constraint::IDType, int, double *, Molecule *, std::vector<tagint> &);
  double get_temperature(std::vector<tagint> &);
  double get_totalcharge(Reaction &, std::vector<tagint> &);
  void customvarnames();                                   // get per-atom variables names used by custom constraint
  void get_customvars();                                   // evaluate local values for variables names used by custom constraint
  bool custom_constraint(const std::string &, Reaction &, std::vector<tagint> &);
  double rxnfunction(const std::string &, const std::string &, const std::string &, Molecule *, std::vector<tagint> &);
  void get_atoms2bond(int);
  int get_chirality(double[12]);                           // get handedness given an ordered set of coordinates

  void readline(char *);
  void parse_keyword(int, char *, char *);

  void far_partner(Reaction &);
  void close_partner(Reaction &);
  void find_landlocked_atoms(Reaction &);
  void glove_ghostcheck();
  void ghost_glovecast();
  void update_everything();
  int insert_atoms_setup(tagint **, int, Reaction &);
  void unlimit_bond();                                     // removes atoms from stabilization, and other post-reaction every-step operations
  void dedup_mega_gloves(Dedup_Modes);                     // dedup global mega_glove
  void write_restart(FILE *) override;
  void restart(char *buf) override;

  // store restart data
  struct Set {
    int nrxns;
    char rxn_name[MAXNAME];
    int reaction_count_total;
    int nratelimits;
  };
  Set *set;

  tagint addatomtag;
  struct AddAtom {
    tagint tag, molecule;
    int type, mask;
    imageint image;
    double rmass, x[3], v[3];
  };
  std::vector<AddAtom> addatoms;

  struct RateLimit {
    int Nrxns, var_flag, var_id, Nlimit, Nsteps;
    std::vector<int> rxnIDs;
    std::vector<std::string> rxn_names;
    std::deque<int> store_rxn_counts;
  };
  std::vector<RateLimit> rate_limits;

  struct MaxRxnLimit {
    int Nrxns, max_rxn;
    std::vector<int> rxnIDs;
    std::vector<std::string> rxn_names;
  };
  std::vector<MaxRxnLimit> max_rxn_limits;

  // DEBUG

  void print_bb();
};

}    // namespace LAMMPS_NS

#endif
#endif
