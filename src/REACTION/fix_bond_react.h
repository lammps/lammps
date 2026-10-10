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
#include "reaction.h"
#include "topology_matcher.h"

#include <array>
#include <deque>
#include <map>
#include <memory>
#include <set>

namespace LAMMPS_NS {

struct json_metadata;                                      // forward declaration. full declaration in json_metadata.h

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
  std::string get_thermo_colname(int) override;
  int modify_param(int, char **) override;
  double memory_usage() override;

  void write_restart(FILE *) override;
  void restart(char *buf) override;

  int image(int *&, double **&) override;

 private:
  static constexpr double BIG = 1.0e20;
  enum class Reset_Mol_IDs { YES, NO, MOLMAP };            // values for reset_mol_ids keyword
  enum class Dedup_Modes { LOCAL, GLOBAL };                // flag for one-proc vs shared reaction sites

  int newton_bond;
  tagint lastcheck;
  FILE *fpout;
  bool outflag;
  int json_init;
  std::unique_ptr<json_metadata> rxn_metadata;
  int stabilization_flag;
  Reset_Mol_IDs molid_mode;
  int custom_exclude_flag;
  int rescale_charges_anyflag;                             // indicates if any reactions do charge rescaling

  std::vector<Reaction> rxns;
  TopologyMatcher *topo_matcher;

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

  // arrays for dump image rendering

  int *imgobjs;
  double **imgparms;
  std::map<tagint, int> vizatoms;  // maps atom IDs to number of steps they have been highlighted
  int vizsteps;                    // number of steps to highlight atoms in reactions

  void superimpose_algorithm();
  bool compare_atomtype(int, Reaction &, int);
  double get_totalcharge(Reaction &, std::vector<tagint> &);

  void far_partner(Reaction &);
  void close_partner(Reaction &);
  void find_landlocked_atoms(Reaction &);
  void glove_ghostcheck();
  void ghost_glovecast();
  void update_everything();
  int insert_atoms_setup(tagint **, int, Reaction &);
  void unlimit_bond();                                     // removes atoms from stabilization, and other post-reaction every-step operations
  void dedup_mega_gloves(Dedup_Modes);                     // dedup global mega_glove

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
    int Nrxns = 0, var_flag = 0, var_id = -1, Nlimit = 0, Nsteps = 0;
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
