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

#ifdef FIX_CLASS
// clang-format off
FixStyle(msevb,FixMSEVB);
// clang-format on
#else

#ifndef LMP_FIX_MSEVB_H
#define LMP_FIX_MSEVB_H

#include "fix.h"

#include <string>
#include <vector>

namespace LAMMPS_NS {

class FixMSEVB : public Fix {
 public:
  FixMSEVB(class LAMMPS *, int, char **);
  ~FixMSEVB() override;

  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void setup(int) override;
  void min_setup(int) override;
  void post_integrate() override;
  void setup_pre_force(int) override;
  void pre_force(int) override;
  void post_force(int) override;
  void min_pre_force(int) override;
  void min_post_force(int) override;
  void post_run() override;
  double compute_scalar() override;
  double compute_vector(int) override;
  void *extract(const char *, int &) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 protected:
  // ---- Coupling type enum -----------------------------------------------
  enum CouplingType {
    COUPLING_NONE,
    COUPLING_RAITERI2011,
    COUPLING_VUILLEUMIER1998,
    COUPLING_GRIMME2015
  };

  // ---- Reaction definition (per-reaction) ------------------------------
  // All reactions are now template-based: pre/post molecule files + map file.
  struct ReactionDef {
    int type_H, type_Y;    // atom types (from pre_mol), used for fast scan
    int type_X;            // type of X (bond-breaking partner of H); -1 if none
    double cutoff_sq;

    // Per-reaction coupling (optional override; defaults from global coupling)
    bool coupling_set;
    CouplingType coupling_type;
    double coupling_lambda, coupling_zeta;
    double coupling_v12, coupling_alpha, coupling_gamma_v;
    double coupling_a, coupling_b;
    double coupling_taper;

    // Per-reaction shell depth; -1 = inherit global default
    int shells;

    // Topology diff sub-structs
    struct TypeChange {
      int pre_idx;
      int new_type;
    };
    struct BondBreak {
      int pre_idx1, pre_idx2;
    };
    struct BondCreate {
      int pre_idx1, pre_idx2;
      int bond_type;
    };
    struct BondRetype {
      int pre_idx1, pre_idx2;
      int new_bond_type;
    };    // same pair, type changes

    // Names / paths used for init() loading
    std::string pre_mol_id, post_mol_id, map_file_path;

    // Borrowed pointers (owned by lmp->atom->molecules[])
    class Molecule *pre_mol;
    class Molecule *post_mol;

    int glove_n;       // mol->natoms for pre-template (0 until init())
    int ibonding;      // 0-based pre-template index of H
    int jbonding;      // 0-based pre-template index of Y
    int ix_bonding;    // 0-based pre-template index of X (-1 if no breaking bond)

    std::vector<int> is_edge;        // [glove_n]: 1=edge atom (not matched by BFS)
    std::vector<int> pre_to_post;    // [glove_n]: pre_to_post[i]=0-based post idx (-1=none)

    std::vector<TypeChange> type_changes;
    std::vector<BondBreak> bond_breaks;
    std::vector<BondCreate> bond_creates;
    std::vector<BondRetype> bond_retypes;    // bonds that persist but change type

    ReactionDef() :
        type_H(0), type_Y(0), type_X(-1), cutoff_sq(0.0), coupling_set(false),
        coupling_type(COUPLING_NONE), coupling_lambda(0.0), coupling_zeta(0.0), coupling_v12(0.0),
        coupling_alpha(0.0), coupling_gamma_v(0.0), coupling_a(0.0), coupling_b(0.0),
        coupling_taper(0.0), shells(-1), pre_mol(nullptr), post_mol(nullptr),
        glove_n(0), ibonding(-1), jbonding(-1), ix_bonding(-1)
    {
    }
  };

  // ---- Reactive site --------------------------------------------------
  static const int MAX_COMPONENTS = 8;
  struct ReactiveSite {
    tagint tag_H, tag_X, tag_Y;
    int idx_H, idx_X, idx_Y;
    double dist_sq;
    int rxn_idx;
    int parent_state;
    int shell;
    int chain_len;
    int n_components;                  // 0 = single-site; >0 = product state
    int components[MAX_COMPONENTS];    // indices into sites[] for each component
    // For a product state, neighbor_state[j] is the STATE index (site+1, or 0
    // for the reference) of the configuration obtained by removing component j
    // -- i.e. the state that differs from this one by exactly component j's
    // transfer.  Used to place the off-diagonal coupling of an n-way product
    // onto its (n-1)-component neighbours.  Recomputed locally each step.
    int neighbor_state[MAX_COMPONENTS];
  };

  // ---- Partition / communicator ---------------------------------------
  int npartitions;
  int ipartition;
  MPI_Comm samerank;

  // ---- Communication buffer -------------------------------------------
  double *commbuf;
  int nmax;

  // ---- Computes -------------------------------------------------------
  class Compute *pe;
  std::string id_pe;
  class Compute *temp_compute;
  class Compute *press_compute;
  std::string id_temp, id_press;
  double pressure[6];

  // ---- Per-partition potential energies -------------------------------
  double *epot;
  double *my_epot;
  MPI_Request epot_request;
  int epot_nmax;
  double epot_ground;

  // ---- Eigensystem arrays (nstates x nstates, grown dynamically) ------
  double *hamiltonian;
  double *H_work;
  double *eigenvalues;
  double *eigenvectors;
  double *amplitudes;
  double *fd_occ;    // FD occupancy per eigenstate; null when FD disabled
  int eigensys_nmax;

  // ---- Consistency check scratch arrays (npartitions, fixed) ----------
  int *my_nlocal;
  int *all_nlocal;

  // ---- Weight-based force mixing --------------------------------------
  double *weights;    // [nstates] scalar mixing weights per EVB state
  int weights_nmax;

  // Sparse coupling correction for position-based (Raiteri/Vuilleumier) edges
  struct CouplingCorrection {
    tagint tag_H, tag_X, tag_Y;
    double fH[3], fX[3], fY[3];
    double rho_factor;    // 2 * rho[parent][child]
  };
  std::vector<CouplingCorrection> coupling_corrections;
  double coupling_virial[6];    // virial from coupling correction forces

  // Result of evaluating coupling for a single edge
  struct CouplingResult {
    double V;                      // coupling value (Hamiltonian off-diagonal)
    double g;                      // Grimme weight scalar (-2*b*dE*V); 0 for position-based
    double fH[3], fX[3], fY[3];    // sparse forces; valid when has_forces=true
    bool valid;                    // false if needed atoms not found on this rank
    bool is_grimme;                // true for COUPLING_GRIMME2015
    bool has_forces;               // true if fH/fX/fY are set
  };

  CouplingResult evaluate_coupling(const ReactionDef &rxn, double E_parent, double E_child,
                                   tagint tag_H, tagint tag_X, tagint tag_Y);

  // ---- Excess state force storage (Approach A: save/restore) ----------
  double *saved_forces;             // [nmax*3] save/restore atom->f around excess eval
  double *excess_forces;            // [nsites_serial * ef_nmax * 3] per-excess-state forces
  int excess_forces_nmax;           // max atoms dimension of excess_forces
  int excess_forces_nmax_serial;    // max nsites_serial dimension

  // ---- Excess force mode: 0=save/restore (A), 1=split (B) -----------
  int excess_force_mode;

  // ---- Reactive sites -------------------------------------------------
  std::vector<ReactionDef> rxndefs;
  int reaction_enabled;
  int enumerate_product_states;
  ReactiveSite *sites;
  int sites_nmax;
  int nsites;
  int nsites_parallel;
  int nsites_serial;
  int nstates;
  int max_shells;

  // Per-site chain data: index = site_idx * max_shells + d
  tagint *chain_H_flat;
  tagint *chain_X_flat;
  tagint *chain_Y_flat;
  int *chain_rxn_flat;    // per-depth rxn_idx: chain_rxn_flat[site * max_shells + d]
  tagint *
      chain_glove_flat;    // per-depth glove: chain_glove_flat[(site * max_shells + d) * glove_nmax + g]

  // Per-site glove for template reactions: glove_flat[site_idx * glove_nmax +
  // i] = real global tag for pre-template atom i. Zero for unmatched atoms.
  tagint *glove_flat;
  int glove_nmax;    // max pre-template natoms across all reactions; 0 if none

  // ---- Reference topology snapshot ------------------------------------
  int *ref_type;
  double *ref_charge;
  tagint *ref_molecule;
  int *ref_num_bond;
  int *ref_bond_type_flat;
  tagint *ref_bond_atom_flat;
  int *ref_num_angle;
  int *ref_angle_type_flat;
  tagint *ref_angle_atom1_flat;
  tagint *ref_angle_atom2_flat;
  tagint *ref_angle_atom3_flat;
  int *ref_num_dihedral;
  int *ref_dihedral_type_flat;
  tagint *ref_dihedral_atom1_flat;
  tagint *ref_dihedral_atom2_flat;
  tagint *ref_dihedral_atom3_flat;
  tagint *ref_dihedral_atom4_flat;
  int *ref_num_improper;
  int *ref_improper_type_flat;
  tagint *ref_improper_atom1_flat;
  tagint *ref_improper_atom2_flat;
  tagint *ref_improper_atom3_flat;
  tagint *ref_improper_atom4_flat;
  int *ref_nspecial_flat;
  tagint *ref_special_flat;
  int ref_maxspecial;
  tagint ref_maxtag;
  int ref_bond_per_atom;
  int ref_angle_per_atom;
  int ref_dihedral_per_atom;
  int ref_improper_per_atom;
  int ref_snapshot_valid;

  // ---- Neighbor list --------------------------------------------------
  class NeighList *list;

  // ---- Coupling -------------------------------------------------------
  CouplingType coupling_type;
  double coupling_lambda, coupling_zeta;
  double coupling_v12, coupling_alpha, coupling_gamma_v;
  double coupling_a, coupling_b;
  double coupling_taper;
  int coupling_enabled;

  // ---- Fermi-Dirac ----------------------------------------------------
  int fermi_dirac_enabled;
  double fd_temperature;
  double fd_RT;

  // ---- Reactive-atom group tracking -----------------------------------
  // A LAMMPS group named "<fix_id>_atoms" is created in init() and kept
  // current after every detect_reactive_sites() call.  It contains every
  // atom touched by any active reactive site (H, X, Y for every chain
  // depth of every non-product-state site).  Use this group in dump
  // commands to write only MSEVB-relevant atoms.
  std::string reactive_group_name;    // "<id>_atoms"
  int reactive_group_bit;             // bitmask bit for that group; 0 = not yet created

  void update_reactive_group();    // called after detect_reactive_sites()

  // ---- SCF topology minimisation --------------------------------------
  bool scf_topology;    // default: true — re-evaluate EVB until no more transfers
  int scf_max_iter;     // max SCF iterations per step (default: 10)

  // ---- Minimization termination handshake ----------------------------
  bool min_terminate;    // set on non-zero partitions when P0 sends stop signal

  // ---- Misc -----------------------------------------------------------
  bool partition_warning{false};
  bool product_states_multishell_warning{false};

  // ---- Structured (JSON) output ---------------------------------------
  std::string file_name;    // path given to the 'file' keyword
  bool file_flag;           // true if the 'file' keyword was given
  int file_every;           // >0: also write every N steps; 0: on reaction only
  FILE *fpout;              // open only on the universe root (ipartition 0, me 0)
  int json_init;            // 0 until the first record is written (comma control)

  void write_msevb_json(bigint timestep, int dominant_state, double dominant_amp);

  // ---- Helper: dynamic array growth -----------------------------------
  void grow_eigensystem_arrays();
  void grow_excess_forces(int nlocal_max);
  void grow_epot_arrays();
  void grow_sites(int n);

  // ---- Helper: reference topology memory management ------------------
  void allocate_ref_topology(tagint sz, int bpa, int apa, int dpa, int ipa, int mspc);
  void destroy_ref_topology();
  void broadcast_reference_topology();

  // ---- Helper: apply state changes for one site (all chain depths) ---
  // Calls apply_state_change for each depth of site sk's transfer chain,
  // with forward_comm(this) between consecutive steps within a component
  // and between components of a product state.
  void apply_site_changes(int sk);

  void init_site(int site_idx, tagint tag_H, tagint tag_X, tagint tag_Y, double dist_sq);

  void compute_excess_states();      // Approach A: save/restore + buffer
  void compute_excess_energies();    // Approach B: energy-only pass
  void apply_excess_forces();        // Approach B: deferred force application

  // ---- Methods --------------------------------------------------------
  void check_consistency_atoms();
  void build_hamiltonian();
  void solve_eigensystem();
  void compute_coupling_values();
  void compute_mixing_weights();
  void apply_coupling_corrections();
  void weight_based_hellmann_feynman_forces();
  void fermi_dirac_occupancies(double *evals, int ns, double *occupancies);

  void sync_positions();
  void apply_per_partition_state_changes();
  void gather_potential_energies();
  bool do_permanent_transfer(int &out_max_state, double &out_max_amp);
  void snapshot_reference_topology();
  void restore_reference_topology();
  void detect_reactive_sites();

  // Apply one template state change using the given glove mapping.
  void apply_state_change(const tagint *glove, int idx_X, int idx_H, int idx_Y, tagint tag_X,
                          tagint tag_H, tagint tag_Y, const ReactionDef &rxn);

  // Template-reaction helpers called from init().
  void parse_template_mapfile(ReactionDef &rd);
  void compute_topology_diff(ReactionDef &rd);
  void check_reaction_topology_limits(const ReactionDef &rd);

  // ---- Reactive-site atom-set primitives ------------------------------
  // The set of real atom tags that a site's transfer chain *modifies*: the
  // non-edge glove atoms across every chain depth (edge atoms are matched but
  // left untouched).  This is the atom_set() of the Chain abstraction — the
  // basis for deciding whether two chains can coexist in one (product) state.
  std::vector<tagint> chain_atom_set(int site) const;
  // True iff two sites' chains modify disjoint atom sets (i.e. they may be
  // combined into a product state without clobbering each other's atoms).
  bool chains_disjoint(int site_a, int site_b) const;
  // After sites[] (incl. product states) is finalized and broadcast, fill each
  // product state's neighbor_state[] by mapping its (n-1)-component subsets to
  // their state indices, so the coupling routines can reach them.
  void build_product_neighbors();
  // Collective shell depth (sum of component chain lengths) and aggregate
  // transfer distance of a site -- the (primary, secondary) ranking keys used
  // to prune to max_states.
  int site_collective_shells(int site) const;
  double site_aggregate_dist(int site) const;
  // Prune sites[] (on partition 0, before broadcast) to keep at most
  // max_states-1 reactive states, ranked by (collective shells, aggregate
  // distance) ascending, and re-index the survivors.  The ranking is closed
  // under component/parent removal, so kept products keep their neighbours.
  void prune_states_to_max();

  // Maximum total EVB states (reference + reactive) evaluated per step.
  // 0 = unlimited.  When exceeded, states are pruned keeping the lowest
  // collective shell depth first, then the smallest aggregate distance.
  int max_states;

  // ---- Kokkos sync hooks (no-ops in base class, overridden in /kk) ------
  virtual void sync_before_force_compute();
  virtual void modified_after_force_compute();
  virtual void sync_forces_to_host();
  virtual void modified_forces_on_host();
  virtual void modified_topology_on_host();
  virtual void sync_before_neighbor_build();
};

}    // namespace LAMMPS_NS

#endif
#endif
