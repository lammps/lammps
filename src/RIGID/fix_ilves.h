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
FixStyle(ilves,FixIlves);
// clang-format on
#else

#ifndef LMP_FIX_ILVES_H
#define LMP_FIX_ILVES_H

#include "fix.h"

#include <unordered_map>
#include <vector>

namespace LAMMPS_NS {

class FixRespa;

// ILVES bond/angle constraint solver.  Gathers the complete bond/angle
// topology onto every MPI rank once at init via MPI_Allgatherv, then
// builds a per-rank constraint list that includes every constraint
// cluster intersecting at least one local atom -- even constraints
// between two ghost atoms.  Supports arbitrary cluster topologies and
// any spatial extent (clusters spanning many subdomains, e.g. an
// all-backbone-constrained polymer).  The replicated topology is the
// memory price for this generality.

class FixIlves : public Fix {

 public:
  FixIlves(class LAMMPS *, int, char **);
  ~FixIlves() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void pre_neighbor() override;
  void post_neighbor() override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void end_of_step() override;
  void min_pre_reverse(int, int) override;
  void min_post_force(int) override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;

  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

  bigint dof(int) override;
  void reset_dt() override;
  double compute_scalar() override;

 protected:
  // command-line settings
  double tolerance;      // convergence tolerance
  int max_iter;          // max Newton iterations
  int output_every;      // statistics print interval (0 = never)
  int variant;           // 0 = full (asymmetric Newton, LU)
                         // 1 = fast (symmetric quasi-Newton, banded Cholesky)
  bigint next_output;

  // Near-linear angle handling.  When an angle's equilibrium theta_0 is
  // close to 180 degrees, constraining the three legs {AB, BC, AC} of
  // the angle's atom triple makes the constraint Jacobian rank-deficient
  // (the triangle inequality saturates: |AC| = |AB|+|BC|, so the three
  // constraints become linearly dependent at exactly 180 deg, and very
  // ill-conditioned nearby).  For these angles we drop the bond between
  // B and the higher-tag endpoint of {A, C}, and instead add a 3-atom
  // "virtual" constraint |B - (A+C)/2|^2 = L^2 with L the median length
  // sqrt((|AB|^2 + |BC|^2)/2 - |AC|^2/4).  The retained set
  // {AB-or-BC, AC, B-M} keeps |B-C| within solver tolerance via the
  // triangle geometry while remaining well-conditioned at near-180 deg.
  // Controlled by the 'linearangle <deg>' keyword; default 165 deg.
  // Set to >= 180 to disable.
  double linear_threshold;        // degrees; angle types with theta_0 >=
                                  // this value use the alternate set
  double linear_Lmin;             // minimum |B-M| target; the natural
                                  // median length is clamped up to this
                                  // floor when the geometry would give
                                  // a near-zero target (symmetric ~180
                                  // case).  Default: 0.01 (length units).
  int *angle_linear;              // angle_linear[at] = 1 if angle type at
                                  // is near-linear (size nangletypes+1)
  double *angle_dBM;              // target |B-M| length for near-linear
                                  // angle type at (size nangletypes+1)

  // selection flags (allocated 1..n)
  int *bond_flag;     // bond_flag[bt]   = 1 if bond  type bt is constrained
  int *angle_flag;    // angle_flag[at]  = 1 if angle type at is constrained
  int *type_flag;     // type_flag[at]   = 1 if any atom of type at triggers
  double *mass_list;  // list of atom masses (with tolerance) that trigger
  int nmass;

  // equilibrium distances (1..n)
  double *bond_distance;     // bond_distance[bt]   = r0(bt)
  double *angle_distance;    // angle_distance[at]  = end-to-end distance for
                             //   angle type at (computed from the two flanking
                             //   bond types + angle, like fix shake)
  double *angle_r1;          // flanking bond equilibrium lengths for angle type at
  double *angle_r2;          //   (r_AB and r_BC); stored so stats() can invert d_AC -> angle
  bool has_angle;            // true if any angle type is selected

  // restraint mode (minimization)
  double kbond;     // spring constant for harmonic restraint substitute
  double ebond;     // accumulated restraint energy

  // store constraint forces per atom (fix_modify store yes)
  int store_flag;
  int maxstore;
  double **fstore;

  // per-atom data (migrates with atoms)
  int *ilves_flag;    // 1 if atom participates in any constrained bond/angle
  double **xshake;    // working buffer for unconstrained-then-projected coords

  // flat constraint list (rebuilt each reneighbor)
  // c_atom1[k] is normally a local atom (< nlocal) for ownership uniqueness;
  // it may be a ghost for "cluster completion" entries where every atom of
  // the constraint is a ghost on this rank.
  //
  // c_atom3[k] = -1 for ordinary 2-atom constraints.  For "B-M virtual"
  // constraints (added for near-linear angles, see linear_threshold), it
  // holds the index of the third atom -- by convention c_atom1 = B (the
  // central atom, gradient weight +1), and c_atom2, c_atom3 are A and C
  // (the angle endpoints, gradient weight -1/2 each).  The constraint
  // vector r_k is then x[B] - (x[A]+x[C])/2 (closest image) with
  // |r_k|^2 = c_dist[k]^2 as the target.
  int n_constr;
  int max_constr;
  int *c_atom1;
  int *c_atom2;
  int *c_atom3;      // -1 for 2-atom constraints; third atom index for B-M
  int *c_type;       // bond type for bond entries, or -atype for virtual
                     // entries from an angle (A-C or B-M)
  double *c_dist;    // target distance (not squared) for clarity
  double *c_lambda;  // accumulated Lagrange multiplier (Phase 2+)

  // connected-component label for each constraint
  int *c_cluster;
  int n_clusters;

  // precomputed per-constraint reference quantities (rebuilt each reneighbor)
  // For 2-atom constraints (c_atom3 == -1): r_k = x[a1]-x[a2], c_invma/c_invmb
  //   are 1/mass[a1] and 1/mass[a2], and c_invmc is unused.
  // For 3-atom B-M constraints (c_atom3 >= 0): r_k = x[a1] - (x[a2]+x[a3])/2,
  //   c_invma = 1/mass[a1] (the central atom B), c_invmb = 1/mass[a2] and
  //   c_invmc = 1/mass[a3] (the endpoints A and C).
  double *c_rx, *c_ry, *c_rz;
  double *c_rsq;
  double *c_invma, *c_invmb, *c_invmc;

  // cluster grouping (built alongside connected components in build_constraint_list)
  int *cluster_offset;           // size n_clusters+1, CSR-style ranges into c_perm
  int *c_perm;                   // size n_constr; c_perm[s] = global constraint id at slot s
  int *c_slot;                   // size n_constr; inverse map: c_slot[k] = slot of constraint k

  // Newton workspace (sized to largest cluster)
  int lu_alloc;
  double *lu_A;                  // row-major (lu_alloc * lu_alloc)
  double *lu_b;                  // rhs / dlambda
  int *lu_pivot;
  double *cl_sx, *cl_sy, *cl_sz; // s_k cache for the current cluster
  int largest_cluster;           // size of the largest cluster (informational)

  // Per-cluster Cholesky factor cache (variant "fast" only).  The cluster
  // Jacobian in the "fast" variant is built entirely from r_k, c_rsq, and
  // inverse masses -- all step-constant within solve_constraints() -- so
  // its Cholesky factor can be computed once at the first Newton iteration
  // and reused (forward+back substitution only) on subsequent iterations.
  //
  // Storage is banded: cluster_bw[c] is the bandwidth of cluster c's
  // matrix after the RCM reordering applied at constraint-list build time
  // (smaller is faster).  The lower-band-packed L for cluster c lives at
  // chol_pool[chol_pool_offset[c] ..] with n_c rows of (bw_c + 1) entries
  // each, row-major.  L[i][j] (j <= i) is at offset i*(bw_c+1) + (i-j);
  // the diagonal lives at offset 0 of each row.  cluster_cached[c] is 1
  // when slot c holds a valid factor for the current step (cleared at
  // start of every solve_constraints and on Cholesky fallback to LU).
  double *chol_pool;
  bigint chol_pool_alloc;        // current allocation in doubles
  bigint *chol_pool_offset;      // size n_clusters+1
  int chol_offset_alloc;
  int *cluster_bw;               // bandwidth per cluster (after RCM)
  int cluster_bw_alloc;
  char *cluster_cached;          // size n_clusters
  int cluster_cached_alloc;
  int largest_bw;                // max bw across all clusters

  // pointers to atom-class quantities, refreshed at pre_neighbor / post_force
  double **x, **v, **f;
  double *mass, *rmass;
  int *type;
  int nlocal;

  // step constants
  double dtv, dtfsq;

  // rRESPA bookkeeping (mirrors fix shake).  respa=0 selects the Verlet
  // path (dtfsq = dt^2*ftm2v); respa=1 selects the multi-level path where
  // dtfsq is recomputed per level from step_respa[ilevel] and the saved
  // per-level forces in fix_respa->f_level are folded into xshake.
  int respa;
  int nlevels_respa;
  int *loop_respa;
  double *step_respa;
  class FixRespa *fix_respa;
  double dtf_inner;

  // statistics, sized [nbondtypes+1] / [nangletypes+1]
  bigint *b_count, *b_count_all;
  double *b_ave, *b_max, *b_min;
  double *b_ave_all, *b_max_all, *b_min_all;
  bigint *a_count, *a_count_all;
  double *a_ave, *a_max, *a_min;
  double *a_ave_all, *a_max_all, *a_min_all;

  // remember vflag from post_force for end-of-step bookkeeping
  int vflag_post_force;
  int eflag_pre_reverse;

  // selector for pack_forward_comm: 0 = pack xshake, 1 = pack atom->v
  int comm_mode;

  // counters incremented per cluster-solve; reported on rank 0 when
  // output_every > 0 so the user sees how often we fall back to LU.
  bigint chol_calls, chol_fallbacks;

  // -----------------------------------------------------------------
  // Replicated global bond / angle topology (gathered once at init via
  // MPI_Allgatherv).  Bonds: (lower-tag, higher-tag, type), sorted by
  // (a, b) and deduped.  Angles: (atom1, mid, atom3, type), sorted by
  // middle then outer atoms; only angle types listed in angle_flag[]
  // are gathered.
  std::vector<tagint> gb_a, gb_b;
  std::vector<int> gb_type;
  std::vector<tagint> ga1, ga2, ga3;
  std::vector<int> ga_type;
  bool global_topology_ready;

  // Topology-change safeguard.  Fix ilves gathers the bond/angle topology
  // once at init() and assumes it is frozen for the duration of the run.
  // To catch silent mismatch when a concurrent fix (e.g. fix bond/create,
  // fix bond/break, delete_atoms) modifies the constrained topology, we
  // record three globally-reduced counts at init and recheck them at every
  // reneighbor; mismatch aborts with a clear error.  The bond/angle counts
  // are raw "sum over local atoms of constrained-type slots" -- not
  // deduplicated -- which is enough for change detection and matches
  // what gather_global_topology() scans.
  bigint natoms_at_init;
  bigint nconstrbonds_sum_at_init;
  bigint nconstrangles_sum_at_init;

  // Global cluster-ID for every tag involved in any bond/angle.  Maps
  // a tag to its cluster's representative tag.  Sparse: only tags that
  // participate in at least one constrained bond or angle appear in the
  // map.  Per-rank size is ~24-40 bytes per involved tag (unordered_map
  // overhead), keeping memory tractable for partial-constraint systems.
  std::unordered_map<tagint, tagint> tag_cluster;

  // -----------------------------------------------------------------
  // build helpers
  void grow_constraint_list(int);
  void add_constraint(int a, int b, int btype, double dist);
  // 3-atom variant: a = B (central), b = A, c = C; constraint vector is
  // x[B] - (x[A]+x[C])/2; btype is conventionally -atype (same as the
  // angle-virtual A-C entry it accompanies).
  void add_constraint(int a, int b, int c, int btype, double dist);
  void group_by_cluster();
  void precompute_constraint_data();
  int masscheck(double massone);
  bool bond_selected_for_atoms(int ia, int ib, int bt);

  void negate_constrained_topology();

  // global topology setup
  void init_topology();
  void gather_global_topology();
  void build_constraint_list();
  bool bond_is_constrained(tagint ta, tagint tb);

  // Topology-change detection.  count_constrained_topology() walks local
  // atoms once and MPI_Allreduces the constrained bond / angle slot counts
  // (plus reports atom->natoms).  record_topology_baseline() saves the
  // result into the *_at_init members.  check_topology_unchanged() compares
  // a fresh recount against the baseline and aborts with error->all() on
  // mismatch; it is called from pre_neighbor() so any change made by a
  // sibling fix between two steps is caught at the next reneighbor.
  void count_constrained_topology(bigint &natoms,
                                  bigint &nconstrbonds,
                                  bigint &nconstrangles) const;
  void record_topology_baseline();
  void check_topology_unchanged();

  void unconstrained_update();
  void unconstrained_update_respa(int ilevel);
  void stats();

  // solver
  void grow_lu_workspace(int n);
  int lu_factor_solve(int n);    // 0 = success, !=0 = singular
  // Cholesky factor + solve for symmetric positive-definite A.  Used for
  // the "fast" position-constraint matrix (truly symmetric: r_k.r_l
  // off-diagonals) and for the RATTLE-style velocity-projection matrix
  // (also symmetric).  Returns 0 on success, 1 if a non-positive pivot
  // is found (matrix is not SPD -- typically a degenerate constraint
  // cluster); the caller falls back to lu_factor_solve in that case.
  int chol_factor_solve(int n);
  // Split variants used by the cached-factor path in solve_constraints():
  // chol_factor does the in-place A = L L^T factorization on a caller-
  // supplied buffer; chol_solve does the forward+back substitution given
  // an already-computed L.  The merged chol_factor_solve above remains for
  // the velocity-projection path which has no factor to reuse.
  int chol_factor(int n, double *A);
  void chol_solve(int n, const double *L, double *b);

  // Banded Cholesky on the lower-band packed storage layout described at
  // the chol_pool member: rows of (bw+1) entries; row i column j (j <= i)
  // lives at AB[i*(bw+1) + (i-j)].  Diagonal at column-offset 0.  These
  // operate in-place on AB (factor) or on b (solve).  Cost: O(n*bw^2/2)
  // for factor, O(n*bw) for solve.
  int band_chol_factor(int n, int bw, double *AB);
  void band_chol_solve(int n, int bw, const double *AB, double *b);

  // Resize the per-cluster banded Cholesky factor cache.  Called from
  // precompute_constraint_data() once cluster_bw[] is filled.  Only
  // meaningful for the "fast" variant; the dense LU path uses lu_A directly.
  void grow_factor_cache();

  // Run per-cluster reverse Cuthill-McKee on the constraint-adjacency
  // graph (constraints share an edge iff they share an atom), and apply
  // the resulting permutation to c_perm and c_slot.  Fills cluster_bw[]
  // with the post-RCM bandwidth of each cluster.  Called from
  // precompute_constraint_data() before grow_factor_cache().
  void rcm_reorder_clusters();
  bool solve_constraints();
  void apply_constraint_forces(int vflag);
  void correct_coordinates(int vflag);
  void correct_velocities();
};

}    // namespace LAMMPS_NS

#endif
#endif
