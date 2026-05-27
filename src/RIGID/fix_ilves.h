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
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

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

  // Near-linear angle handling.  Angle constraints are implemented via the
  // law-of-cosines AC virtual-bond (|AC|^2 = r1^2 + r2^2 - 2 r1 r2 cos(theta_0)).
  // When theta_0 approaches 180 deg the {AB, BC, AC} triplet becomes
  // rank-deficient (|AC| -> |AB|+|BC|), so we silently DECLINE to add the
  // AC constraint for angle types with theta_0 >= linear_threshold (default
  // 165 deg) and leave the angle's force-field term in place (no sign
  // negation of angle_type[i][m] for these types).  The two flanking bond
  // constraints (AB, BC) are added normally when their bond types are
  // selected via the 'b' selector.  Set linearangle 180 to disable the
  // bailout; the AC constraint will then be added for every selected
  // angle, including ill-conditioned ones, and may fail to converge.
  double linear_threshold;        // degrees; angle types with theta_0 >=
                                  // this value get NO AC constraint
  double linear_angle_K;          // force constant of the optional stiff
                                  // angle potential applied to near-linear
                                  // angles in place of the user's
                                  // angle_style.  E = K*(1 + cos(theta));
                                  // no 1/sin(theta) singularity at 180.
                                  // Default 0 = disabled (angle_style
                                  // continues to handle the angle).
  int *angle_linear;              // angle_linear[at] = 1 if angle type at
                                  // is near-linear (size nangletypes+1)

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

  // flat constraint list (rebuilt each reneighbor); all entries are
  // 2-atom distance constraints.  c_atom1[k] is normally a local atom
  // (< nlocal) for ownership uniqueness; it may be a ghost for "cluster
  // completion" entries where every atom of the constraint is a ghost
  // on this rank.
  int n_constr;
  int max_constr;
  int *c_atom1;
  int *c_atom2;
  int *c_type;       // bond type for bond entries, or -atype for the
                     // angle-derived AC virtual entry
  double *c_dist;    // target distance (not squared) for clarity
  double *c_lambda;  // accumulated Lagrange multiplier (Phase 2+)

  // connected-component label for each constraint
  int *c_cluster;
  int n_clusters;

  // precomputed per-constraint reference quantities (rebuilt each reneighbor)
  // r_k = x[a1]-x[a2]; c_invma/c_invmb = 1/mass[a1] and 1/mass[a2].
  double *c_rx, *c_ry, *c_rz;
  double *c_rsq;
  double *c_invma, *c_invmb;

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

  // Newton iteration count diagnostics, accumulated between stats() calls.
  // newton_iter_sum = sum of iterations across all solve_constraints calls
  // newton_iter_max = max iterations in any single solve
  // newton_solve_count = number of solve_constraints calls
  // The next stats() prints (ave, max) over this interval and resets.
  bigint newton_iter_sum, newton_iter_max, newton_solve_count;

  // Topology-change safeguard.  Fix ilves builds the constraint list
  // from local bond/angle storage at init() and assumes it is frozen for
  // the duration of the run.  To catch silent mismatch when a concurrent
  // fix (e.g. fix bond/create, fix bond/break, delete_atoms) modifies the
  // constrained topology, we record three globally-reduced counts at init
  // and recheck them at every reneighbor; mismatch aborts with a clear
  // error.  The bond/angle counts are raw "sum over local atoms of
  // constrained-type slots" -- not deduplicated -- which is enough for
  // change detection.
  bigint natoms_at_init;
  bigint nconstrbonds_sum_at_init;
  bigint nconstrangles_sum_at_init;
  bool baseline_ready;

  // Per-angle-type flanking bond types, populated at init.  Each angle of
  // type at has flanking bond types (b1, b2) sorted so b1 <= b2.  Computed
  // per-rank from local angles and consensus-checked via MPI_Allreduce so
  // every rank agrees.  Used in build_constraint_list to verify that an
  // angle's flanking bonds are themselves constrained without consulting a
  // global table.
  int *angle_btype1, *angle_btype2;     // size nangletypes+1

  // -----------------------------------------------------------------
  // build helpers
  void grow_constraint_list(int);
  void add_constraint(int a, int b, int btype, double dist);
  void group_by_cluster();
  void precompute_constraint_data();
  int masscheck(double massone);
  bool bond_selected_for_atoms(int ia, int ib, int bt);

  void negate_constrained_topology();

  // global topology setup
  void init_topology();
  void build_constraint_list();
  bool bond_is_constrained(tagint ta, tagint tb);
  int  lookup_local_bond_type(tagint ta, tagint tb);

  // Refresh ilv_bond_* per-fix arrays from atom->bond_* and forward_comm
  // the ghost-side entries.  Called from build_constraint_list at every
  // reneighbor (cheap once-per-reneighbor cost, not per Newton iter).
  void refresh_ilv_bond_data();

  // Optional stiff angle force applied to near-linear constrained angle
  // types when linear_angle_K > 0.  Uses E = K*(1 + cos(theta)), the
  // standard "cosine" angle form -- no 1/sin singularity at theta=180.
  // Replaces the user's angle_style for these angle types (the
  // angle_type is negated so angle_style->compute skips them).
  void apply_linear_angle_forces(int vflag);

  // Schwarz overlap: ghost-side bond storage.  Without overlap each rank's
  // local Jacobian misses off-diagonal couplings to constraints stored at
  // remote endpoints, so Newton degenerates to slow Schwarz iteration on
  // cross-rank clusters.  We refresh these arrays at every reneighbor by
  // copying from atom->bond_* for local atoms and forward_comming to the
  // ghost copies.  build_constraint_list then walks ALL atoms in halo
  // (local + ghost) via these arrays, so each rank's constraint list
  // includes bonds fully inside its halo (even ghost-only bonds) -- the
  // Jacobian Schwarz subdomain overlaps by one ghost-shell layer with
  // neighbor ranks.
  int    *ilv_num_bond;        // per-atom bond count for local + ghosts
  tagint *ilv_bond_atom;       // [nmax * ilv_bond_per_atom] partner tags
  int    *ilv_bond_type;       // [nmax * ilv_bond_per_atom] bond types
  int     ilv_bond_per_atom;   // size of inner dimension
  int     ilv_nmax_alloc;      // current outer (nmax) allocation
  void grow_ilv_bond(int nmax);

  // Warm-start of c_lambda from the previous step.  Constraint indexing
  // is stable between reneighbor calls, so storing c_lambda by index is
  // safe within a reneighbor window.  build_constraint_list() resets
  // lambda_warm_valid to 0 on reneighbor; solve_constraints() sets it to
  // 1 after a successful Newton convergence.  When valid, the next call
  // initializes c_lambda from the saved values (skip the zeroing) and
  // applies the corresponding xshake correction before the Newton loop.
  //
  // Off by default (preserves bit-reproducibility against the YAML
  // reference set).  Enable with 'warmstart yes' in the fix command:
  // changes the iteration path and produces ULP-level (~5e-11) trajectory
  // differences vs cold-start over 100 steps, but cuts iter counts in
  // multi-rank Schwarz runs.
  int lambda_warm_valid;
  int warmstart_enabled;

  // Newton-update damping factor.  Default 1.0 (no damping).  Lower
  // values reduce |dlambda| applied per iteration -- useful when
  // Schwarz iterations between ranks oscillate.  Set via the optional
  // 'relax' keyword in the fix command.
  double newton_relax;

  // Per-fix force buffer used by the stiff angle path under newton on
  // bond.  Sized to nmax * 3.  In apply_linear_angle_forces we write the
  // angle force trio (for all three atoms, some of which may be ghosts)
  // into this buffer, then reverse_comm to send ghost contributions to
  // owner ranks (adding into the owner's local entries), and finally
  // add the local-sum buffer into atom->f.  Avoids touching atom->f's
  // ghost portion (which other code may rely on being zero after the
  // standard force reverse_comm).
  double *lang_fbuf;
  int lang_fbuf_alloc;
  void grow_lang_fbuf(int nmax);

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
