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

#ifndef LMP_FIX_ILVES_H
#define LMP_FIX_ILVES_H

#include "fix.h"

namespace LAMMPS_NS {

// Abstract base class for the ILVES bond/angle constraint solver.
// Concrete subclasses (FixIlvesLocal, FixIlvesGlobal) supply the
// constraint-discovery strategy via the three protected pure-virtual
// hooks (init_topology, build_constraint_list, bond_is_constrained).
// The base owns the flat per-rank constraint list plus all solver,
// communication, and bookkeeping code.

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
  int variant;           // 0 = ilves (symmetric), 1 = ilvesf (struct. symmetric)
  bigint next_output;

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

  // flat constraint list (rebuilt each reneighbor by the derived class)
  // c_atom1[k] is always a local atom (< nlocal) for ownership uniqueness
  // c_atom2[k] is local or ghost
  int n_constr;
  int max_constr;
  int *c_atom1;
  int *c_atom2;
  int *c_type;       // bond type for bond entries, or -atype for A-C "virtual"
                     // entries from an angle (negative tags angle origin)
  double *c_dist;    // target distance (not squared) for clarity
  double *c_lambda;  // accumulated Lagrange multiplier (Phase 2+)

  // connected-component label for each constraint
  int *c_cluster;
  int n_clusters;

  // precomputed per-constraint reference quantities (rebuilt each reneighbor)
  double *c_rx, *c_ry, *c_rz;    // r_k = x[a_k] - x[b_k] (closest image)
  double *c_rsq;
  double *c_invma, *c_invmb;     // inverse masses of a_k, b_k

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

  // Per-cluster Cholesky factor cache (ILVES_FAST only).  The cluster
  // Jacobian in ILVES_FAST is built entirely from r_k, c_rsq, and inverse
  // masses -- all step-constant within solve_constraints() -- so its
  // Cholesky factor can be computed once at the first Newton iteration
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

  // if false (local variant), stats() counts every constraint in the list
  // without deduplication; if true (global variant), it deduplicates by
  // lower-tag ownership because the same constraint lives on multiple ranks.
  bool stats_dedup;

  // counters incremented per cluster-solve; reported on rank 0 when
  // output_every > 0 so the user sees how often we fall back to LU.
  bigint chol_calls, chol_fallbacks;

  // ---- pure virtual hooks (variant-specific) ----

  // Variant-specific one-time topology setup at init.  The local variant
  // performs a local sanity scan and computes angle_distance from local
  // angle storage; the global variant gathers the full bond/angle table
  // and computes angle_distance from it.  Called once from FixIlves::init()
  // after bond_distance[] has been filled.
  virtual void init_topology() = 0;

  // Variant-specific construction of the per-rank flat constraint list
  // (c_atom1/c_atom2/c_type/c_dist).  Called from post_neighbor.  The
  // implementation is expected to call grow_constraint_list / add_constraint
  // and then group_by_cluster + precompute_constraint_data at the end.
  virtual void build_constraint_list() = 0;

  // Variant-specific query: is the bond between global tags (ta, tb)
  // selected for constraint?  Used by negate_constrained_topology to
  // decide which angle entries to negate (both flanking bonds must be
  // themselves constrained).  Implementations may assume both atoms are
  // reachable locally (own or ghost).
  virtual bool bond_is_constrained(tagint ta, tagint tb) = 0;

  // ---- shared helpers used by derived classes ----

  void grow_constraint_list(int);
  void add_constraint(int a, int b, int btype, double dist);
  void group_by_cluster();
  void precompute_constraint_data();
  int masscheck(double massone);
  bool bond_selected_for_atoms(int ia, int ib, int bt);

  void negate_constrained_topology();

  // ilves variant: 0 = symmetric (use r.r in off-diagonals -> approx Jacobian)
  //                1 = structurally-symmetric (use s.r -> exact Newton)
  // Both converge to the exact solution of the constraint equations.
  virtual void unconstrained_update();
  virtual void stats();

  // solver
  void grow_lu_workspace(int n);
  int lu_factor_solve(int n);    // 0 = success, !=0 = singular
  // Cholesky factor + solve for symmetric positive-definite A.  Used for the
  // ILVES_FAST position-constraint matrix (truly symmetric: r_k.r_l off-
  // diagonals) and for the RATTLE-style velocity-projection matrix (also
  // symmetric).  Returns 0 on success, 1 if a non-positive pivot is found
  // (matrix is not SPD -- typically a degenerate constraint cluster); the
  // caller falls back to lu_factor_solve in that case.
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
  // meaningful for ILVES_FAST; the dense LU path uses lu_A directly.
  void grow_factor_cache();

  // Run per-cluster reverse Cuthill-McKee on the constraint-adjacency
  // graph (constraints share an edge iff they share an atom), and apply
  // the resulting permutation to c_perm and c_slot.  Fills cluster_bw[]
  // with the post-RCM bandwidth of each cluster.  Called from
  // precompute_constraint_data() before grow_factor_cache().
  void rcm_reorder_clusters();
  bool solve_constraints();
  virtual void apply_constraint_forces(int vflag);
  virtual void correct_coordinates(int vflag);
  virtual void correct_velocities();

  // Per-iteration xshake sync hook inside solve_constraints.  The global
  // variant calls comm->forward_comm to keep ghost xshake in sync across
  // ranks doing redundant cluster solves; the local variant (which uses
  // single-rank cluster ownership) overrides this to a no-op, because
  // each owner maintains the ghost xshake locally and the leaves' owners
  // don't need to be updated.
  virtual void sync_xshake();

  // Should the cluster owner update ghost atoms' xshake during Newton
  // iteration?  True for the local variant (single-rank cluster
  // ownership: the owner maintains the cluster's ghost positions
  // locally because the ghost-owners don't process the cluster).
  // False for the global variant (redundant solve: ghost xshake gets
  // overwritten by forward_comm anyway).
  virtual bool update_ghost_xshake() const { return false; }
};

}    // namespace LAMMPS_NS

#endif
