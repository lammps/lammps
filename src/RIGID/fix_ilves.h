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

#include <vector>

namespace LAMMPS_NS {

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
  int *cluster_offset;           // size n_clusters+1, CSR-style ranges into c_perm[]
  int *c_perm;                   // size n_constr; c_perm[s] = global constraint id at slot s
  int *c_slot;                   // size n_constr; inverse map: c_slot[k] = slot of constraint k

  // Newton workspace (sized to largest cluster)
  int lu_alloc;
  double *lu_A;                  // row-major (lu_alloc * lu_alloc)
  double *lu_b;                  // rhs / dlambda
  int *lu_pivot;
  double *cl_sx, *cl_sy, *cl_sz; // s_k cache for the current cluster
  int largest_cluster;           // size of the largest cluster (informational)

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

  // remember vflag from post_force for end-of-step bookkeeping (Phase 2+)
  int vflag_post_force;
  int eflag_pre_reverse;

  // selector for pack_forward_comm: 0 = pack xshake, 1 = pack atom->v
  int comm_mode;

  // Global bond / angle topology, gathered once at init.  Used to build the
  // local constraint list so that every rank with at least one atom of a
  // bond / angle picks it up -- this makes the solver correct under MPI
  // for both newton on and newton off, even when bond storage is one-sided.
  // Encoded as (lower-tag, higher-tag, type) for bonds and as
  // (atom1_tag, mid_tag, atom3_tag, type) for angles, sorted and deduplicated.
  std::vector<tagint> gb_a, gb_b;
  std::vector<int> gb_type;
  std::vector<tagint> ga1, ga2, ga3;
  std::vector<int> ga_type;
  bool global_topology_ready;

  // global cluster-ID for every tag involved in any bond/angle.
  // tag_to_cluster_tag[i] = representative tag of the cluster that
  // global tag (i+1) belongs to.  Tags that don't belong to any cluster
  // get -1.  Sized to atom->natoms+1 (indexed 1..natoms).
  // This lets every rank quickly find which clusters intersect its
  // local atoms, and pull in ALL constraints from those clusters --
  // including constraints between two ghosts (needed for cluster
  // completeness when the cluster spans multiple ranks).
  std::vector<tagint> tag_cluster;     // size natoms+1; -1 if no cluster

  // helpers
  void gather_global_topology();
  void build_constraint_list();
  void group_by_cluster();
  void precompute_constraint_data();
  void grow_constraint_list(int);
  void add_constraint(int a, int b, int btype, double dist);
  int masscheck(double massone);
  void negate_constrained_topology();
  int bondtype_findset(int i, tagint n1, tagint n2, int setflag);
  int angletype_findset(int i, tagint n1, tagint n2, int setflag);
  bool bond_selected_for_atoms(int ia, int ib, int bt);
  bool bond_is_constrained_global(tagint ta, tagint tb);
  // ilves variant: 0 = symmetric (use r.r in off-diagonals -> approx Jacobian)
  //                1 = structurally-symmetric (use s.r -> exact Newton)
  // Both converge to the exact solution of the constraint equations.
  virtual void unconstrained_update();
  virtual void stats();

  // solver
  void grow_lu_workspace(int n);
  int lu_factor_solve(int n);    // 0 = success, !=0 = singular
  bool solve_constraints();
  void apply_constraint_forces(int vflag);
  void correct_coordinates(int vflag);
  void correct_velocities();
};

}    // namespace LAMMPS_NS

#endif
#endif
