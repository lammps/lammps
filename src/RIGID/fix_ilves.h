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

#include "ilves.h"

#include <vector>

namespace LAMMPS_NS {

class FixIlves : public Fix {
 public:
  FixIlves(class LAMMPS *, int, char **);
  ~FixIlves() override;
  int setmask() override;
  void init() override;
  void setup_pre_neighbor() override;
  void setup(int) override;
  void min_setup(int) override;
  void pre_neighbor() override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

  bigint dof(int) override;
  double compute_scalar() override;
  double memory_usage() override;

 protected:
  // user settings
  double tolerance;        // convergence tolerance on relative bond-length error
  int max_iter;            // max # of Newton iterations per step
  int output_every;        // print constraint statistics every this many steps (0 = never)
  bigint next_output;      // next timestep for statistics output
  bool fixed_iter_flag;    // true = run exactly max_iter iterations (no convergence test)

  // near-linear angle handling.  the A-C "virtual bond" of an angle becomes
  // rank-deficient as the equilibrium angle approaches 180 degrees, so angle
  // types whose theta0 is at or above linear_threshold (in degrees) are treated
  // specially according to linear_mode.
  int linear_mode;            // LINEAR_ERROR, LINEAR_SKIP, or LINEAR_RESTRAIN
  double linear_threshold;    // degrees; angle types with theta0 >= this are near-linear
  double kbond;               // force constant for the restrain substitute (< 0 = auto)
  double erestraint;          // potential energy of the restrain substitute (this rank)

  // selectors (which bonds/angles to constrain), as in fix shake
  std::vector<int> bond_flag;       // [nbondtypes+1]  constrain these bond types
  std::vector<int> angle_flag;      // [nangletypes+1] constrain these angle types (deferred)
  std::vector<int> type_flag;       // [ntypes+1] constrain bonds touching these atom types
  std::vector<double> mass_list;    // constrain bonds touching atoms of these masses

  int molecular;                         // copy of atom->molecular
  std::vector<double> bond_distance;     // [nbondtypes+1] equilibrium bond lengths
  std::vector<double> angle_distance;    // [nangletypes+1] equilibrium A-C distances
  std::vector<int> angle_linear;    // [nangletypes+1] 1 if type is near-linear (theta0>=threshold)
  int types_negated;                // 1 once constrained bond/angle types have been negated

  int store_flag;     // 1 to expose per-atom constraint forces via array_atom
  double **fstore;    // per-atom constraint forces (when store_flag)
  int maxstore;       // current allocated length of fstore
  int niter_max;      // max Newton iterations used since the last stats output

  // local constraint list, rebuilt every reneighbor.
  // constraint k joins local/ghost atoms clist_a[k] and clist_b[k] with target
  // distance clist_d[k].  clist_a/clist_b are GEOMETRY indices (nearest periodic
  // image, for the bond vector); clist_node_a/clist_node_b are the matching
  // canonical (owner) NODE ids (an atom and its ghost image share one node id),
  // used so a periodically wrapped bond stays a single graph edge.
  // clist_btype[k] > 0 is the bond type; clist_btype[k] < 0 marks an angle A-C
  // "virtual bond" of angle type -clist_btype[k], whose vertex (center) atom is
  // clist_vertex[k] (used only to report the actual angle in the statistics
  // output, as fix shake does); clist_vertex[k] is -1 for bonds.
  int nconstraints;
  std::vector<int> clist_a, clist_b, clist_node_a, clist_node_b, clist_btype, clist_vertex;
  std::vector<double> clist_d;

  // near-linear angles handled by the restrain substitute, rebuilt every
  // reneighbor: a stiff harmonic bond on the A-C "virtual bond" between the
  // outer atoms rlist_a/rlist_c (local or ghost) with target distance rlist_d.
  std::vector<int> rlist_a, rlist_c;
  std::vector<double> rlist_d;

  // the ported ILVES solver, rebuilt every reneighbor from the constraint list
  ILVES::Ilves *ilves_solver;

  // per-atom inverse mass (1/m) handed to the solver; sized nlocal+nghost
  std::vector<double> invmass;
  // predicted (unconstrained) positions, iterated by the solver (home+ghost),
  // the saved unconstrained prediction (home), and the per-iteration increment
  // accumulator (home+ghost) used for the reverse-sum across ranks
  double **xpred;
  double **xpred0;
  double **dx;
  int maxatom;      // current allocated length of the above per-atom arrays
  int commstage;    // 0 = forward-comm positions (PBC shift), 1 = velocities (no shift)

  // timestep factors, as in fix shake.  for r-RESPA dtfsq/inv_dtfsq are
  // recomputed per level in post_force_respa (= dtf_inner * step_respa[ilevel]).
  double dtv;          // = dt (or step_respa[0] for r-RESPA)
  double dtfsq;        // = dt^2 * ftm2v
  double inv_dtfsq;    // = 1 / dtfsq

  // r-RESPA support, mirroring fix shake (SHAKE-form coupling, so the full-step
  // dtf_inner is used, not the half-step rattle form)
  int respa;                    // 0 = velocity Verlet, 1 = r-RESPA
  int nlevels_respa;            // number of r-RESPA levels
  int *loop_respa;              // sub-cycling factor at each level (Respa-owned)
  double *step_respa;           // timestep at each level (Respa-owned)
  class FixRespa *fix_respa;    // holds the per-level forces (f_level)
  double dtf_inner;             // = step_respa[0] * ftm2v
  double dtf_innerhalf;         // = 0.5 * step_respa[0] * ftm2v

  // cached local ptrs to atom-class quantities
  double **x, **v, **f;
  double *mass, *rmass;
  int *type, *mask;
  int nlocal;

  void build_constraint_list();
  int run_newton();           // drive prepared xpred onto the constraint manifold; returns # iters
  int solve_constraints();    // run_newton, then turn the displacement into forces; returns # iters
  void apply_linear_restraint();
  double min_harmonic_bond(int a, int b, double d, double k);
  void project_velocities();
  void grow_arrays_local();
  void stats();
  void negate_bond_types(int sign);     // sign < 0 negate, sign > 0 restore
  void negate_angle_types(int sign);    // sign < 0 negate, sign > 0 restore
  int bond_selected(int i, int j, int btype);
  int angle_selected(int i, int m, int &a, int &c, int &atype);
  int find_bond_type(int i, int j);
  int masscheck(double massone);
};

}    // namespace LAMMPS_NS

#endif
#endif
