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
  void pre_neighbor() override;
  void post_force(int) override;
  void end_of_step() override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

  bigint dof(int) override;

 protected:
  // user settings
  double tolerance;       // convergence tolerance on relative bond-length error
  int max_iter;           // max # of Newton iterations per step
  int output_every;       // print constraint statistics every this many steps (0 = never)
  bigint next_output;     // next timestep for statistics output
  int variant;            // ILVES_FAST (symmetric) or ILVES_FULL (asymmetric)
  int fixed_iter;         // 1 = run exactly max_iter iterations (no convergence test)

  // selectors (which bonds/angles to constrain), as in fix shake
  std::vector<int> bond_flag;     // [nbondtypes+1]  constrain these bond types
  std::vector<int> angle_flag;    // [nangletypes+1] constrain these angle types (deferred)
  std::vector<int> type_flag;     // [ntypes+1] constrain bonds touching these atom types
  std::vector<double> mass_list;  // constrain bonds touching atoms of these masses

  int molecular;                     // copy of atom->molecular
  std::vector<double> bond_distance; // [nbondtypes+1] equilibrium bond lengths
  int types_negated;                 // 1 once constrained bond types have been negated

  int store_flag;     // 1 to expose per-atom constraint forces via array_atom
  double **fstore;    // per-atom constraint forces (when store_flag)
  int maxstore;       // current allocated length of fstore
  int niter_max;      // max Newton iterations used since the last stats output

  // local constraint list, rebuilt every reneighbor.
  // constraint k joins local/ghost atoms clist_a[k] and clist_b[k] with target
  // distance clist_d[k]; clist_btype[k] is the (positive) bond type for stats.
  int nconstraints;
  std::vector<int> clist_a, clist_b, clist_btype;
  std::vector<double> clist_d;

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
  int maxatom;    // current allocated length of the above per-atom arrays
  int commstage;  // 0 = forward-comm positions (PBC shift), 1 = velocities (no shift)

  // timestep factors, as in fix shake
  double dtv;          // = dt
  double dtfsq;        // = dt^2 * ftm2v
  double inv_dtfsq;    // = 1 / dtfsq

  // cached local ptrs to atom-class quantities
  double **x, **v, **f;
  double *mass, *rmass;
  int *type, *mask;
  int nlocal;

  void build_constraint_list();
  void project_velocities();
  void grow_arrays_local();
  void stats();
  void negate_bond_types(int sign);    // sign < 0 negate, sign > 0 restore
  int bond_selected(int i, int j, int btype);
  int masscheck(double massone);
};

}    // namespace LAMMPS_NS

#endif
#endif
