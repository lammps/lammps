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

namespace ILVES {
class Ilves;
}

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

 protected:
  // user settings
  double tolerance;       // convergence tolerance on relative bond-length error
  int max_iter;           // max # of Newton iterations per step
  int output_every;       // print constraint statistics every this many steps (0 = never)
  bigint next_output;     // next timestep for statistics output
  int variant;            // ILVES_FAST (symmetric) or ILVES_FULL (asymmetric)

  // selectors (which bonds/angles to constrain), as in fix shake
  std::vector<int> bond_flag;     // [nbondtypes+1]  constrain these bond types
  std::vector<int> angle_flag;    // [nangletypes+1] constrain these angle types (deferred)
  std::vector<int> type_flag;     // [ntypes+1] constrain bonds touching these atom types
  std::vector<double> mass_list;  // constrain bonds touching atoms of these masses

  int molecular;                     // copy of atom->molecular
  std::vector<double> bond_distance; // [nbondtypes+1] equilibrium bond lengths
  int types_negated;                 // 1 once constrained bond types have been negated

  // local constraint list, rebuilt every reneighbor.
  // constraint k joins local/ghost atoms clist_a[k] and clist_b[k] with target
  // distance clist_d[k]; clist_btype[k] is the (positive) bond type for stats.
  int nconstraints;
  std::vector<int> clist_a, clist_b, clist_btype;
  std::vector<double> clist_d;

  // the ported ILVES solver, rebuilt every reneighbor from the constraint list
  ILVES::Ilves *ilves_solver;

  // per-atom inverse mass (1/m) handed to the solver; kept alive across steps
  std::vector<double> invmass;
  // predicted (unconstrained) positions, iterated by the solver
  double **xpred;
  int maxatom;    // current allocated length of xpred

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
  void negate_bond_types(int sign);    // sign < 0 negate, sign > 0 restore
  int bond_selected(int i, int j, int btype);
  int masscheck(double massone);
};

}    // namespace LAMMPS_NS

#endif
#endif
