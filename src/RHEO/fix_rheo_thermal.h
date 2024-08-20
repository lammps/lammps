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
FixStyle(rheo/thermal,FixRHEOThermal)
// clang-format on
#else

#ifndef LMP_FIX_RHEO_THERMAL_H
#define LMP_FIX_RHEO_THERMAL_H

#include "fix.h"

#include <vector>

namespace LAMMPS_NS {

class FixRHEOThermal : public Fix {
 public:
  FixRHEOThermal(class LAMMPS *, int, char **);
  ~FixRHEOThermal() override;
  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void setup_pre_force(int) override;
  void initial_integrate(int) override;
  void post_integrate() override;
  void post_neighbor() override;
  void pre_force(int) override;
  void final_integrate() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  void reset_dt() override;
  double calc_cv(int);
  double calc_Tc(int);
  double calc_L(int);

 private:
  double *cv, *Tc, *kappa, *L;
  double dt, dth;
  double cut_kernel, cut_bond, cutsq_bond;
  int *cv_style, *Tc_style, *kappa_style, *L_style;
  int btype;
  class NeighList *list;

  int n_histories;
  std::vector<Fix *> histories;

  class FixRHEO *fix_rheo;
  class ComputeRHEOGrad *compute_grad;
  class ComputeRHEOVShift *compute_vshift;
  class FixUpdateSpecialBonds *fix_update_special_bonds;

  void grow_array(int);
  void break_bonds();
  void create_bonds();
};

}    // namespace LAMMPS_NS

#endif
#endif
