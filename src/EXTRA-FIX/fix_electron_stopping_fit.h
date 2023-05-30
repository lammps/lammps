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
   Contributing authors: Stephen M. Foiles (SNL)
                         James A. Stewart (SNL)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(electron/stopping/fit,FixElectronStoppingFit);
// clang-format on
#else

#ifndef LMP_FIX_ELECTRON_STOPPING_FIT_H
#define LMP_FIX_ELECTRON_STOPPING_FIT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixElectronStoppingFit : public Fix {
 public:
  FixElectronStoppingFit(class LAMMPS *, int, char **);
  ~FixElectronStoppingFit() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  double compute_scalar() override;

 private:
  double *energy_coh_in, *v_min_sq, *v_max_sq, *drag_fac_in_1, *drag_fac_in_2, *drag_fac_1,
      *drag_fac_2;
  double electronic_loss, electronic_loss_this_node;
  double f_dot_v_prior, f_dot_v_current;
  int last_step, this_step;
  int nlevels_respa;
};

}    // namespace LAMMPS_NS

#endif
#endif
