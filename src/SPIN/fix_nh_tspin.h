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

/* ------------------------------------------------------------------------
   Contributing author: Zhengtao Huang (The University of Hong Kong)
                        hzt990224@gmail.com
------------------------------------------------------------------------- */

#ifndef LMP_FIX_NH_TSPIN_H
#define LMP_FIX_NH_TSPIN_H

#include "fix_nh.h"

namespace LAMMPS_NS {

class FixNHTSpin : public FixNH {
 public:
  FixNHTSpin(class LAMMPS *, int, char **);
  ~FixNHTSpin() override;

  void post_constructor() override;
  void init() override;
  void setup(int) override;
  double compute_scalar() override;
  void reset_dt() override;
  void write_restart(FILE *) override;
  int size_restart_global() override;
  void restart(char *) override;

  // read by PairSpin::init_style()

  int lattice_flag;

 protected:
  int spin_flag;        // 1 if the spin degrees of freedom are integrated
  int spinmass_flag;    // 1 if the spinmass keyword was used, then this fix owns s_mass
  double spinmass;      // spin mass in units of the atomic mass of the type
  int index_vs;         // custom per-atom property holding the spin velocity
  int index_sm;         // custom per-atom property holding the spin mass
  double hbar;          // Planck constant over 2 pi, in eV/(rad.THz)

  double t_current_spin, ke_target_spin, tdof_spin;

  // Nose-Hoover chain acting on the spin velocities, parallel to the chain
  // of the parent class and driven by the same target temperature

  double *etas, *etas_dot, *etas_dotdot, *etas_mass;

  void nve_v() override;
  void nve_x() override;
  void nh_v_temp() override;
  void nhc_temp_integrate() override;

  void nhc_spin_integrate();
  double compute_spin_temp();
  void nh_vs_scale(double);
  void set_spin_mass();
  void spin_kick();
};

}    // namespace LAMMPS_NS

#endif
