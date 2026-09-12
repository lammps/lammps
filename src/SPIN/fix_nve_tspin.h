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

#ifdef FIX_CLASS
// clang-format off
FixStyle(nve/tspin,FixNVETSpin);
// clang-format on
#else

#ifndef LMP_FIX_NVE_TSPIN_H
#define LMP_FIX_NVE_TSPIN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNVETSpin : public Fix {
 public:
  FixNVETSpin(class LAMMPS *, int, char **);

  void post_constructor() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void reset_dt() override;

  // read by PairSpin::init_style()

  int lattice_flag;

 protected:
  int spin_flag;      // 1 if the spin degrees of freedom are integrated
  int spinmass_flag;  // 1 if the spinmass keyword was used, then this fix owns s_mass
  double spinmass;    // spin mass in units of the atomic mass of the type
  int index_vs;       // custom per-atom property holding the spin velocity
  int index_sm;       // custom per-atom property holding the spin mass
  double dtv, dtf;
  double hbar;    // Planck constant over 2 pi, in eV/(rad.THz)

  void set_spin_mass();
  void check_spin_mass();
  void spin_kick();
};

}    // namespace LAMMPS_NS

#endif
#endif
