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
FixStyle(press/langevin,FixPressLangevin);
// clang-format on
#else

#ifndef LMP_FIX_PRESS_LANGEVIN_H
#define LMP_FIX_PRESS_LANGEVIN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPressLangevin : public Fix {
 public:
  FixPressLangevin(class LAMMPS *, int, char **);
  ~FixPressLangevin() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void pre_exchange() override;
  void initial_integrate(int) override;
  void post_integrate() override;
  void post_force(int) override;
  void end_of_step() override;
  void reset_dt() override;
  int modify_param(int, char **) override;

 protected:
  int dimension;
  int pstyle, pcouple, allremap;
  int p_flag[6];    // 1 if control P on this dim, 0 if not
  double t_start, t_stop, t_target;
  double p_fric[6], p_ltime;    // Friction and Langevin charac. time
  double p_alpha[6];
  double p_start[6], p_stop[6], p_period[6];
  double p_mass[6], p_target[6], p_current[6];
  double p_deriv[6], dilation[6];
  double f_piston[6], f_old_piston[6];
  double gjfa[6], gjfb[6], fran[6];
  int kspace_flag;            // 1 if KSpace invoked, 0 if not
  std::vector<Fix *> rfix;    // indices of rigid fixes

  char *id_temp, *id_press;
  class Compute *temperature, *pressure;
  int pflag;

  int flipflag;
  int pre_exchange_flag;         // set if pre_exchange needed for box flips
  class Irregular *irregular;    // for migrating atoms after box flips

  class RanMars *random;
  int seed;

  void couple_pressure();
  void couple_kinetic();
  void couple_beta();
  void remap();
};

}    // namespace LAMMPS_NS

#endif
#endif
