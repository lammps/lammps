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
  void initial_integrate(int) override;
  void post_force(int) override;
  void end_of_step() override;
  int modify_param(int, char **) override;

 protected:
  int dimension, which;

  int pstyle, pcouple, allremap;
  int p_flag[3];    // 1 if control P on this dim, 0 if not
  double t_start, t_stop, t_target;
  double p_fric;
  double p_start[3], p_stop[3], p_current[3];
  double p_period[3], p_target[3];
  double p_deriv[3], dilation[3];
  double f_piston[3], f_old_piston[3];
  double gjfa[3], gjfb[3], fran[3];
  int kspace_flag;    // 1 if KSpace invoked, 0 if not
  int nrigid;         // number of rigid fixes
  int *rfix;          // indices of rigid fixes

  char *id_temp, *id_press;
  class Compute *temperature, *pressure;
  int pflag;

  class RanMars *random;
  int seed;

  void couple_pressure();
  void couple_kinetic(double );
  void couple_beta(double );
  void remap();
};

}    // namespace LAMMPS_NS

#endif
#endif
