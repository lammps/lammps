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
FixStyle(deform/pressure,FixDeformPressure);
// clang-format on
#else

#ifndef LMP_FIX_DEFORM_PRESSURE_H
#define LMP_FIX_DEFORM_PRESSURE_H

#include "fix_deform.h"

namespace LAMMPS_NS {

class FixDeformPressure : public FixDeform {
 public:
  FixDeformPressure(class LAMMPS *, int, char **);
  ~FixDeformPressure() override;
  void init() override;
  void setup(int) override;
  void end_of_step() override;
  void write_restart(FILE *) override;
  void restart(char *buf) override;
  int modify_param(int, char **) override;

 protected:
  int pcouple;
  double max_h_rate;
  int strain_flag;               // 1 if strain-based option is used, 0 if not
  int pressure_flag;             // 1 if pressure tensor used, 0 if not
  int volume_flag;               // 1 if VOLUME option is used, 0 if not
  int normalize_pressure_flag;   // 1 if normalize pressure deviation by target
  int vol_balance_flag;          // 1 if pressures balanced when maintaining const vol

  char *id_temp, *id_press;
  class Compute *temperature, *pressure;
  int tflag, pflag;

  struct SetExtra {
    double ptarget, pgain;
    double prior_pressure, prior_rate;
    double cumulative_shift;
    double cumulative_vshift[3];
    double cumulative_remap;
    int saved;
    char *pstr;
    int pvar, pvar_flag;
    int coupled_flag;
  };
  SetExtra *set_extra;
  Set set_box;

  void options(int, int, char **);
  void apply_volume() override;
  void apply_pressure();
  void apply_box();
  void couple();
  void adjust_linked_rates(double&, double&, double, double, double);
};

}    // namespace LAMMPS_NS

#endif
#endif
