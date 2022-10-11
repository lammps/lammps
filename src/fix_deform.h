/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(deform,FixDeform);
// clang-format on
#else

#ifndef LMP_FIX_DEFORM_H
#define LMP_FIX_DEFORM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDeform : public Fix {
 public:
  int remapflag;     // whether x,v are remapped across PBC
  int dimflag[6];    // which dims are deformed

  FixDeform(class LAMMPS *, int, char **);
  ~FixDeform() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void pre_exchange() override;
  void end_of_step() override;
  void write_restart(FILE *) override;
  void restart(char *buf) override;
  double memory_usage() override;
  int modify_param(int, char **) override;

 protected:
  int dimension, triclinic, scaleflag, flipflag, pcouple;
  int flip, flipxy, flipxz, flipyz;
  double *h_rate, *h_ratelo;
  int varflag;                   // 1 if VARIABLE option is used, 0 if not
  int strain_flag;               // 1 if strain-based option is used, 0 if not
  int volume_flag;               // 1 if VOLUME option is used, 0 if not
  int pressure_flag;             // 1 if pressure tensor used, 0 if not
  int kspace_flag;               // 1 if KSpace invoked, 0 if not
  std::vector<Fix *> rfix;       // pointers to rigid fixes
  class Irregular *irregular;    // for migrating atoms after box flips

  double TWOPI;

  char *id_temp, *id_press;
  class Compute *temperature, *pressure;
  int tflag, pflag;

  struct Set {
    int style, substyle;
    double flo, fhi, ftilt;
    double dlo, dhi, dtilt;
    double scale, vel, rate;
    double amplitude, tperiod;
    double lo_initial, hi_initial;
    double lo_start, hi_start, lo_stop, hi_stop, lo_target, hi_target;
    double tilt_initial, tilt_start, tilt_stop, tilt_target, tilt_flip;
    double tilt_min, tilt_max;
    double vol_initial, vol_start;
    double ptarget, pgain;
    double prior_pressure, prior_rate;
    double box_length;
    int saved;
    int fixed, dynamic1, dynamic2;
    char *hstr, *hratestr, *pstr;
    int hvar, hratevar;
    int pvar, pvar_flag;
    int coupled_flag;
  };
  Set *set;

  void options(int, char **);
  void set_strain();
  void set_pressure();
  void set_volume();
  void couple();
};

}    // namespace LAMMPS_NS

#endif
#endif
