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

  enum { NONE, FINAL, DELTA, SCALE, VEL, ERATE, TRATE, VOLUME, WIGGLE, VARIABLE, PRESSURE, PMEAN };
  enum { ONE_FROM_ONE, ONE_FROM_TWO, TWO_FROM_ONE };

  FixDeform(class LAMMPS *, int, char **);
  ~FixDeform() override;
  int setmask() override;
  void init() override;
  void pre_exchange() override;
  void end_of_step() override;
  void virtual write_restart(FILE *) override;
  void virtual restart(char *buf) override;
  double memory_usage() override;

 protected:
  int triclinic, scaleflag, flipflag;
  int flip, flipxy, flipxz, flipyz;
  double *h_rate, *h_ratelo;
  int varflag;                   // 1 if VARIABLE option is used, 0 if not
  int kspace_flag;               // 1 if KSpace invoked, 0 if not
  std::vector<Fix *> rfix;       // pointers to rigid fixes
  class Irregular *irregular;    // for migrating atoms after box flips

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
    int fixed, dynamic1, dynamic2;
    char *hstr, *hratestr;
    int hvar, hratevar;
  };
  Set *set;

  std::vector<int> leftover_iarg;
  int iarg_options_start;

  void options(int, char **);
  void virtual apply_volume();
  void apply_strain();
  void update_domain();
};

}    // namespace LAMMPS_NS

#endif
#endif
