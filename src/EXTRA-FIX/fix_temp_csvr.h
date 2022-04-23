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
FixStyle(temp/csvr,FixTempCSVR);
// clang-format on
#else

#ifndef LMP_FIX_TEMP_CSVR_H
#define LMP_FIX_TEMP_CSVR_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTempCSVR : public Fix {
 public:
  FixTempCSVR(class LAMMPS *, int, char **);
  ~FixTempCSVR() override;
  int setmask() override;
  void init() override;
  void end_of_step() override;
  int modify_param(int, char **) override;
  void reset_target(double) override;
  double compute_scalar() override;
  void write_restart(FILE *) override;
  void restart(char *buf) override;
  void *extract(const char *, int &) override;

 private:
  double t_start, t_stop, t_period, t_target;
  double energy;
  int nmax, which;
  int tstyle, tvar;
  char *tstr;

  char *id_temp;
  class Compute *temperature;
  int tflag;

  class RanMars *random;

 private:
  double resamplekin(double, double);
  double sumnoises(int);
  double gamdev(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
