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
FixStyle(langevin/drude,FixLangevinDrude);
// clang-format on
#else

#ifndef LMP_FIX_LANGEVIN_DRUDE_H
#define LMP_FIX_LANGEVIN_DRUDE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLangevinDrude : public Fix {
 public:
  FixLangevinDrude(class LAMMPS *, int, char **);
  ~FixLangevinDrude() override;
  int setmask() override;
  void init() override;
  void setup(int vflag) override;
  void post_force(int vflag) override;
  void reset_target(double) override;
  void *extract(const char *, int &) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  int modify_param(int, char **) override;

 protected:
  double t_start_core, t_period_core, t_target_core;
  double t_start_drude, t_period_drude, t_target_drude;
  int tstyle_core, tstyle_drude;
  int tvar_core, tvar_drude;
  char *tstr_core, *tstr_drude;
  double energy;
  int tflag;

  class RanMars *random_core, *random_drude;
  int zero;
  bigint ncore;
  class FixDrude *fix_drude;
  class Compute *temperature;
  char *id_temp;
};

}    // namespace LAMMPS_NS

#endif
#endif
