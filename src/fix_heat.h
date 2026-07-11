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
FixStyle(heat,FixHeat);
// clang-format on
#else

#ifndef LMP_FIX_HEAT_H
#define LMP_FIX_HEAT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixHeat : public Fix {
 public:
  FixHeat(class LAMMPS *, int, char **);
  ~FixHeat() override;
  int setmask() override;
  void init() override;
  void end_of_step() override;
  double compute_scalar() override;
  double memory_usage() override;

 private:
  double heat_input;
  double masstotal;
  double scale;
  char *idregion;
  class Region *region;
  char *hstr;
  int hstyle, hvar;

  int maxatom;
  double *vheat;
  double *vscale;
};

}    // namespace LAMMPS_NS

#endif
#endif
