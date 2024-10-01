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
FixStyle(add/heat,FixAddHeat);
// clang-format on
#else

#ifndef LMP_FIX_ADD_HEAT_H
#define LMP_FIX_ADD_HEAT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAddHeat : public Fix {
 public:
  FixAddHeat(class LAMMPS *, int, char **);
  ~FixAddHeat() override;
  int setmask() override;
  void init() override;
  void post_force(int) override;

 protected:
  double value, prefactor;
  int var, vstyle, maxatom, style, overwrite_flag;
  char *varstr;
  double *vatom;
};

}    // namespace LAMMPS_NS

#endif
#endif
