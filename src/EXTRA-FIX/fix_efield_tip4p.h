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
FixStyle(efield/tip4p,FixEfieldTIP4P);
// clang-format on
#else

#ifndef LMP_FIX_EFIELD_TIP4P_H
#define LMP_FIX_EFIELD_TIP4P_H

#include "fix_efield.h"

namespace LAMMPS_NS {

class FixEfieldTIP4P : public FixEfield {

 public:
  FixEfieldTIP4P(class LAMMPS *, int, char **);

  void init() override;
  void post_force(int) override;

 protected:
  double alpha;
  int typeO, typeH;    // atom types for TIP4P molecule
  void find_M(double *, double *, double *, double *);
};
}    // namespace LAMMPS_NS

#endif
#endif
