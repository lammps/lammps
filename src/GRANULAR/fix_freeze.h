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
FixStyle(freeze,FixFreeze);
// clang-format on
#else

#ifndef LMP_FIX_FREEZE_H
#define LMP_FIX_FREEZE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixFreeze : public Fix {
 public:
  FixFreeze(class LAMMPS *, int, char **);
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  double compute_vector(int) override;

 protected:
  int force_flag;
  double foriginal[3], foriginal_all[3];
};

}    // namespace LAMMPS_NS

#endif
#endif
