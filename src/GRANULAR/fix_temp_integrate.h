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
FixStyle(temp/integrate,FixTempIntegrate);
// clang-format on
#else

#ifndef LMP_FIX_TEMP_INTEGRATE_H
#define LMP_FIX_TEMP_INTEGRATE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTempIntegrate : public Fix {
 public:
  FixTempIntegrate(class LAMMPS *, int, char **);

  int setmask() override;
  void init() override;
  void final_integrate() override;
  void final_integrate_respa(int, int) override;
  void reset_dt() override;

 protected:
  double dt;
  double cp, *cp_type;
  int cp_style;

  int mass_require;

  double calc_cp(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
