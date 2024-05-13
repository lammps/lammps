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
FixStyle(addtorque,FixAddTorque);
// clang-format on
#else

#ifndef LMP_FIX_ADDTORQUE_H
#define LMP_FIX_ADDTORQUE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAddTorque : public Fix {
 public:
  FixAddTorque(class LAMMPS *, int, char **);
  ~FixAddTorque() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double compute_scalar() override;
  double compute_vector(int) override;

 private:
  double xvalue, yvalue, zvalue;
  int varflag;
  char *xstr, *ystr, *zstr;
  int xvar, yvar, zvar, xstyle, ystyle, zstyle;
  double foriginal[4], foriginal_all[4];
  int force_flag;
  int ilevel_respa;
};

}    // namespace LAMMPS_NS

#endif
#endif
