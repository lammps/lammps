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
FixStyle(aveforce,FixAveForce);
// clang-format on
#else

#ifndef LMP_FIX_AVEFORCE_H
#define LMP_FIX_AVEFORCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAveForce : public Fix {
 public:
  FixAveForce(class LAMMPS *, int, char **);
  ~FixAveForce() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double compute_vector(int) override;

 private:
  double xvalue, yvalue, zvalue;
  int varflag;
  char *xstr, *ystr, *zstr;
  char *idregion;
  class Region *region;
  int xvar, yvar, zvar, xstyle, ystyle, zstyle;
  double foriginal_all[4];
  int nlevels_respa, ilevel_respa;
};

}    // namespace LAMMPS_NS

#endif
#endif
