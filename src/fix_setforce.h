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
FixStyle(setforce,FixSetForce);
// clang-format on
#else

#ifndef LMP_FIX_SET_FORCE_H
#define LMP_FIX_SET_FORCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSetForce : public Fix {
 public:
  FixSetForce(class LAMMPS *, int, char **);
  ~FixSetForce() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double compute_vector(int) override;

  double memory_usage() override;

 protected:
  double xvalue, yvalue, zvalue;
  int varflag;
  char *xstr, *ystr, *zstr;
  char *idregion;
  class Region *region;
  int xvar, yvar, zvar, xstyle, ystyle, zstyle;
  double foriginal[3], foriginal_all[3], foriginal_saved[3];
  int force_flag;
  int nlevels_respa, ilevel_respa;

  int maxatom;
  double **sforce;
};

}    // namespace LAMMPS_NS

#endif
#endif
