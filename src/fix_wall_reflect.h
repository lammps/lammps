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
FixStyle(wall/reflect,FixWallReflect);
// clang-format on
#else

#ifndef LMP_FIX_WALL_REFLECT_H
#define LMP_FIX_WALL_REFLECT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallReflect : public Fix {
 public:
  enum { XLO = 0, XHI = 1, YLO = 2, YHI = 3, ZLO = 4, ZHI = 5 };
  enum { NONE = 0, EDGE, CONSTANT, VARIABLE };

  FixWallReflect(class LAMMPS *, int, char **);
  ~FixWallReflect() override;
  int setmask() override;
  void init() override;
  void post_integrate() override;

 protected:
  int nwall;
  int wallwhich[6], wallstyle[6];
  double coord0[6];
  char *varstr[6];
  int varindex[6];
  int varflag;
  double xscale, yscale, zscale;

  virtual void wall_particle(int m, int which, double coord);
};

}    // namespace LAMMPS_NS

#endif
#endif
