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
FixStyle(gravity,FixGravity);
// clang-format on
#else

#ifndef LMP_FIX_GRAVITY_H
#define LMP_FIX_GRAVITY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGravity : public Fix {
  friend class FixPour;

 public:
  FixGravity(class LAMMPS *, int, char **);
  ~FixGravity() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  double compute_scalar() override;
  void *extract(const char *, int &) override;
  enum { CONSTANT, EQUAL };

 protected:
  int style, disable;
  double magnitude;
  double vert, phi, theta;
  double xdir, ydir, zdir;
  double xgrav, ygrav, zgrav, xacc, yacc, zacc;
  int ilevel_respa;
  int time_origin;
  double gvec[3];

  int eflag;
  double egrav, egrav_all;

  int varflag;
  int mstyle, vstyle, pstyle, tstyle, xstyle, ystyle, zstyle;
  int mvar, vvar, pvar, tvar, xvar, yvar, zvar;
  char *mstr, *vstr, *pstr, *tstr, *xstr, *ystr, *zstr;

  void set_acceleration();
};

}    // namespace LAMMPS_NS

#endif
#endif
