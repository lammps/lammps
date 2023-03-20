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
FixStyle(efield,FixEfield);
// clang-format on
#else

#ifndef LMP_FIX_EFIELD_H
#define LMP_FIX_EFIELD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEfield : public Fix {
  friend class FixQEqReaxFF;

 public:
  FixEfield(class LAMMPS *, int, char **);
  ~FixEfield() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double memory_usage() override;
  double compute_scalar() override;
  double compute_vector(int) override;

  enum { NONE, CONSTANT, EQUAL, ATOM };

 protected:
  double ex, ey, ez;
  int varflag;
  char *xstr, *ystr, *zstr, *estr;
  char *idregion;
  class Region *region;
  int xvar, yvar, zvar, evar, xstyle, ystyle, zstyle, estyle;
  int ilevel_respa;
  double qe2f;
  int qflag, muflag;

  int maxatom, maxatom_energy;
  double **efield;

  int force_flag;
  double fsum[4], fsum_all[4];
};
}    // namespace LAMMPS_NS
#endif
#endif
