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
FixStyle(append/atoms,FixAppendAtoms);
// clang-format on
#else

#ifndef FIX_APPEND_ATOMS_H
#define FIX_APPEND_ATOMS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAppendAtoms : public Fix {
 public:
  FixAppendAtoms(class LAMMPS *, int, char **);
  ~FixAppendAtoms() override;
  int setmask() override;
  void setup(int) override;
  void pre_exchange() override;
  void initial_integrate(int) override;
  void post_force(int) override;

 private:
  int get_spatial();
  int spatflag, xloflag, xhiflag, yloflag, yhiflag, zloflag, zhiflag;
  int ranflag, tempflag, xseed, tseed;
  double ranx, rany, ranz, t_target, t_period, t_extent;
  class RanMars *randomx;
  class RanMars *randomt;
  int scaleflag, freq;
  int nbasis;
  int *basistype;
  int advance, advance_sum;
  double size, spatlead;
  char *spatialid;
  double tfactor;
  double *gfactor1, *gfactor2;
};

}    // namespace LAMMPS_NS

#endif
#endif
