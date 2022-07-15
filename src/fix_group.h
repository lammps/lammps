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
FixStyle(GROUP,FixGroup);
// clang-format on
#else

#ifndef LMP_FIX_GROUP_H
#define LMP_FIX_GROUP_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGroup : public Fix {
 public:
  FixGroup(class LAMMPS *, int, char **);
  ~FixGroup() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  void *extract(const char *, int &) override;

 private:
  int gbit, gbitinverse;
  int regionflag, varflag, propflag, proptype;
  int ivar, iprop;
  char *idregion, *idvar, *idprop;
  class Region *region;

  int nlevels_respa;

  void set_group();
};

}    // namespace LAMMPS_NS

#endif
#endif
