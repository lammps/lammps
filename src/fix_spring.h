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
FixStyle(spring,FixSpring);
// clang-format on
#else

#ifndef LMP_FIX_SPRING_H
#define LMP_FIX_SPRING_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSpring : public Fix {
 public:
  FixSpring(class LAMMPS *, int, char **);
  ~FixSpring() override;
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
  double xc, yc, zc, r0;
  double k_spring;
  int xflag, yflag, zflag;
  int styleflag;
  char *group2;
  int igroup2, group2bit;
  double masstotal, masstotal2;
  int ilevel_respa;
  double espring, ftotal[4];

  void spring_tether();
  void spring_couple();
};

}    // namespace LAMMPS_NS

#endif
#endif
