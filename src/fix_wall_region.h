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
FixStyle(wall/region,FixWallRegion);
// clang-format on
#else

#ifndef LMP_FIX_WALL_REGION_H
#define LMP_FIX_WALL_REGION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallRegion : public Fix {
 public:
  FixWallRegion(class LAMMPS *, int, char **);
  ~FixWallRegion() override;
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
  int style;
  double epsilon, sigma, cutoff;
  double alpha;
  int eflag;
  double ewall[4], ewall_all[4];
  int ilevel_respa;
  char *idregion;
  class Region *region;

  double coeff1, coeff2, coeff3, coeff4, offset;
  double coeff5, coeff6, coeff7;
  double eng, fwall;

  void lj93(double);
  void lj126(double);
  void lj1043(double);
  void morse(double);
  void colloid(double, double);
  void harmonic(double);
};

}    // namespace LAMMPS_NS

#endif
#endif
