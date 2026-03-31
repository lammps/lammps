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

#ifdef ANGLE_CLASS
// clang-format off
AngleStyle(class2xe,AngleClass2xe);
// clang-format on
#else

#ifndef LMP_ANGLE_CLASS2XE_H
#define LMP_ANGLE_CLASS2XE_H

#include "angle.h"

namespace LAMMPS_NS {

class AngleClass2xe : public Angle {
 public:
  AngleClass2xe(class LAMMPS *);
  ~AngleClass2xe() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  double equilibrium_angle(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;
  double single(int, int, int, int) override;

 protected:
  double *theta0, *k2, *k3, *k4;
  double *bb_d0, *bb_alpha, *bb_r1, *bb_r2;
  double *ba_d1, *ba_d2, *ba_alpha1, *ba_alpha2, *ba_r1, *ba_r2;
  int *setflag_a, *setflag_bb, *setflag_ba;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
