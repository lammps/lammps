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

#ifdef ANGLE_CLASS
// clang-format off
AngleStyle(cosine/buck6d, AngleCosineBuck6d);
// clang-format on
#else

#ifndef LMP_ANGLE_COSINE_BUCK6D_H
#define LMP_ANGLE_COSINE_BUCK6D_H

#include "angle.h"

namespace LAMMPS_NS {

class AngleCosineBuck6d : public Angle {
 public:
  AngleCosineBuck6d(class LAMMPS *);
  ~AngleCosineBuck6d() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double equilibrium_angle(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;
  double single(int, int, int, int) override;

 protected:
  double *k, *th0;
  double *eps, *d0;
  double **buck6d1, **buck6d2, **buck6d3, **buck6d4, **cut_ljsq;
  double **c0, **c1, **c2, **c3, **c4, **c5, **rsmooth_sq, **offset;
  int *multiplicity;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
