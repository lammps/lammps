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
AngleStyle(spica,AngleSPICA);
AngleStyle(sdk,AngleSPICA);
// clang-format on
#else

#ifndef LMP_ANGLE_SPICA_H
#define LMP_ANGLE_SPICA_H

#include "angle.h"

namespace LAMMPS_NS {

class AngleSPICA : public Angle {
 public:
  AngleSPICA(class LAMMPS *);
  ~AngleSPICA() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double equilibrium_angle(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;
  double single(int, int, int, int) override;

 protected:
  double *k, *theta0;

  // scaling factor for repulsive 1-3 interaction
  double *repscale;
  // parameters from SPICA pair style
  int **lj_type;
  double **lj1, **lj2, **lj3, **lj4;
  double **rminsq, **emin;

  int repflag;    // 1 if we have to handle 1-3 repulsion

  void ev_tally13(int, int, int, int, double, double, double, double, double);

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
