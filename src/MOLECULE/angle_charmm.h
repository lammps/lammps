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
AngleStyle(charmm,AngleCharmm);
// clang-format on
#else

#ifndef LMP_ANGLE_CHARMM_H
#define LMP_ANGLE_CHARMM_H

#include "angle.h"

namespace LAMMPS_NS {

class AngleCharmm : public Angle {
 public:
  AngleCharmm(class LAMMPS *);
  virtual ~AngleCharmm();
  virtual void compute(int, int);
  virtual void coeff(int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  virtual void read_restart(FILE *);
  void write_data(FILE *);
  double single(int, int, int, int);

 protected:
  double *k, *theta0, *k_ub, *r_ub;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for angle coefficients

Self-explanatory.  Check the input script or data file.

*/
