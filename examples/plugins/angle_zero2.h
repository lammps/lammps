/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_ANGLE_ZERO2_H
#define LMP_ANGLE_ZERO2_H

#include "angle.h"

namespace LAMMPS_NS {

class AngleZero2 : public Angle {
 public:
  AngleZero2(class LAMMPS *);
  ~AngleZero2() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void settings(int, char **) override;

  double equilibrium_angle(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;

  double single(int, int, int, int) override;

 protected:
  double *theta0;
  int coeffflag;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

UNDOCUMENTED

E: Incorrect args for angle coefficients

Self-explanatory.  Check the input script or data file.

*/
