/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Philipp Kloza (University of Cambridge)
                        pak37@cam.ac.uk
------------------------------------------------------------------------- */

#ifdef ANGLE_CLASS
// clang-format off
AngleStyle(mesocnt, AngleMesoCNT);
// clang-format on
#else

#ifndef LMP_ANGLE_MESOCNT_H
#define LMP_ANGLE_MESOCNT_H

#include "angle.h"

namespace LAMMPS_NS {

class AngleMesoCNT : public Angle {
 public:
  AngleMesoCNT(class LAMMPS *);
  ~AngleMesoCNT() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double equilibrium_angle(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;
  double single(int, int, int, int) override;

 protected:
  int *buckling;
  double *kh, *kb, *thetab, *theta0;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
