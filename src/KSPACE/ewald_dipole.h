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

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(ewald/dipole,EwaldDipole);
// clang-format on
#else

#ifndef LMP_EWALD_DIPOLE_H
#define LMP_EWALD_DIPOLE_H

#include "ewald.h"

namespace LAMMPS_NS {

class EwaldDipole : public Ewald {
 public:
  EwaldDipole(class LAMMPS *);
  ~EwaldDipole() override;
  void init() override;
  void setup() override;
  void compute(int, int) override;

 protected:
  double musum, musqsum, mu2;
  double **tk;    // field for torque
  double **vc;    // virial per k

  void musum_musq();
  double rms_dipole(int, double, bigint);
  void eik_dot_r() override;
  void slabcorr();
  double NewtonSolve(double, double, bigint, double, double);
  double f(double, double, bigint, double, double);
  double derivf(double, double, bigint, double, double);
};

}    // namespace LAMMPS_NS

#endif
#endif
