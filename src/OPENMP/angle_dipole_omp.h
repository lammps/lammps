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

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef ANGLE_CLASS
// clang-format off
AngleStyle(dipole/omp,AngleDipoleOMP);
// clang-format on
#else

#ifndef LMP_ANGLE_DIPOLE_OMP_H
#define LMP_ANGLE_DIPOLE_OMP_H

#include "angle_dipole.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class AngleDipoleOMP : public AngleDipole, public ThrOMP {

 public:
  AngleDipoleOMP(class LAMMPS *lmp);
  void compute(int, int) override;

 private:
  template <int EFLAG> void eval(int ifrom, int ito, ThrData *const thr);
};

}    // namespace LAMMPS_NS

#endif
#endif
