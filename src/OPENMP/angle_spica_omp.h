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
AngleStyle(spica/omp,AngleSPICAOMP);
AngleStyle(sdk/omp,AngleSPICAOMP);
// clang-format on
#else

#ifndef LMP_ANGLE_SPICA_OMP_H
#define LMP_ANGLE_SPICA_OMP_H

#include "angle_spica.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class AngleSPICAOMP : public AngleSPICA, public ThrOMP {

 public:
  AngleSPICAOMP(class LAMMPS *lmp);
  void compute(int, int) override;

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_BOND>
  void eval(int ifrom, int ito, ThrData *const thr);
};

}    // namespace LAMMPS_NS

#endif
#endif
