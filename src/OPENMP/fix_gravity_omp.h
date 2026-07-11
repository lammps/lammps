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

#ifdef FIX_CLASS
// clang-format off
FixStyle(gravity/omp,FixGravityOMP);
// clang-format on
#else

#ifndef LMP_FIX_GRAVITY_OMP_H
#define LMP_FIX_GRAVITY_OMP_H

#include "fix_gravity.h"

namespace LAMMPS_NS {

class FixGravityOMP : public FixGravity {

 public:
  FixGravityOMP(class LAMMPS *, int, char **);
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
