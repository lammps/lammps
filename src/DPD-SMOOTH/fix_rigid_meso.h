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

#ifdef FIX_CLASS
// clang-format off
FixStyle(rigid/meso,FixRigidMeso);
// clang-format on
#else

#ifndef LMP_FIX_RIGID_MESO_H
#define LMP_FIX_RIGID_MESO_H

#include "fix_rigid.h"

namespace LAMMPS_NS {

class FixRigidMeso : public FixRigid {
 public:
  FixRigidMeso(class LAMMPS *, int, char **);
  ~FixRigidMeso() override;
  int setmask() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void final_integrate() override;
  double compute_scalar() override { return 0.0; }
  double compute_array(int, int) override;

 protected:
  void set_xv();
  void set_v();
  double **conjqm;    // conjugate quaternion momentum
};

}    // namespace LAMMPS_NS

#endif
#endif
