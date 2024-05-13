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
FixStyle(rigid/omp,FixRigidOMP);
// clang-format on
#else

#ifndef LMP_FIX_RIGID_OMP_H
#define LMP_FIX_RIGID_OMP_H

#include "fix_rigid.h"

namespace LAMMPS_NS {

class FixRigidOMP : public FixRigid {
 public:
  FixRigidOMP(class LAMMPS *lmp, int narg, char **args) : FixRigid(lmp, narg, args) {}

  void initial_integrate(int) override;
  void final_integrate() override;

 protected:
  void compute_forces_and_torques() override;

 private:
  template <int, int, int> void set_xv_thr();
  template <int, int, int> void set_v_thr();
};

}    // namespace LAMMPS_NS

#endif
#endif
