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

#ifndef LMP_FIX_RIGID_NH_OMP_H
#define LMP_FIX_RIGID_NH_OMP_H

#include "fix_rigid_nh.h"

namespace LAMMPS_NS {

class FixRigidNHOMP : public FixRigidNH {
 public:
  FixRigidNHOMP(class LAMMPS *lmp, int narg, char **args) : FixRigidNH(lmp, narg, args) {}

  void initial_integrate(int) override;
  void final_integrate() override;
  void remap() override;

 protected:
  void compute_forces_and_torques() override;

 private:    // copied from FixRigidOMP
  template <int, int, int> void set_xv_thr();
  template <int, int, int> void set_v_thr();
};
}    // namespace LAMMPS_NS

#endif
