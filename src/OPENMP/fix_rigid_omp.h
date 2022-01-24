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
  ~FixRigidOMP() {}

  virtual void initial_integrate(int);
  virtual void final_integrate();

 protected:
  virtual void compute_forces_and_torques();

 private:
  template <int, int> void set_xv_thr();
  template <int, int> void set_v_thr();
};

}    // namespace LAMMPS_NS

#endif
#endif
