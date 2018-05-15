/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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
  FixRigidNHOMP(class LAMMPS *lmp, int narg, char **args)
    : FixRigidNH(lmp,narg,args) {}
  virtual ~FixRigidNHOMP() {}

  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual void remap();

 protected:
  virtual void compute_forces_and_torques();

 private: // copied from FixRigidOMP
  template <int, int> void set_xv_thr();
  template <int, int> void set_v_thr();
};
}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Target temperature for fix rigid nvt/npt cannot be 0.0

Self-explanatory.

E: Invalid fix rigid npt/nph command for a 2d simulation

Cannot control z dimension in a 2d model.

E: Fix rigid npt/nph dilate group ID does not exist

Self-explanatory.

E: Invalid fix rigid npt/nph command pressure settings

If multiple dimensions are coupled, those dimensions must be
specified.

E: Cannot use fix rigid npt/nph on a non-periodic dimension

When specifying a diagonal pressure component, the dimension must be
periodic.

E: Invalid fix rigid npt/nph pressure settings

Settings for coupled dimensions must be the same.

E: Fix rigid nvt/npt/nph damping parameters must be > 0.0

Self-explanatory.

E: Cannot use fix rigid npt/nph and fix deform on same component of stress tensor

This would be changing the same box dimension twice.

E: Temperature ID for fix rigid npt/nph does not exist

Self-explanatory.

E: Pressure ID for fix rigid npt/nph does not exist

Self-explanatory.

E: Could not find fix_modify temperature ID

The compute ID for computing temperature does not exist.

E: Fix_modify temperature ID does not compute temperature

The compute ID assigned to the fix must compute temperature.

W: Temperature for fix modify is not for group all

The temperature compute is being used with a pressure calculation
which does operate on group all, so this may be inconsistent.

E: Pressure ID for fix modify does not exist

Self-explanatory.

E: Could not find fix_modify pressure ID

The compute ID for computing pressure does not exist.

E: Fix_modify pressure ID does not compute pressure

The compute ID assigned to the fix must compute pressure.

*/
