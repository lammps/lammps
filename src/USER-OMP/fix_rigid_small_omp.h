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
FixStyle(rigid/small/omp,FixRigidSmallOMP);
// clang-format on
#else

#ifndef LMP_FIX_RIGID_SMALL_OMP_H
#define LMP_FIX_RIGID_SMALL_OMP_H

#include "fix_rigid_small.h"

namespace LAMMPS_NS {

class FixRigidSmallOMP : public FixRigidSmall {
 public:
  FixRigidSmallOMP(class LAMMPS *lmp, int narg, char **args) : FixRigidSmall(lmp, narg, args){};
  virtual ~FixRigidSmallOMP(){};

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

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix rigid molecule requires atom attribute molecule

Self-explanatory.

E: Could not find fix rigid group ID

A group ID used in the fix rigid command does not exist.

E: One or more atoms belong to multiple rigid bodies

Two or more rigid bodies defined by the fix rigid command cannot
contain the same atom.

E: No rigid bodies defined

The fix specification did not end up defining any rigid bodies.

E: Fix rigid z force cannot be on for 2d simulation

Self-explanatory.

E: Fix rigid xy torque cannot be on for 2d simulation

Self-explanatory.

E: Fix rigid langevin period must be > 0.0

Self-explanatory.

E: One or zero atoms in rigid body

Any rigid body defined by the fix rigid command must contain 2 or more
atoms.

W: More than one fix rigid

It is not efficient to use fix rigid more than once.

E: Rigid fix must come before NPT/NPH fix

NPT/NPH fix must be defined in input script after all rigid fixes,
else the rigid fix contribution to the pressure virial is
incorrect.

W: Computing temperature of portions of rigid bodies

The group defined by the temperature compute does not encompass all
the atoms in one or more rigid bodies, so the change in
degrees-of-freedom for the atoms in those partial rigid bodies will
not be accounted for.

E: Fix rigid atom has non-zero image flag in a non-periodic dimension

You cannot set image flags for non-periodic dimensions.

E: Insufficient Jacobi rotations for rigid body

Eigensolve for rigid body was not sufficiently accurate.

E: Fix rigid: Bad principal moments

The principal moments of inertia computed for a rigid body
are not within the required tolerances.

E: Cannot open fix rigid infile %s

UNDOCUMENTED

E: Unexpected end of fix rigid file

UNDOCUMENTED

E: Incorrect rigid body format in fix rigid file

UNDOCUMENTED

E: Invalid rigid body ID in fix rigid file

UNDOCUMENTED

*/
