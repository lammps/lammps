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

/* ----------------------------------------------------------------------
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#ifndef LMP_FIX_NH_GPU_H
#define LMP_FIX_NH_GPU_H

#include "fix_nh.h"

namespace LAMMPS_NS {

class FixNHGPU : public FixNH {
 public:
  FixNHGPU(class LAMMPS *, int, char **);
  virtual ~FixNHGPU();
  virtual void setup(int vflag);
  void reset_dt();
  virtual void final_integrate();
  virtual double memory_usage();

 protected:
  double *_dtfm;
  int _nlocal3, _nlocal_max, _respa_on;

  virtual void remap();
  virtual void nve_x();
  virtual void nve_v();
  virtual void nh_v_press();
  virtual void nh_v_temp();
};

}    // namespace LAMMPS_NS

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Target temperature for fix nvt/npt/nph cannot be 0.0

Self-explanatory.

E: Invalid fix nvt/npt/nph command for a 2d simulation

Cannot control z dimension in a 2d model.

E: Fix nvt/npt/nph dilate group ID does not exist

Self-explanatory.

E: Invalid fix nvt/npt/nph command pressure settings

If multiple dimensions are coupled, those dimensions must be
specified.

E: Cannot use fix nvt/npt/nph on a non-periodic dimension

When specifying a diagonal pressure component, the dimension must be
periodic.

E: Cannot use fix nvt/npt/nph on a 2nd non-periodic dimension

When specifying an off-diagonal pressure component, the 2nd of the two
dimensions must be periodic.  E.g. if the xy component is specified,
then the y dimension must be periodic.

E: Cannot use fix nvt/npt/nph with yz scaling when z is non-periodic dimension

The 2nd dimension in the barostatted tilt factor must be periodic.

E: Cannot use fix nvt/npt/nph with xz scaling when z is non-periodic dimension

The 2nd dimension in the barostatted tilt factor must be periodic.

E: Cannot use fix nvt/npt/nph with xy scaling when y is non-periodic dimension

The 2nd dimension in the barostatted tilt factor must be periodic.

E: Cannot use fix nvt/npt/nph with both yz dynamics and yz scaling

Self-explanatory.

E: Cannot use fix nvt/npt/nph with both xz dynamics and xz scaling

Self-explanatory.

E: Cannot use fix nvt/npt/nph with both xy dynamics and xy scaling

Self-explanatory.

E: Can not specify Pxy/Pxz/Pyz in fix nvt/npt/nph with non-triclinic box

Only triclinic boxes can be used with off-diagonal pressure components.
See the region prism command for details.

E: Invalid fix nvt/npt/nph pressure settings

Settings for coupled dimensions must be the same.

E: Fix nvt/npt/nph damping parameters must be > 0.0

Self-explanatory.

E: Cannot use fix npt and fix deform on same component of stress tensor

This would be changing the same box dimension twice.

E: Temperature ID for fix nvt/npt does not exist

Self-explanatory.

E: Pressure ID for fix npt/nph does not exist

Self-explanatory.

E: Fix npt/nph has tilted box too far in one step - periodic cell is too far from equilibrium state

Self-explanatory.  The change in the box tilt is too extreme
on a short timescale.

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
