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

#ifdef COMMAND_CLASS

CommandStyle(velocity,Velocity)

#else

#ifndef LMP_VELOCITY_H
#define LMP_VELOCITY_H

#include "pointers.h"

namespace LAMMPS_NS {

class Velocity : protected Pointers {
 public:
  Velocity(class LAMMPS *);
  void command(int, char **);
  void init_external(const char *);
  void options(int, char **);
  void create(double, int);

 private:
  int igroup,groupbit;
  int style;
  int dist_flag,sum_flag,momentum_flag,rotation_flag;
  int bias_flag,loop_flag,scale_flag,rfix;
  double xscale,yscale,zscale;
  class Compute *temperature;

  void set(int, char **);
  void scale(int, char **);
  void ramp(int, char **);
  void zero(int, char **);

  void rescale(double, double);
  void zero_momentum();
  void zero_rotation();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Velocity command before simulation box is defined

The velocity command cannot be used before a read_data, read_restart,
or create_box command.

E: Velocity command with no atoms existing

A velocity command has been used, but no atoms yet exist.

E: Could not find velocity group ID

A group ID used in the velocity command does not exist.

W: Mismatch between velocity and compute groups

The temperature computation used by the velocity command will not be
on the same group of atoms that velocities are being set for.

E: Too big a problem to use velocity create loop all

The system size must fit in a 32-bit integer to use this option.

E: Cannot use velocity create loop all unless atoms have IDs

Atoms in the simulation to do not have IDs, so this style
of velocity creation cannot be performed.

E: Atom IDs must be consecutive for velocity create loop all

Self-explanatory.

E: Variable name for velocity set does not exist

Self-explanatory.

E: Variable for velocity set is invalid style

Only atom-style variables can be used.

E: Cannot set non-zero z velocity for 2d simulation

Self-explanatory.

E: Cannot set variable z velocity for 2d simulation

Self-explanatory.

E: Velocity ramp in z for a 2d problem

Self-explanatory.

E: Velocity rigid used with non-rigid fix-ID

Self-explanatory.

E: Attempting to rescale a 0.0 temperature

Cannot rescale a temperature that is already 0.0.

E: Cannot zero momentum of no atoms

Self-explanatory.

E: Could not find velocity temperature ID

The compute ID needed by the velocity command to compute temperature
does not exist.

E: Velocity temperature ID does not compute temperature

The compute ID given to the velocity command must compute
temperature.

E: Fix ID for velocity does not exist

Self-explanatory.

*/
