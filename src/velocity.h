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

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(velocity,Velocity);
// clang-format on
#else

#ifndef LMP_VELOCITY_H
#define LMP_VELOCITY_H

#include "command.h"

namespace LAMMPS_NS {

class Velocity : public Command {
 public:
  Velocity(class LAMMPS *);
  void command(int, char **) override;
  void init_external(const char *);
  void options(int, char **);
  void create(double, int);

 private:
  int igroup, groupbit;
  int style;
  int dist_flag, sum_flag, momentum_flag, rotation_flag;
  int bias_flag, loop_flag, scale_flag;
  double xscale, yscale, zscale;
  class Fix *rigid_fix;
  class Compute *temperature;

  void set(int, char **);
  void scale(int, char **);
  void ramp(int, char **);
  void zero(int, char **);

  void rescale(double, double);
  void zero_momentum();
  void zero_rotation();
};

}    // namespace LAMMPS_NS

#endif
#endif
