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
CommandStyle(temper/grem,TemperGrem);
// clang-format on
#else

#ifndef LMP_TEMPER_GREM_H
#define LMP_TEMPER_GREM_H

#include "command.h"

namespace LAMMPS_NS {

class TemperGrem : public Command {
 public:
  TemperGrem(class LAMMPS *);
  ~TemperGrem();
  void command(int, char **);

 private:
  int me, me_universe;                  // my proc ID in world and universe
  int iworld, nworlds;                  // world info
  double boltz;                         // copy from output->boltz
  MPI_Comm roots;                       // MPI comm with 1 root proc from each world
  class RanPark *ranswap, *ranboltz;    // RNGs for swapping and Boltz factor
  int nevery;                           // # of timesteps between swaps
  int nswaps;                           // # of tempering swaps to perform
  int seed_swap;                        // 0 = toggle swaps, n = RNG for swap direction
  int seed_boltz;                       // seed for Boltz factor comparison
  int whichfix;                         // index of temperature fix to use
  int fixstyle;                         // what kind of temperature fix is used

  int my_set_lambda;     // which set lambda I am simulating
  double *set_lambda;    // static list of replica set lambdas
  int *lambda2world;     // lambda2world[i] = world simulating set lambda i
  int *world2lambda;     // world2lambda[i] = lambda simulated by world i
  int *world2root;       // world2root[i] = root proc of world i

  void print_status();

  class FixGrem *fix_grem;

 protected:
  char *id_nh;
  int pressflag;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Must have more than one processor partition to grem

Cannot use the grem command with only one processor partition.  Use
the -partition command-line option.

E: Grem command before simulation box is defined

The grem command cannot be used before a read_data, read_restart, or
create_box command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Tempering fix ID is not defined

The fix ID specified by the grem command does not exist.

E: Invalid frequency in grem command

Nevery must be > 0.

E: Non integer # of swaps in grem command

Swap frequency in grem command must evenly divide the total # of
timesteps.

E: Grem temperature fix is not valid

The fix specified by the grem command is not one that controls
temperature (nvt or npt).

E: Too many timesteps

The cumulative timesteps must fit in a 64-bit integer.

E: Grem could not find thermo_pe compute

This compute is created by the thermo command.  It must have been
explicitly deleted by a uncompute command.

*/
