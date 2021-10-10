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
CommandStyle(temper/npt,TemperNPT);
// clang-format on
#else

#ifndef LMP_TEMPERNPT_H
#define LMP_TEMPERNPT_H

#include "command.h"

namespace LAMMPS_NS {

class TemperNPT : public Command {
 public:
  TemperNPT(class LAMMPS *);
  ~TemperNPT();
  void command(int, char **);

 private:
  int me, me_universe;    // my proc ID in world and universe
  int iworld, nworlds;    // world info
  double boltz;           // copy from output->boltz
  double nktv2p;
  MPI_Comm roots;                       // MPI comm with 1 root proc from each world
  class RanPark *ranswap, *ranboltz;    // RNGs for swapping and Boltz factor
  int nevery;                           // # of timesteps between swaps
  int nswaps;                           // # of tempering swaps to perform
  int seed_swap;                        // 0 = toggle swaps, n = RNG for swap direction
  int seed_boltz;                       // seed for Boltz factor comparison
  int whichfix;                         // index of temperature fix to use
  int fixstyle;                         // what kind of temperature fix is used

  int my_set_temp;     // which set temp I am simulating
  double *set_temp;    // static list of replica set temperatures
  int *temp2world;     // temp2world[i] = world simulating set temp i
  int *world2temp;     // world2temp[i] = temp simulated by world i
  int *world2root;     // world2root[i] = root proc of world i

  void scale_velocities(int, int);
  void print_status();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Must have more than one processor partition to temper

Cannot use the temper command with only one processor partition.  Use
the -partition command-line option.

E: temper/npt command before simulation box is defined

The temper/npt command cannot be used before a read_data, read_restart,
or create_box command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Tempering fix ID is not defined

The fix ID specified by the temper/npt command does not exist.

E: Invalid frequency in temper/npt command

Nevery must be > 0.

E: Non integer # of swaps in temper/npt command

Swap frequency in temper/npt command must evenly divide the total
# of timesteps.

E: Tempering temperature fix is not valid

The fix specified by the temper command is not one that controls
temperature and pressure (npt).

E: Too many timesteps

The cummulative timesteps must fit in a 64-bit integer.

E: Tempering could not find thermo_pe compute

This compute is created by the thermo command.  It must have been
explicitly deleted by a uncompute command.

*/
