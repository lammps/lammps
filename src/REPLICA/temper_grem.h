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
  ~TemperGrem() override;
  void command(int, char **) override;

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
