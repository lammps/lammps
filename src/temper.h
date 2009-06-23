/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef TEMPER_H
#define TEMPER_H

#include "pointers.h"

namespace LAMMPS_NS {

class Temper : protected Pointers {
 public:
  Temper(class LAMMPS *);
  ~Temper();
  void command(int, char **);

 private:
  int me,me_universe;          // my proc ID in world and universe
  int iworld,nworlds;          // world info
  double boltz;                // copy from output->boltz
  MPI_Comm roots;              // Comm with 1 root proc from each world
  class RanPark *ranswap,*ranboltz;  // RNGs for swapping and Boltz factor
  int nevery;                  // # of timesteps between swaps
  int nswaps;                  // # of tempering swaps to perform
  int seed_swap;               // 0 = toggle swaps, n = RNG for swap direction
  int seed_boltz;              // seed for Boltz factor comparison
  int whichfix;                // index of temperature fix to use
  int fixstyle;                // what kind of temperature fix is used

  int my_set_temp;             // which set temp I am simulating
  double *set_temp;            // static list of replica set temperatures
  int *temp2world;             // temp2world[i] = world simulating set temp i
  int *world2temp;             // world2temp[i] = temp simulated by world i
  int *world2root;             // world2root[i] = root proc of world i

  void scale_velocities(int, int);
  void print_status();
};

}

#endif
