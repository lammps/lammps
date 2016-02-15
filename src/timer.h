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

#ifndef LMP_TIMER_H
#define LMP_TIMER_H

#include "pointers.h"


namespace LAMMPS_NS {

class Timer : protected Pointers {
 public:

  enum ttype  {RESET=-2,START=-1,TOTAL=0,PAIR,BOND,KSPACE,NEIGH,COMM,
               MODIFY,OUTPUT,SYNC,ALL,DEPHASE,DYNAMICS,QUENCH,NEB,REPCOMM,
               REPOUT,NUM_TIMER};
  enum tlevel {OFF=0,LOOP,NORMAL,FULL};

  Timer(class LAMMPS *);
  ~Timer() {};
  void init();

  // inline function to reduce overhead if we want no detailed timings

  void stamp(enum ttype which=START) {
    if (_level > LOOP) _stamp(which);
  }

  void barrier_start();
  void barrier_stop();

  // accessor methods for supported level of detail

  bool has_loop()   const { return (_level >= LOOP); }
  bool has_normal() const { return (_level >= NORMAL); }
  bool has_full()   const { return (_level >= FULL); }
  bool has_sync()   const { return (_sync  != OFF); }

  double elapsed(enum ttype);
  double cpu(enum ttype);

  double get_cpu(enum ttype which) const {
    return cpu_array[which]; };
  double get_wall(enum ttype which) const {
    return wall_array[which]; };

  void set_wall(enum ttype, double);


  void modify_params(int, char **);

 private:
  double cpu_array[NUM_TIMER];
  double wall_array[NUM_TIMER];
  double previous_cpu;
  double previous_wall;
  int _level;  // level of detail: off=0,loop=1,normal=2,full=3
  int _sync;   // if nonzero, synchronize tasks before setting the timer

  // update requested timer array
  void _stamp(enum ttype);
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
