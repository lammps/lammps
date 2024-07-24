/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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
  enum ttype {
    RESET = -2,
    START = -1,
    TOTAL = 0,
    PAIR,
    BOND,
    KSPACE,
    NEIGH,
    COMM,
    MODIFY,
    OUTPUT,
    SYNC,
    ALL,
    DEPHASE,
    DYNAMICS,
    QUENCH,
    NEB,
    REPCOMM,
    REPOUT,
    NUM_TIMER
  };
  enum tlevel { OFF = 0, LOOP, NORMAL, FULL };

  Timer(class LAMMPS *);

  void init();

  // inline function to reduce overhead if we want no detailed timings

  void stamp(enum ttype which = START)
  {
    if (_level > LOOP) _stamp(which);
  }

  void barrier_start();
  void barrier_stop();

  // accessor methods for supported level of detail

  bool has_loop() const { return (_level >= LOOP); }
  bool has_normal() const { return (_level >= NORMAL); }
  bool has_full() const { return (_level >= FULL); }
  bool has_sync() const { return (_sync != OFF); }
  bool has_timeout() const { return (_timeout >= 0.0); }

  // flag if wallclock time is expired
  bool is_timeout() const { return (_timeout == 0.0); }

  double elapsed(enum ttype);
  double cpu(enum ttype);

  double get_cpu(enum ttype which) const { return cpu_array[which]; };
  double get_wall(enum ttype which) const { return wall_array[which]; };

  void set_wall(enum ttype, double);

  // initialize timeout timer
  void init_timeout();

  // trigger enforced timeout
  void force_timeout() { _timeout = 0.0; }

  // restore original timeout setting after enforce timeout
  void reset_timeout() { _timeout = _s_timeout; }

  // get remaining time in seconds. 0.0 if inactive, negative if expired
  double get_timeout_remain();

  // print timeout message
  void print_timeout(FILE *);

  // check for timeout. inline wrapper around internal
  // function to reduce overhead in case there is no check.
  bool check_timeout(int step)
  {
    if (_timeout == 0.0) return true;
    if (_nextcheck != step)
      return false;
    else
      return _check_timeout();
  }

  void modify_params(int, char **);

 private:
  double cpu_array[NUM_TIMER];
  double wall_array[NUM_TIMER];
  double previous_cpu;
  double previous_wall;
  double timeout_start;
  double _timeout;      // max allowed wall time in seconds. infinity if negative
  double _s_timeout;    // copy of timeout for restoring after a forced timeout
  int _level;           // level of detail: off=0,loop=1,normal=2,full=3
  int _sync;            // if nonzero, synchronize tasks before setting the timer
  int _checkfreq;       // frequency of timeout checking
  int _nextcheck;       // loop number of next timeout check

  // update one specific timer array
  void _stamp(enum ttype);

  // check for timeout
  bool _check_timeout();
};

}    // namespace LAMMPS_NS

#endif
