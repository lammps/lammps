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

#ifndef TIMER_H
#define TIMER_H

#include "pointers.h"

enum{TIME_LOOP,TIME_PAIR,TIME_BOND,TIME_KSPACE,TIME_NEIGHBOR,
     TIME_COMM,TIME_OUTPUT,TIME_N};

namespace LAMMPS_NS {

class Timer : protected Pointers {
 public:
  double *array;

  Timer(class LAMMPS *);
  ~Timer();
  void init();
  void stamp();
  void stamp(int);
  void barrier_start(int);
  void barrier_stop(int);
  double elapsed(int);

 private:
  double previous_time;
};

}

#endif
