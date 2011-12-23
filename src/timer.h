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
  enum ttype {LOOP=0,PAIR,BOND,KSPACE,NEIGHBOR,COMM,MODIFY,OUTPUT,TIME_N};

  Timer(class LAMMPS *);
  ~Timer();
  void init();
  void stamp();
  void stamp(enum ttype);
  void barrier_start(enum ttype);
  void barrier_stop(enum ttype);
  double elapsed(enum ttype);
  double cpu(enum ttype);
  double sys(enum ttype);

  double get_cpu(enum ttype which) const {
    return cpu_array[which]; };
  double get_sys(enum ttype which) const {
    return sys_array[which]; };
  double get_wall(enum ttype which) const {
    return wall_array[which]; };

  void set_wall(enum ttype, double);

 private:
  double *cpu_array;
  double *sys_array;
  double *wall_array;
  double previous_cpu;
  double previous_sys;
  double previous_wall;
};

}

#endif
