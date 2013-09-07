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

#ifdef FIX_CLASS

FixStyle(countdown,FixCountdown)

#else

#ifndef LMP_FIX_COUNTDOWN_H
#define LMP_FIX_COUNTDOWN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixCountdown : public Fix {
 public:
  FixCountdown(class LAMMPS *, int, char **);
  ~FixCountdown();
  int setmask();
  void setup(int);
  void end_of_step();
  double compute_scalar();
  double compute_vector(int);

 private:
  double *_times;                // list of times per time step
  double _last_wall, _complete;  // last wall time; completion fraction
  double _eta[4];                // estimated time: days, hours, mins, secs
  bigint _outfreq;               // output frequency
  int    _ntimes, _idx;          // number of recorded times, index in array
};

}

#endif
#endif
