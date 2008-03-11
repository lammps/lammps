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

#ifndef FIX_TEMP_RESCALE_H
#define FIX_TEMP_RESCALE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTempRescale : public Fix {
 public:
  FixTempRescale(class LAMMPS *, int, char **);
  ~FixTempRescale();
  int setmask();
  void init();
  void end_of_step();
  int modify_param(int, char **);
  void reset_target(double);
  double compute_scalar();

 private:
  int which;
  double t_start,t_stop,t_window;
  double fraction,energy,efactor;

  char *id_temp;
  class Compute *temperature;
  int tflag;
};

}

#endif
