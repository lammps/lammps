/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef FIX_TEMP_RESCALE_H
#define FIX_TEMP_RESCALE_H

#include "fix.h"

class Temperature;

class FixTempRescale : public Fix {
 public:
  FixTempRescale(int, char **);
  ~FixTempRescale() {}
  int setmask();
  void init();
  void end_of_step();
  int modify_param(int, char **);
  int thermo_fields(int, int *, char **);
  int thermo_compute(double *);

 private:
  int iregion;
  double t_start,t_end,t_window;
  double fraction,energy,efactor;
  Temperature *temperature;
};

#endif
