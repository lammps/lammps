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

#ifndef TEMP_RAMP_H
#define TEMP_RAMP_H

#include "temperature.h"

class TempRamp : public Temperature {
 public:
  TempRamp(int, char **);
  ~TempRamp() {}
  void init();
  double compute();
  void tensor();

 private:
  int coord_dim;
  double coord_lo,coord_hi;
  int v_dim;
  double v_lo,v_hi;

  void recount();
};

#endif
