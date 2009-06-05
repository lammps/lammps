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

#ifndef FIX_DT_RESET_H
#define FIX_DT_RESET_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDtReset : public Fix {
 public:
  FixDtReset(class LAMMPS *, int, char **);
  ~FixDtReset() {}
  int setmask();
  void init();
  void setup(int);
  void end_of_step();
  double compute_scalar();
  double compute_vector(int);

 private:
  int minbound,maxbound,laststep;
  double tmin,tmax,xmax;
  double ftm2v;
  double t_elapsed;
  int respaflag;
};

}

#endif
