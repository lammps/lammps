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

#ifndef FIX_HEAT_H
#define FIX_HEAT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixHeat : public Fix {
 public:
  FixHeat(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void end_of_step();
  double compute_scalar();

 private:
  double heat_input;
  double masstotal;
  double scale;
};

}

#endif
