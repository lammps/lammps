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

#ifndef FIX_GYRATION_H
#define FIX_GYRATION_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixGyration : public Fix {
 public:
  FixGyration(class LAMMPS *, int, char **);
  ~FixGyration();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();

 private:
  int me,first;
  FILE *fp;
  double masstotal;
};

}

#endif
