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

#ifdef COMMAND_CLASS

CommandStyle(set,Set)

#else

#ifndef LMP_SET_H
#define LMP_SET_H

#include "pointers.h"

namespace LAMMPS_NS {

class Set : protected Pointers {
 public:
  Set(class LAMMPS *lmp) : Pointers(lmp) {};
  void command(int, char **);

 private:
  char *id;
  int *select;
  int style,ivalue,newtype,count;
  int ximage,yimage,zimage,ximageflag,yimageflag,zimageflag;
  double dvalue,xvalue,yvalue,zvalue,wvalue,fraction;
  double PI;

  void selection(int);
  void set(int);
  void setrandom(int);
  void topology(int);
};

}

#endif
#endif
