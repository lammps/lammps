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

#ifdef FIX_CLASS

FixStyle(GROUP,FixGroup)

#else

#ifndef LMP_FIX_GROUP_H
#define LMP_FIX_GROUP_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGroup : public Fix {
 public:
  FixGroup(class LAMMPS *, int, char **);
  ~FixGroup();
  int setmask();
  void init();
  void setup(int);
  void post_integrate();

 private:
  int gbit,gbitinverse;
  int regionflag,varflag;
  int iregion,ivar;
  char *idregion,*idvar;
  class Region *region;

  void set_group();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
