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
  void post_integrate_respa(int,int);
  void *extract(const char *,int &);

 private:
  int gbit,gbitinverse;
  int regionflag,varflag,propflag,typeflag;
  int iregion,ivar,iprop;
  char *idregion,*idvar,*idprop;
  class Region *region;

  int nlevels_respa;

  void set_group();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for group dynamic does not exist

Self-explanatory.

E: Variable name for group dynamic does not exist

Self-explanatory.

E: Per atom property for group dynamic does not exist

Self-explanatory.

E: Group dynamic parent group cannot be dynamic

Self-explanatory.

E: Variable for group dynamic is invalid style

The variable must be an atom-style variable.

W: One or more dynamic groups may not be updated at correct point in timestep

If there are other fixes that act immediately after the initial stage
of time integration within a timestep (i.e. after atoms move), then
the command that sets up the dynamic group should appear after those
fixes.  This will insure that dynamic group assignments are made
after all atoms have moved.

*/
