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

FixStyle(deposit,FixDeposit)

#else

#ifndef LMP_FIX_DEPOSIT_H
#define LMP_FIX_DEPOSIT_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixDeposit : public Fix {
 public:
  FixDeposit(class LAMMPS *, int, char **);
  ~FixDeposit();
  int setmask();
  void init();
  void pre_exchange();
  void write_restart(FILE *);
  void restart(char *);

 private:
  int ninsert,ntype,nfreq,seed;
  int iregion,globalflag,localflag,maxattempt,rateflag,scaleflag,targetflag;
  char *idregion;
  double lo,hi,deltasq,nearsq,rate;
  double vxlo,vxhi,vylo,vyhi,vzlo,vzhi;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double tx,ty,tz;
  int nfirst,ninserted;
  int idnext,maxtag_all;
  class RanPark *random;

  void options(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Must specify a region in fix deposit

The region keyword must be specified with this fix.

E: Fix deposit region does not support a bounding box

Not all regions represent bounded volumes.  You cannot use
such a region with the fix deposit command.

E: Fix deposit region cannot be dynamic

Only static regions can be used with fix deposit.

E: Deposition region extends outside simulation box

Self-explanatory.

E: Use of fix deposit with undefined lattice

Must use lattice command with compute fix deposit command if units
option is set to lattice.

E: Region ID for fix deposit does not exist

Self-explanatory.

W: Particle deposition was unsuccessful

The fix deposit command was not able to insert as many atoms as
needed.  The requested volume fraction may be too high, or other atoms
may be in the insertion region.

*/
