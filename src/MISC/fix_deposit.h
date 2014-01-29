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
  int ntype;      // type of deposited atom, visible to PairGran

  FixDeposit(class LAMMPS *, int, char **);
  ~FixDeposit();
  int setmask();
  void init();
  void pre_exchange();
  void write_restart(FILE *);
  void restart(char *);

 private:
  int ninsert,nfreq,seed;
  int iregion,globalflag,localflag,maxattempt,rateflag,scaleflag,targetflag;
  int mode,rigidflag,shakeflag,idnext;
  double lo,hi,deltasq,nearsq,rate;
  double vxlo,vxhi,vylo,vyhi,vzlo,vzhi;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double tx,ty,tz;
  char *idregion;
  char *idrigid,*idshake;

  class Molecule *onemol;
  int natom;
  double **coords;
  imageint *imageflags;
  class Fix *fixrigid,*fixshake;

  int nfirst,ninserted;
  tagint maxtag_all,maxmol_all;
  class RanPark *random;

  void find_maxid();
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

E: Invalid atom type in fix deposit command

UNDOCUMENTED

E: Must specify a region in fix deposit

The region keyword must be specified with this fix.

E: Fix deposit region does not support a bounding box

Not all regions represent bounded volumes.  You cannot use
such a region with the fix deposit command.

E: Fix deposit region cannot be dynamic

Only static regions can be used with fix deposit.

E: Deposition region extends outside simulation box

Self-explanatory.

E: Cannot use fix_deposit unless atoms have IDs

UNDOCUMENTED

E: Fix deposit molecule must have coordinates

UNDOCUMENTED

E: Fix deposit molecule must have atom types

UNDOCUMENTED

E: Invalid atom type in fix deposit mol command

UNDOCUMENTED

E: Fix deposit molecule template ID must be same as atom style template ID

UNDOCUMENTED

E: Cannot use fix deposit rigid and not molecule

UNDOCUMENTED

E: Cannot use fix deposit shake and not molecule

UNDOCUMENTED

E: Cannot use fix deposit rigid and shake

UNDOCUMENTED

E: Region ID for fix deposit does not exist

Self-explanatory.

E: Fix pour rigid fix does not exist

UNDOCUMENTED

E: Fix deposit and fix rigid/small not using same molecule template ID

UNDOCUMENTED

E: Fix deposit shake fix does not exist

UNDOCUMENTED

E: Fix deposit and fix shake not using same molecule template ID

UNDOCUMENTED

W: Particle deposition was unsuccessful

The fix deposit command was not able to insert as many atoms as
needed.  The requested volume fraction may be too high, or other atoms
may be in the insertion region.

E: Too many total atoms

See the setting for bigint in the src/lmptype.h file.

E: New atom IDs exceed maximum allowed ID

UNDOCUMENTED

E: Molecule template ID for fix deposit does not exist

UNDOCUMENTED

W: Molecule template for fix deposit has multiple molecules

UNDOCUMENTED

U: Use of fix deposit with undefined lattice

Must use lattice command with compute fix deposit command if units
option is set to lattice.

*/
