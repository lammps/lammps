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

#ifndef FIX_DEPOSIT_H
#define FIX_DEPOSIT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDeposit : public Fix {
 public:
  FixDeposit(class LAMMPS *, int, char **);
  ~FixDeposit();
  int setmask();
  void pre_exchange();
  void write_restart(FILE *);
  void restart(char *);

 private:
  int ninsert,ntype,nfreq,seed;
  int iregion,globalflag,localflag,maxattempt,rateflag,scaleflag;
  double lo,hi,deltasq,nearsq,rate;
  double vxlo,vxhi,vylo,vyhi,vzlo,vzhi;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  int nfirst,ninserted;
  class RanPark *random;

  void options(int, char **);
};

}

#endif
