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

#ifndef FIX_DEPOSIT_H
#define FIX_DEPOSIT_H

#include "fix.h"

class RanPark;

class FixDeposit : public Fix {
 public:
  FixDeposit(int, char **);
  ~FixDeposit();
  int setmask();
  void pre_exchange();

 private:
  int ninsert,ntype,nfreq,seed;
  int iregion,globalflag,localflag,maxattempt,rateflag,scaleflag;
  double lo,hi,deltasq,nearsq,rate;
  double vxlo,vxhi,vylo,vyhi,vzlo,vzhi;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  int nfirst,ninserted;
  RanPark *random;

  void options(int, char **);
};

#endif
