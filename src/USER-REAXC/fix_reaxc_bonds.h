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

FixStyle(reax/c/bonds,FixReaxCBonds)

#else

#ifndef LMP_FIX_REAXC_BONDS_H
#define LMP_FIX_REAXC_BONDS_H

#include "stdio.h"
#include "fix.h"
#include "pair_reax_c.h"
#include "reaxc_defs.h"

namespace LAMMPS_NS {

class FixReaxCBonds : public Fix {
 public:
  FixReaxCBonds(class LAMMPS *, int, char **);
  ~FixReaxCBonds();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();

 private:
  int me;
  int nfreq;
  FILE *fp;

  void OutputReaxCBonds(bigint, FILE*);
  int nint(const double&);

  class PairReaxC *reaxc;
};

}

#endif
#endif
