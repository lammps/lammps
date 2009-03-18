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

#ifndef FIX_EVAPORATE_H
#define FIX_EVAPORATE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEvaporate : public Fix {
 public:
  FixEvaporate(class LAMMPS *, int, char **);
  ~FixEvaporate();
  int setmask();
  void init();
  void pre_exchange();
  double compute_scalar();
  double memory_usage();

 private:
  int nevery,nflux,iregion;
  int ndeleted;

  int nmax;
  int *list,*mark;

  class RanPark *random;
};

}

#endif
