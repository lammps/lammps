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

#ifndef FIX_RESPA_H
#define FIX_RESPA_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRespa : public Fix {
  friend class Respa;
  friend class FixShake;

 public:
  FixRespa(class LAMMPS *, int, char **);
  ~FixRespa();
  int setmask();
  void init() {}

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

 private:
  int nlevels;
  double ***f_level;            // force at each rRESPA level
};

}

#endif
