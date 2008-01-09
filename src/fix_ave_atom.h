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

#ifndef FIX_AVE_ATOM_H
#define FIX_AVE_ATOM_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixAveAtom : public Fix {
 public:
  FixAveAtom(class LAMMPS *, int, char **);
  ~FixAveAtom();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

 private:
  int nvalues;
  int nrepeat,nvalid,irepeat;
  int *which,*argindex,*value2index;
  char **ids;

  double **vector;
};

}

#endif
