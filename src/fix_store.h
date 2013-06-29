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

FixStyle(STORE,FixStore)

#else

#ifndef LMP_FIX_STORE_H
#define LMP_FIX_STORE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixStore : public Fix {
 public:
  double *vstore;        // vector storage if nvalues = 1
  double **astore;       // array storage if nvalues > 1

  FixStore(class LAMMPS *, int, char **);
  ~FixStore();
  int setmask();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();

 private:
  int nvalues;                  // total # of values per atom
  int vecflag;                  // 1 if nvalues = 1
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

UNDOCUMENTED

*/
