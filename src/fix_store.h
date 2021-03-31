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

FixStyle(STORE,FixStore)

#else

#ifndef LMP_FIX_STORE_H
#define LMP_FIX_STORE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixStore : public Fix {
 public:
  int nrow,ncol;         // size of global data array
  int nvalues;           // number of per-atom values
  double *vstore;        // vector storage for GLOBAL or PERATOM
  double **astore;       // array storage for GLOBAL or PERATOM
  int disable;        // 1 if operations (except grow) are currently disabled

  FixStore(class LAMMPS *, int, char **);
  ~FixStore();
  int setmask();
  void reset_global(int, int);

  void write_restart(FILE *);
  void restart(char *);

  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();

  double memory_usage();

 private:
  int flavor;                   // GLOBAL or PERATOM
  int vecflag;                  // 1 if ncol=1 or nvalues=1

  double *rbuf;                 // restart buffer for GLOBAL vec/array
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
