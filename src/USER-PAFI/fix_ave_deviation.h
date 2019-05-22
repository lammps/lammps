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

FixStyle(ave/deviation,FixAveDeviation)

#else

#ifndef LMP_FIX_AVE_DEVIATION_H
#define LMP_FIX_AVE_DEVIATION_H

#include <cstdio>
#include "fix.h"

namespace LAMMPS_NS {

class FixAveDeviation : public Fix {
 public:
  FixAveDeviation(class LAMMPS *, int, char **);
  ~FixAveDeviation();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

  double compute_array(int,int);

 private:
  int nvalues;
  int nrepeat,irepeat;
  bigint nvalid,nvalid_last;
  double **array;
  bigint nextvalid();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid timestep reset for fix ave/deviation

Resetting the timestep has invalidated the sequence of timesteps this
fix needs to process.

*/
