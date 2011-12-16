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

#ifdef COMPUTE_CLASS

ComputeStyle(slice,ComputeSlice)

#else

#ifndef LMP_COMPUTE_SLICE_H
#define LMP_COMPUTE_SLICE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSlice : public Compute {
 public:
  ComputeSlice(class LAMMPS *, int, char **);
  virtual ~ComputeSlice();
  void init();
  void compute_vector();
  void compute_array();

 private:
  int me;
  int nstart,nstop,nskip,nvalues;
  int *which,*argindex,*value2index;
  char **ids;

  void extract_one(int, double *, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute ID for compute slice does not exist

UNDOCUMENTED

E: Compute slice compute does not calculate a global array

UNDOCUMENTED

E: Compute slice compute vector is accessed out-of-range

UNDOCUMENTED

E: Compute slice compute does not calculate a global vector

UNDOCUMENTED

E: Compute slice compute array is accessed out-of-range

UNDOCUMENTED

E: Compute slice compute does not calculate global vector or array

UNDOCUMENTED

E: Fix ID for compute slice does not exist

UNDOCUMENTED

E: Compute slice fix does not calculate a global array

UNDOCUMENTED

E: Compute slice fix vector is accessed out-of-range

UNDOCUMENTED

E: Compute slice fix does not calculate a global vector

UNDOCUMENTED

E: Compute slice fix array is accessed out-of-range

UNDOCUMENTED

E: Compute slice fix does not calculate global vector or array

UNDOCUMENTED

E: Fix used in compute slice not computed at compatible time

UNDOCUMENTED

*/
