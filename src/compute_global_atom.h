/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(global/atom,ComputeGlobalAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_GLOBAL_ATOM_H
#define LMP_COMPUTE_GLOBAL_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeGlobalAtom : public Compute {
 public:
  ComputeGlobalAtom(class LAMMPS *, int, char **);
  virtual ~ComputeGlobalAtom();
  void init();
  void compute_peratom();
  double memory_usage();

 protected:
  int whichref, indexref, ref2index;
  char *idref;

  int nvalues;
  int *which, *argindex, *value2index;
  char **ids;

  int nmax, maxvector;
  int *indices;
  double *varatom;
  double *vecglobal;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute ID for compute global/atom does not exist

UNDOCUMENTED

E: Compute global/atom compute does not calculate a per-atom vector or array

UNDOCUMENTED

E: Compute global/atom compute does not calculate a per-atom vector

UNDOCUMENTED

E: Compute global/atom compute does not calculate a per-atom array

UNDOCUMENTED

E: Compute global/atom compute array is accessed out-of-range

UNDOCUMENTED

E: Fix ID for compute global/atom does not exist

UNDOCUMENTED

E: Compute global/atom fix does not calculate a per-atom vector or array

UNDOCUMENTED

E: Compute global/atom fix does not calculate a per-atom vector

UNDOCUMENTED

E: Compute global/atom fix does not calculate a per-atom array

UNDOCUMENTED

E: Compute global/atom fix array is accessed out-of-range

UNDOCUMENTED

E: Variable name for compute global/atom does not exist

UNDOCUMENTED

E: Compute global/atom variable is not atom-style variable

UNDOCUMENTED

E: Compute global/atom compute does not calculate a global vector

UNDOCUMENTED

E: Compute global/atom compute does not calculate a global array

UNDOCUMENTED

E: Compute global/atom fix does not calculate a global vector

UNDOCUMENTED

E: Compute global/atom fix does not calculate a global array

UNDOCUMENTED

E: Compute global/atom variable is not vector-style variable

UNDOCUMENTED

E: Fix used in compute global/atom not computed at compatible time

UNDOCUMENTED

U: Region ID for compute reduce/region does not exist

Self-explanatory.

U: Compute reduce replace requires min or max mode

Self-explanatory.

U: Invalid replace values in compute reduce

Self-explanatory.

U: Compute ID for compute reduce does not exist

Self-explanatory.

U: Compute reduce compute does not calculate a per-atom vector

Self-explanatory.

U: Compute reduce compute does not calculate a per-atom array

Self-explanatory.

U: Compute reduce compute array is accessed out-of-range

An index for the array is out of bounds.

U: Compute reduce compute does not calculate a local vector

Self-explanatory.

U: Compute reduce compute does not calculate a local array

Self-explanatory.

U: Compute reduce compute calculates global values

A compute that calculates peratom or local values is required.

U: Fix ID for compute reduce does not exist

Self-explanatory.

U: Compute reduce fix does not calculate a per-atom vector

Self-explanatory.

U: Compute reduce fix does not calculate a per-atom array

Self-explanatory.

U: Compute reduce fix array is accessed out-of-range

An index for the array is out of bounds.

U: Compute reduce fix does not calculate a local vector

Self-explanatory.

U: Compute reduce fix does not calculate a local array

Self-explanatory.

U: Compute reduce fix calculates global values

A fix that calculates peratom or local values is required.

U: Variable name for compute reduce does not exist

Self-explanatory.

U: Compute reduce variable is not atom-style variable

Self-explanatory.

U: Fix used in compute reduce not computed at compatible time

Fixes generate their values on specific timesteps.  Compute reduce is
requesting a value on a non-allowed timestep.

*/
