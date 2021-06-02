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
ComputeStyle(chunk/spread/atom,ComputeChunkSpreadAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_CHUNK_SPREAD_ATOM_H
#define LMP_COMPUTE_CHUNK_SPREAD_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeChunkSpreadAtom : public Compute {
 public:
  ComputeChunkSpreadAtom(class LAMMPS *, int, char **);
  ~ComputeChunkSpreadAtom();
  void init();
  void compute_peratom();
  double memory_usage();

 protected:
  int mode, nvalues;
  char *idchunk;
  char **ids;
  int *which, *argindex, *value2index;

  int nmax;
  class ComputeChunkAtom *cchunk;

  void init_chunk();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
