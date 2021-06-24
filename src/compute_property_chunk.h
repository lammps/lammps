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
ComputeStyle(property/chunk,ComputePropertyChunk);
// clang-format on
#else

#ifndef LMP_COMPUTE_CHUNK_MOLECULE_H
#define LMP_COMPUTE_CHUNK_MOLECULE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePropertyChunk : public Compute {
 public:
  ComputePropertyChunk(class LAMMPS *, int, char **);
  ~ComputePropertyChunk();
  void init();
  void compute_vector();
  void compute_array();

  void lock_enable();
  void lock_disable();
  int lock_length();
  void lock(class Fix *, bigint, bigint);
  void unlock(class Fix *);

  double memory_usage();

 private:
  int nchunk, maxchunk;
  char *idchunk;
  class ComputeChunkAtom *cchunk;
  int *ichunk;

  int nvalues, countflag;
  double *buf;
  int *count_one, *count_all;

  void allocate();

  typedef void (ComputePropertyChunk::*FnPtrPack)(int);
  FnPtrPack *pack_choice;    // ptrs to pack functions

  void pack_count(int);
  void pack_id(int);
  void pack_coord1(int);
  void pack_coord2(int);
  void pack_coord3(int);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute chunk/atom stores no IDs for compute property/chunk

It will only store IDs if its compress option is enabled.

E: Compute chunk/atom stores no coord1 for compute property/chunk

Only certain binning options for compute chunk/atom store coordinates.

E: Compute chunk/atom stores no coord2 for compute property/chunk

Only certain binning options for compute chunk/atom store coordinates.

E: Compute chunk/atom stores no coord3 for compute property/chunk

Only certain binning options for compute chunk/atom store coordinates.

E: Invalid keyword in compute property/chunk command

Self-explanatory.

E: Chunk/atom compute does not exist for compute property/chunk

Self-explanatory.

E: Compute property/chunk does not use chunk/atom compute

The style of the specified compute is not chunk/atom.

*/
