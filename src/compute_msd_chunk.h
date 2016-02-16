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

#ifdef COMPUTE_CLASS

ComputeStyle(msd/chunk,ComputeMSDChunk)

#else

#ifndef LMP_COMPUTE_MSD_CHUNK_H
#define LMP_COMPUTE_MSD_CHUNK_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeMSDChunk : public Compute {
 public:
  ComputeMSDChunk(class LAMMPS *, int, char **);
  ~ComputeMSDChunk();
  void init();
  void setup();
  void compute_array();

  void lock_enable();
  void lock_disable();
  int lock_length();
  void lock(class Fix *, bigint, bigint);
  void unlock(class Fix *);

  double memory_usage();

 private:
  int nchunk;
  char *idchunk;
  class ComputeChunkAtom *cchunk;
  char *id_fix;
  class FixStore *fix;
  int firstflag;

  double *massproc,*masstotal;
  double **com,**comall;
  double **msd;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Chunk/atom compute does not exist for compute msd/chunk

Self-explanatory.

E: Compute msd/chunk does not use chunk/atom compute

The style of the specified compute is not chunk/atom.

E: Could not find compute msd/chunk fix ID

The compute creates an internal fix, which has been deleted.

E: Compute msd/chunk nchunk is not static

This is required because the MSD cannot be computed consistently if
the number of chunks is changing.  Compute chunk/atom allows setting
nchunk to be static.

*/
