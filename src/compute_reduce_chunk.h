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

ComputeStyle(reduce/chunk,ComputeReduceChunk)

#else

#ifndef LMP_COMPUTE_REDUCE_CHUNK_H
#define LMP_COMPUTE_REDUCE_CHUNK_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeReduceChunk : public Compute {
 public:
  ComputeReduceChunk(class LAMMPS *, int, char **);
  ~ComputeReduceChunk();
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
  int mode,nvalues;
  int *which,*argindex,*value2index;
  char *idchunk;
  char **ids;

  int nchunk;
  int maxchunk,maxatom;
  double initvalue;
  double *vlocal,*vglobal;
  double **alocal,**aglobal;
  double *varatom;

  class ComputeChunkAtom *cchunk;
  int *ichunk;

  void init_chunk();
  void compute_one(int, double *, int);
  void combine(double &, double);
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
