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

FixStyle(spring/chunk,FixSpringChunk)

#else

#ifndef LMP_FIX_SPRING_CHUNK_H
#define LMP_FIX_SPRING_CHUNK_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSpringChunk : public Fix {
 public:
  FixSpringChunk(class LAMMPS *, int, char **);
  ~FixSpringChunk();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();

 private:
  int ilevel_respa;
  double k_spring;
  double esprings;
  char *idchunk,*idcom;

  int nchunk;
  double **com0,**fcom;

  class ComputeChunkAtom *cchunk;
  class ComputeCOMChunk *ccom;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Chunk/atom compute does not exist for fix spring/chunk

UNDOCUMENTED

E: Fix spring/chunk does not use chunk/atom compute

UNDOCUMENTED

E: Com/chunk compute does not exist for fix spring/chunk

UNDOCUMENTED

E: Fix spring/chunk does not use com/chunk compute

UNDOCUMENTED

E: Fix spring chunk chunkID not same as comID chunkID

UNDOCUMENTED

U: R0 < 0 for fix spring command

Equilibrium spring length is invalid.

U: Fix spring couple group ID does not exist

Self-explanatory.

U: Two groups cannot be the same in fix spring couple

Self-explanatory.

*/
