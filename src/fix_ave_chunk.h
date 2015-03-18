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

FixStyle(ave/chunk,FixAveChunk)

#else

#ifndef LMP_FIX_AVE_CHUNK_H
#define LMP_FIX_AVE_CHUNK_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixAveChunk : public Fix {
 public:
  FixAveChunk(class LAMMPS *, int, char **);
  ~FixAveChunk();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();
  double compute_array(int,int);
  double memory_usage();
  void reset_timestep(bigint);

 private:
  int me,nvalues;
  int nrepeat,nfreq,irepeat;
  int normflag,scaleflag,overwrite,biasflag,colextra;
  bigint nvalid,nvalid_last;
  double adof,cdof;
  char *tstring,*sstring,*id_bias;
  int *which,*argindex,*value2index;
  char **ids;
  class Compute *tbias;     // ptr to additional bias compute
  FILE *fp;

  int ave,nwindow;
  int normcount,iwindow,window_limit;

  int nchunk,maxchunk;
  char *idchunk;
  class ComputeChunkAtom *cchunk;

  long filepos;

  int maxvar;
  double *varatom;

  // one,many,sum vecs/arrays are used with a single Nfreq epoch
  // total,list vecs/arrays are used across epochs

  double *count_one,*count_many,*count_sum;
  double **values_one,**values_many,**values_sum;
  double *count_total,**count_list;
  double **values_total,***values_list;

  void allocate();
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

E: Cannot open fix ave/chunk file %s

UNDOCUMENTED

E: Could not find compute ID for temperature bias

UNDOCUMENTED

E: Bias compute does not calculate temperature

UNDOCUMENTED

E: Bias compute does not calculate a velocity bias

UNDOCUMENTED

E: Compute ID for fix ave/chunk does not exist

UNDOCUMENTED

E: Fix ave/chunk compute does not calculate per-atom values

UNDOCUMENTED

E: Fix ave/chunk compute does not calculate a per-atom vector

UNDOCUMENTED

E: Fix ave/chunk compute does not calculate a per-atom array

UNDOCUMENTED

E: Fix ave/chunk compute vector is accessed out-of-range

UNDOCUMENTED

E: Fix ID for fix ave/chunk does not exist

UNDOCUMENTED

E: Fix ave/chunk fix does not calculate per-atom values

UNDOCUMENTED

E: Fix ave/chunk fix does not calculate a per-atom vector

UNDOCUMENTED

E: Fix ave/chunk fix does not calculate a per-atom array

UNDOCUMENTED

E: Fix ave/chunk fix vector is accessed out-of-range

UNDOCUMENTED

E: Variable name for fix ave/chunk does not exist

UNDOCUMENTED

E: Fix ave/chunk variable is not atom-style variable

UNDOCUMENTED

E: Chunk/atom compute does not exist for fix ave/chunk

UNDOCUMENTED

E: Fix ave/chunk does not use chunk/atom compute

UNDOCUMENTED

E: Fix for fix ave/chunk not computed at compatible time

UNDOCUMENTED

E: Fix ave/chunk missed timestep

UNDOCUMENTED

U: Cannot use fix ave/spatial z for 2 dimensional model

Self-explanatory.

U: Same dimension twice in fix ave/spatial

Self-explanatory.

U: Region ID for fix ave/spatial does not exist

Self-explanatory.

U: Cannot open fix ave/spatial file %s

The specified file cannot be opened.  Check that the path and name are
correct.

U: Compute ID for fix ave/spatial does not exist

Self-explanatory.

U: Fix ave/spatial compute does not calculate per-atom values

A compute used by fix ave/spatial must generate per-atom values.

U: Fix ave/spatial compute does not calculate a per-atom vector

A compute used by fix ave/spatial must generate per-atom values.

U: Fix ave/spatial compute does not calculate a per-atom array

Self-explanatory.

U: Fix ave/spatial compute vector is accessed out-of-range

The index for the vector is out of bounds.

U: Fix ID for fix ave/spatial does not exist

Self-explanatory.

U: Fix ave/spatial fix does not calculate per-atom values

A fix used by fix ave/spatial must generate per-atom values.

U: Fix ave/spatial fix does not calculate a per-atom vector

A fix used by fix ave/spatial must generate per-atom values.

U: Fix ave/spatial fix does not calculate a per-atom array

Self-explanatory.

U: Fix ave/spatial fix vector is accessed out-of-range

The index for the vector is out of bounds.

U: Variable name for fix ave/spatial does not exist

Self-explanatory.

U: Fix ave/spatial variable is not atom-style variable

A variable used by fix ave/spatial must generate per-atom values.

U: Fix ave/spatial for triclinic boxes requires units reduced

Self-explanatory.

U: Fix ave/spatial settings invalid with changing box size

If the box size changes, only the units reduced option can be
used.

U: Fix for fix ave/spatial not computed at compatible time

Fixes generate their values on specific timesteps.  Fix ave/spatial is
requesting a value on a non-allowed timestep.

U: Fix ave/spatial missed timestep

You cannot reset the timestep to a value beyond where the fix
expects to next perform averaging.

*/
