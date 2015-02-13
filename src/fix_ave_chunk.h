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
  bigint nvalid;
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

E: Cannot use fix ave/spatial z for 2 dimensional model

Self-explanatory.

E: Same dimension twice in fix ave/spatial

Self-explanatory.

E: Region ID for fix ave/spatial does not exist

Self-explanatory.

E: Cannot open fix ave/spatial file %s

The specified file cannot be opened.  Check that the path and name are
correct.

E: Compute ID for fix ave/spatial does not exist

Self-explanatory.

E: Fix ave/spatial compute does not calculate per-atom values

A compute used by fix ave/spatial must generate per-atom values.

E: Fix ave/spatial compute does not calculate a per-atom vector

A compute used by fix ave/spatial must generate per-atom values.

E: Fix ave/spatial compute does not calculate a per-atom array

Self-explanatory.

E: Fix ave/spatial compute vector is accessed out-of-range

The index for the vector is out of bounds.

E: Fix ID for fix ave/spatial does not exist

Self-explanatory.

E: Fix ave/spatial fix does not calculate per-atom values

A fix used by fix ave/spatial must generate per-atom values.

E: Fix ave/spatial fix does not calculate a per-atom vector

A fix used by fix ave/spatial must generate per-atom values.

E: Fix ave/spatial fix does not calculate a per-atom array

Self-explanatory.

E: Fix ave/spatial fix vector is accessed out-of-range

The index for the vector is out of bounds.

E: Variable name for fix ave/spatial does not exist

Self-explanatory.

E: Fix ave/spatial variable is not atom-style variable

A variable used by fix ave/spatial must generate per-atom values.

E: Fix ave/spatial for triclinic boxes requires units reduced

Self-explanatory.

E: Fix ave/spatial settings invalid with changing box size

If the box size changes, only the units reduced option can be
used.

E: Fix for fix ave/spatial not computed at compatible time

Fixes generate their values on specific timesteps.  Fix ave/spatial is
requesting a value on a non-allowed timestep.

E: Fix ave/spatial missed timestep

You cannot reset the timestep to a value beyond where the fix
expects to next perform averaging.

*/
