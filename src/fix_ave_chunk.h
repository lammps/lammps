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

 private:
  int me,nvalues;
  int nrepeat,nfreq,irepeat;
  int normflag,scaleflag,overwrite,biasflag,colextra;
  bigint nvalid,nvalid_last;
  double adof,cdof;
  char *format,*format_user;
  char *tstring,*sstring,*id_bias;
  int *which,*argindex,*value2index;
  char **ids;
  class Compute *tbias;     // ptr to additional bias compute
  FILE *fp;

  int densityflag;        // 1 if density/number or density/mass requested
  int volflag;            // SCALAR/VECTOR for density normalization by volume
  double chunk_volume_scalar;
  double *chunk_volume_vec;

  int ave,nwindow;
  int normcount,iwindow,window_limit;

  int nchunk,maxchunk;
  char *idchunk;
  class ComputeChunkAtom *cchunk;
  int lockforever;

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

E: No values in fix ave/chunk command

Self-explanatory.

E: Cannot open fix ave/chunk file %s

The specified file cannot be opened.  Check that the path and name are
correct.

E: Could not find compute ID for temperature bias

Self-explanatory.

E: Bias compute does not calculate temperature

The specified compute must compute temperature.

E: Bias compute does not calculate a velocity bias

The specified compute must compute a bias for temperature.

E: Compute ID for fix ave/chunk does not exist

Self-explanatory.

E: Fix ave/chunk compute does not calculate per-atom values

Self-explanatory.

E: Fix ave/chunk compute does not calculate a per-atom vector

Self-explanatory.

E: Fix ave/chunk compute does not calculate a per-atom array

Self-explanatory.

E: Fix ave/chunk compute vector is accessed out-of-range

Self-explanatory.

E: Fix ID for fix ave/chunk does not exist

Self-explanatory.

E: Fix ave/chunk fix does not calculate per-atom values

Self-explanatory.

E: Fix ave/chunk fix does not calculate a per-atom vector

Self-explanatory.

E: Fix ave/chunk fix does not calculate a per-atom array

Self-explanatory.

E: Fix ave/chunk fix vector is accessed out-of-range

Self-explanatory.

E: Variable name for fix ave/chunk does not exist

Self-explanatory.

E: Fix ave/chunk variable is not atom-style variable

Self-explanatory.

E: Chunk/atom compute does not exist for fix ave/chunk

Self-explanatory.

E: Fix ave/chunk does not use chunk/atom compute

The specified compute is not for a compute chunk/atom command.

E: Error writing file header

Something in the output to the file triggered an error.

E: Fix for fix ave/chunk not computed at compatible time

Fixes generate their values on specific timesteps.  Fix ave/chunk is
requesting a value on a non-allowed timestep.

E: Invalid timestep reset for fix ave/chunk

Resetting the timestep has invalidated the sequence of timesteps this
fix needs to process.

E: Error writing averaged chunk data

Something in the output to the file triggered an error.

*/
