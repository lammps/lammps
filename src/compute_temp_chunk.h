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

ComputeStyle(temp/chunk,ComputeTempChunk)

#else

#ifndef LMP_COMPUTE_TEMP_CHUNK_H
#define LMP_COMPUTE_TEMP_CHUNK_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempChunk : public Compute {
 public:
  ComputeTempChunk(class LAMMPS *, int, char **);
  ~ComputeTempChunk();
  void init();
  double compute_scalar();
  void compute_vector();
  void compute_array();

  void remove_bias(int, double *);
  void remove_bias_all();
  void restore_bias(int, double *);
  void restore_bias_all();

  void lock_enable();
  void lock_disable();
  int lock_length();
  void lock(class Fix *, bigint, bigint);
  void unlock(class Fix *);

  double memory_usage();

 private:
  int nchunk,maxchunk,comflag,biasflag;
  int nvalues;
  int *which;
  char *idchunk;
  class ComputeChunkAtom *cchunk;
  double adof,cdof;
  char *id_bias;
  class Compute *tbias;     // ptr to additional bias compute
  bigint comstep;

  double *sum,*sumall;
  int *count,*countall;
  double *massproc,*masstotal;
  double **vcm,**vcmall;

  void vcm_compute();
  void temperature(int);
  void kecom(int);
  void internal(int);
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

E: Could not find compute ID for temperature bias

Self-explanatory.

E: Bias compute does not calculate temperature

The specified compute must compute temperature.

E: Bias compute does not calculate a velocity bias

The specified compute must compute a bias for temperature.

E: Cannot use both com and bias with compute temp/chunk

Self-explanatory.

E: Chunk/atom compute does not exist for compute temp/chunk

Self-explanatory.

E: Compute temp/chunk does not use chunk/atom compute

The style of the specified compute is not chunk/atom.

E: Temperature compute degrees of freedom < 0

This should not happen if you are calculating the temperature
on a valid set of atoms.

*/
