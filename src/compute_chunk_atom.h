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

ComputeStyle(chunk/atom,ComputeChunkAtom)

#else

#ifndef LMP_COMPUTE_CHUNK_ATOM_H
#define LMP_COMPUTE_CHUNK_ATOM_H

#include "compute.h"
#include <map>

namespace LAMMPS_NS {

class ComputeChunkAtom : public Compute {
 public:
  int nchunk,ncoord,compress,idsflag,lockcount;
  int computeflag;    // 1 if this compute invokes other computes
  double chunk_volume_scalar;
  double *chunk_volume_vec;
  double **coord;
  int *ichunk,*chunkID;

  ComputeChunkAtom(class LAMMPS *, int, char **);
  ~ComputeChunkAtom();
  void init();
  void setup();
  void compute_peratom();
  void set_arrays(int);
  double memory_usage();

  void lock(class Fix *, bigint, bigint);
  void unlock(class Fix *);
  int setup_chunks();
  void compute_ichunk();

 private:
  int which,binflag;
  int regionflag,nchunksetflag,nchunkflag,discard;
  int limit,limitstyle,limitfirst;
  int scaleflag;
  double xscale,yscale,zscale;
  int argindex;
  char *cfvid;

  int ndim;
  int dim[3],originflag[3],nlayers[3];
  int minflag[3],maxflag[3];
  double origin[3],delta[3];
  double offset[3],invdelta[3];
  double minvalue[3],maxvalue[3];

  char *idregion;
  class Region *region;

  class Compute *cchunk;
  class Fix *fchunk;
  int vchunk;
  int maxvar;
  double *varatom;

  char *id_fix;
  class FixStore *fixstore;

  class Fix *lockfix;         // ptr to FixAveChunk that is locking out setups
                              // NULL if no lock currently in place
  bigint lockstart,lockstop;  // timesteps for start and stop of locking

  bigint invoked_setup;    // last timestep setup_chunks and nchunk calculated
  bigint invoked_ichunk;   // last timestep ichunk values calculated
  int nmax,nmaxint;
  double *chunk;

  int molcheck;              // one-time check if all molecule atoms in chunk
  int *exclude;              // 1 if atom is not assigned to any chunk
  std::map<int,int> *hash;   // store original chunks IDs before compression

  // static variable for ring communication callback to access class data
  // callback functions for ring communication

  static ComputeChunkAtom *cptr;
  static void idring(int, char *);

  void assign_chunk_ids();
  void compress_chunk_ids();
  void check_molecules();
  int setup_bins();
  void bin_volumes();
  void atom2bin1d();
  void atom2bin2d();
  void atom2bin3d();
  void readdim(int, char **, int, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

W: More than one compute ke/atom

It is not efficient to use compute ke/atom more than once.

*/
