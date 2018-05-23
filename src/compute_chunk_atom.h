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
  int scaleflag,pbcflag;
  double xscale,yscale,zscale;
  int argindex;
  char *cfvid;

  // xyz spatial bins

  int ndim;
  int dim[3],originflag[3],nlayers[3];
  int minflag[3],maxflag[3];
  double origin[3],delta[3];
  double offset[3],invdelta[3];
  double minvalue[3],maxvalue[3];

  // spherical spatial bins

  double sorigin_user[3];
  double sorigin[3];
  double sradmin_user,sradmax_user;
  double sradmin,sradmax,sinvrad;
  int nsbin;

  // cylindrical spatial bins

  double corigin_user[3];
  double corigin[3];
  double cradmin_user,cradmax_user;
  double cradmin,cradmax,cinvrad;
  int cdim1,cdim2;
  int ncbin,ncplane;

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
  std::map<tagint,int> *hash;   // store original chunks IDs before compression

  // callback function for ring communication

  static void idring(int, char *, void *);

  void assign_chunk_ids();
  void compress_chunk_ids();
  void check_molecules();
  int setup_xyz_bins();
  int setup_sphere_bins();
  int setup_cylinder_bins();
  void bin_volumes();
  void atom2bin1d();
  void atom2bin2d();
  void atom2bin3d();
  void atom2binsphere();
  void atom2bincylinder();
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

E: Region ID for compute chunk/atom does not exist

Self-explanatory.

E: Compute chunk/atom molecule for non-molecular system

Self-explanatory.

E: Compute chunk/atom without bins cannot use discard mixed

That discard option only applies to the binning styles.

E: Compute chunk/atom sphere z origin must be 0.0 for 2d

Self-explanatory.

E: Compute chunk/atom cylinder axis must be z for 2d

Self-explanatory.

E: Compute ID for compute chunk /atom does not exist

Self-explanatory.

E: Compute chunk/atom compute does not calculate per-atom values

Self-explanatory.

E: Compute chunk/atom compute does not calculate a per-atom vector

Self-explanatory.

E: Compute chunk/atom compute does not calculate a per-atom array

Self-explanatory.

E: Compute chunk/atom compute array is accessed out-of-range

The index for the array is out of bounds.

E: Fix ID for compute chunk/atom does not exist

Self-explanatory.

E: Compute chunk/atom fix does not calculate per-atom values

Self-explanatory.

E: Compute chunk/atom fix does not calculate a per-atom vector

Self-explanatory.

E: Compute chunk/atom fix does not calculate a per-atom array

Self-explanatory.

E: Compute chunk/atom fix array is accessed out-of-range

the index for the array is out of bounds.

E: Variable name for compute chunk/atom does not exist

Self-explanatory.

E: Compute chunk/atom variable is not atom-style variable

Self-explanatory.

E: Compute chunk/atom for triclinic boxes requires units reduced

Self-explanatory.

E: Compute ID for compute chunk/atom does not exist

Self-explanatory.

E: Molecule IDs too large for compute chunk/atom

The IDs must not be larger than can be stored in a 32-bit integer
since chunk IDs are 32-bit integers.

E: Compute chunk/atom ids once but nchunk is not once

You cannot assign chunks IDs to atom permanently if the number of
chunks may change.

E: Two fix commands using same compute chunk/atom command in incompatible ways

UNDOCUMENTED

E: Fix used in compute chunk/atom not computed at compatible time

The chunk/atom compute cannot query the output of the fix on a timestep
it is needed.

W: One or more chunks do not contain all atoms in molecule

This may not be what you intended.

E: Invalid bin bounds in compute chunk/atom

The lo/hi values are inconsistent.

E: Compute chunk/atom bin/sphere radius is too large for periodic box

Radius cannot be bigger than 1/2 of any periodic dimension.

E: Compute chunk/atom bin/cylinder radius is too large for periodic box

Radius cannot be bigger than 1/2 of a non-axis  periodic dimension.

E: Cannot use compute chunk/atom bin z for 2d model

Self-explanatory.

U: Two fix ave commands using same compute chunk/atom command in incompatible ways

They are both attempting to "lock" the chunk/atom command so that the
chunk assignments persist for some number of timesteps, but are doing
it in different ways.

*/
