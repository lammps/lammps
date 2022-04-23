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
  ~ComputePropertyChunk() override;
  void init() override;
  void compute_vector() override;
  void compute_array() override;

  void lock_enable() override;
  void lock_disable() override;
  int lock_length() override;
  void lock(class Fix *, bigint, bigint) override;
  void unlock(class Fix *) override;

  double memory_usage() override;

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
