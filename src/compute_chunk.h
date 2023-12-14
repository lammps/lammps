/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_COMPUTE_CHUNK_H
#define LMP_COMPUTE_CHUNK_H

#include "compute.h"

namespace LAMMPS_NS {
class ComputeChunkAtom;
class Fix;

class ComputeChunk : public Compute {
 public:
  char *idchunk;    // fields accessed by other classes

  ComputeChunk(class LAMMPS *, int, char **);
  ~ComputeChunk() override;
  void init() override;
  void compute_vector() override;
  void compute_array() override;

  void lock_enable() override;
  void lock_disable() override;
  int lock_length() override;
  void lock(Fix *, bigint, bigint) override;
  void unlock(Fix *) override;

  double memory_usage() override;

 protected:
  int nchunk, maxchunk;
  int firstflag, massneed;
  ComputeChunkAtom *cchunk;

  virtual void allocate(){};
};
}    // namespace LAMMPS_NS
#endif
