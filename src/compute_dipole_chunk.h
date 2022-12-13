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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(dipole/chunk,ComputeDipoleChunk);
// clang-format on
#else

#ifndef LMP_COMPUTE_DIPOLE_CHUNK_H
#define LMP_COMPUTE_DIPOLE_CHUNK_H

#include "compute.h"

namespace LAMMPS_NS {
class Fix;

class ComputeDipoleChunk : public Compute {
 public:
  ComputeDipoleChunk(class LAMMPS *, int, char **);
  ~ComputeDipoleChunk() override;
  void init() override;
  void compute_array() override;

  void lock_enable() override;
  void lock_disable() override;
  int lock_length() override;
  void lock(Fix *, bigint, bigint) override;
  void unlock(Fix *) override;

  double memory_usage() override;

 private:
  int nchunk, maxchunk;
  char *idchunk;
  class ComputeChunkAtom *cchunk;

  double *massproc, *masstotal;
  double *chrgproc, *chrgtotal;
  double **com, **comall;
  double **dipole, **dipoleall;
  int usecenter;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
