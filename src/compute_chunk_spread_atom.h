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
ComputeStyle(chunk/spread/atom,ComputeChunkSpreadAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_CHUNK_SPREAD_ATOM_H
#define LMP_COMPUTE_CHUNK_SPREAD_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeChunkSpreadAtom : public Compute {
 public:
  ComputeChunkSpreadAtom(class LAMMPS *, int, char **);
  ~ComputeChunkSpreadAtom() override;
  void init() override;
  void compute_peratom() override;
  double memory_usage() override;

 protected:
  struct value_t {
    int which;
    int argindex;
    std::string id;
    union {
      class Compute *c;
      class Fix *f;
    } val;
  };
  std::vector<value_t> values;

  char *idchunk;
  class ComputeChunkAtom *cchunk;
  int nmax;

  void init_chunk();
};

}    // namespace LAMMPS_NS

#endif
#endif
