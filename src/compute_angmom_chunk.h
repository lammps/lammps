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
ComputeStyle(angmom/chunk,ComputeAngmomChunk);
// clang-format on
#else

#ifndef LMP_COMPUTE_ANGMOM_CHUNK_H
#define LMP_COMPUTE_ANGMOM_CHUNK_H

#include "compute_chunk.h"

namespace LAMMPS_NS {

class ComputeAngmomChunk : public ComputeChunk {
 public:
  ComputeAngmomChunk(class LAMMPS *, int, char **);
  ~ComputeAngmomChunk() override;
  void compute_array() override;

  double memory_usage() override;

 private:
  double *massproc, *masstotal;
  double **com, **comall;
  double **angmom, **angmomall;

  void allocate() override;
};
}    // namespace LAMMPS_NS
#endif
#endif
