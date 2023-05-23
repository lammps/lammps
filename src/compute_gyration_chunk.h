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
ComputeStyle(gyration/chunk,ComputeGyrationChunk);
// clang-format on
#else

#ifndef LMP_COMPUTE_GYRATION_CHUNK_H
#define LMP_COMPUTE_GYRATION_CHUNK_H

#include "compute_chunk.h"

namespace LAMMPS_NS {

class ComputeGyrationChunk : public ComputeChunk {
 public:
  ComputeGyrationChunk(class LAMMPS *, int, char **);
  ~ComputeGyrationChunk() override;

  void compute_vector() override;
  void compute_array() override;

  double memory_usage() override;

 private:
  int tensor;

  double *massproc, *masstotal;
  double **com, **comall;
  double *rg, *rgall;
  double **rgt, **rgtall;

  void com_chunk();
  void allocate() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
