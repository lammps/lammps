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

#ifdef FIX_CLASS
// clang-format off
FixStyle(momentum/chunk,FixMomentumChunk);
// clang-format on
#else

#ifndef LMP_FIX_MOMENTUM_CHUNK_H
#define LMP_FIX_MOMENTUM_CHUNK_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMomentumChunk : public Fix {
 public:
  FixMomentumChunk(class LAMMPS *, int, char **);
  ~FixMomentumChunk() override;
  int setmask() override;
  void init() override;
  void end_of_step() override;
  void post_run() override;

 protected:
  std::string id_chunk, id_com, id_vcm, id_omega;
  int nchunk, linear, angular, rescale;
  int xflag, yflag, zflag;

  class ComputeChunkAtom *cchunk;
  class Compute *ccom, *cvcm, *comega;
};

}    // namespace LAMMPS_NS

#endif
#endif
