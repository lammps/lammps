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
FixStyle(spring/chunk,FixSpringChunk);
// clang-format on
#else

#ifndef LMP_FIX_SPRING_CHUNK_H
#define LMP_FIX_SPRING_CHUNK_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSpringChunk : public Fix {
 public:
  FixSpringChunk(class LAMMPS *, int, char **);
  ~FixSpringChunk() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  void write_restart(FILE *) override;
  void restart(char *) override;
  double compute_scalar() override;

 private:
  int ilevel_respa;
  double k_spring;
  double esprings;
  char *idchunk, *idcom;

  int nchunk;
  double **com0, **fcom;

  class ComputeChunkAtom *cchunk;
  class ComputeCOMChunk *ccom;
};

}    // namespace LAMMPS_NS

#endif
#endif
