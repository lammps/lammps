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

#ifndef LMP_CHUNK_WRAP_H
#define LMP_CHUNK_WRAP_H

#include "pointers.h"

namespace LAMMPS_NS {

class ChunkWrap : protected Pointers {
 public:
  enum {UNWRAP, ATOM, CHUNK, COM};

  ChunkWrap(LAMMPS *, const char *, const char *);
  virtual ~ChunkWrap();

  void init();
  void wrap(int, double *);

 private:
  int wrapflag, nchunk, maxchunk;
  int *ichunk;
  char *idchunk;
  imageint *imgdiff;
  double **x;
  imageint *image;
};
}

#endif
