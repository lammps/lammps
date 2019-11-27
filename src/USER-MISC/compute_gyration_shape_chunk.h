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

ComputeStyle(gyration/shape/chunk,ComputeGyrationShapeChunk)

#else

#ifndef LMP_COMPUTE_GYRATION_SHAPE_CHUNK_H
#define LMP_COMPUTE_GYRATION_SHAPE_CHUNK_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeGyrationShapeChunk : public Compute {
 public:
  char *id_gyration_chunk;              // fields accessed by other classes

  ComputeGyrationShapeChunk(class LAMMPS *, int, char **);
  ~ComputeGyrationShapeChunk();
  void init();
  void setup();
  void compute_array();

  int lock_length();

  double memory_usage();

 private:
  int current_nchunks, former_nchunks;
  int firstflag;
  double **shape_parameters;
  class Compute *c_gyration_chunk;

  void allocate();

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute gyration/chunk ID does not exist for compute gyration/shape/chunk

Self-explanatory.  Provide a valid compute ID

E: Compute gyration/shape/chunk ID does not point to a gyration/chunk compute

Self-explanatory.  Provide an ID of a compute gyration/chunk command.

E: Compute gyration/chunk does not compute gyration tensor

Self-explanatory. Use keyword tensor in compute gyration/chunk command
*/
