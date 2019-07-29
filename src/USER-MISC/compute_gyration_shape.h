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

ComputeStyle(gyration/shape,ComputeGyrationShape)

#else

#ifndef LMP_COMPUTE_GYRATION_SHAPE_H
#define LMP_COMPUTE_GYRATION_SHAPE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeGyrationShape : public Compute {
 public:
  char *id_gyration;              // fields accessed by other classes

  ComputeGyrationShape(class LAMMPS *, int, char **);
  ~ComputeGyrationShape();
  void init();
  void compute_vector();

 private:
  class Compute *c_gyration;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute gyration ID does not exist for compute gyration/shape

Self-explanatory.  Provide a valid compute ID

E: Compute gyration/shape compute ID does not point to a gyration compute

Self-explanatory.  Provide an ID of a compute gyration command.
*/
