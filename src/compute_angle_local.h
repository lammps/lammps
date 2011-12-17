/* ----------------------------------------------------------------------
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

ComputeStyle(angle/local,ComputeAngleLocal)

#else

#ifndef LMP_COMPUTE_ANGLE_LOCAL_H
#define LMP_COMPUTE_ANGLE_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeAngleLocal : public Compute {
 public:
  ComputeAngleLocal(class LAMMPS *, int, char **);
  ~ComputeAngleLocal();
  void init();
  void compute_local();
  double memory_usage();

 private:
  int nvalues,tflag,eflag;
  int ncount;

  int nmax;
  double *vector;
  double **array;

  int compute_angles(int);
  void reallocate(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute angle/local used when angles are not allowed

The atom style does not support angles.

E: Invalid keyword in compute angle/local command

Self-explanatory.

E: No angle style is defined for compute angle/local

Self-explanatory.

*/
