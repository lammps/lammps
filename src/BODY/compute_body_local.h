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

ComputeStyle(body/local,ComputeBodyLocal)

#else

#ifndef LMP_COMPUTE_BODY_LOCAL_H
#define LMP_COMPUTE_BODY_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeBodyLocal : public Compute {
 public:
  ComputeBodyLocal(class LAMMPS *, int, char **);
  ~ComputeBodyLocal();
  void init();
  void compute_local();
  double memory_usage();

 private:
  int nvalues;
  int *which,*index;

  int nmax;
  double *vector;
  double **array;

  class AtomVecBody *avec;
  class Body *bptr;

  int compute_body(int);
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

E: Compute body/local requires atom style body

Self-explanatory.

E: Invalid index in compute body/local command

Self-explanatory.

E: Invalid index for non-body particles in compute body/local command

Only indices 1,2,3 can be used for non-body particles.

*/
