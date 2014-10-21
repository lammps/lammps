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

ComputeStyle(ti,ComputeTI)

#else

#ifndef COMPUTE_TI_H
#define COMPUTE_TI_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTI : public Compute {
 public:
  ComputeTI(class LAMMPS *, int, char **);
  ~ComputeTI();
  void init();
  double compute_scalar();

 private:
  int nterms;
  int *which;
  int *ivar1,*ivar2;
  int *ilo, *ihi; 
  char **var1,**var2;
  class Pair **pptr;
  char **pstyle;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Variable name for compute ti does not exist

Self-explanatory.

E: Variable for compute ti is invalid style

Self-explanatory.

E: Compute ti pair style does not exist

Self-explanatory.

E: Compute ti tail when pair style does not compute tail corrections

Self-explanatory.

E: Compute ti kspace style does not exist

Self-explanatory.

E: Energy was not tallied on needed timestep

You are using a thermo keyword that requires potentials to
have tallied energy, but they didn't on this timestep.  See the
variable doc page for ideas on how to make this work.

*/
