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

ComputeStyle(rigid/local,ComputeRigidLocal)

#else

#ifndef LMP_COMPUTE_RIGID_LOCAL_H
#define LMP_COMPUTE_RIGID_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeRigidLocal : public Compute {
 public:
  ComputeRigidLocal(class LAMMPS *, int, char **);
  ~ComputeRigidLocal();
  void init();
  void compute_local();
  double memory_usage();

 private:
  int nvalues;
  int ncount;
  int *rstyle;

  char *idrigid;
  class FixRigidSmall *fixrigid;

  int nmax;
  double *vlocal;
  double **alocal;

  int compute_rigid(int);
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

E: Invalid keyword in compute rigid/local command

UNDOCUMENTED

E: FixRigidSmall ID for compute rigid/local does not exist

UNDOCUMENTED

E: Compute rigid/local does not use fix rigid/small fix

UNDOCUMENTED

U: Compute bond/local used when bonds are not allowed

The atom style does not support bonds.

U: Invalid keyword in compute bond/local command

Self-explanatory.

U: No bond style is defined for compute bond/local

Self-explanatory.

*/
