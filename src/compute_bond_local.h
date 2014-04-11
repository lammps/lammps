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

ComputeStyle(bond/local,ComputeBondLocal)

#else

#ifndef LMP_COMPUTE_BOND_LOCAL_H
#define LMP_COMPUTE_BOND_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeBondLocal : public Compute {
 public:
  ComputeBondLocal(class LAMMPS *, int, char **);
  ~ComputeBondLocal();
  void init();
  void compute_local();
  double memory_usage();

 private:
  int nvalues;
  int ncount;
  int *bstyle;
  int singleflag;

  int nmax;
  double *vector;
  double **array;

  int compute_bonds(int);
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

E: Compute bond/local used when bonds are not allowed

The atom style does not support bonds.

E: Invalid keyword in compute bond/local command

Self-explanatory.

E: No bond style is defined for compute bond/local

Self-explanatory.

*/
