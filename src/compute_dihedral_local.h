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

ComputeStyle(dihedral/local,ComputeDihedralLocal)

#else

#ifndef LMP_COMPUTE_DIHEDRAL_LOCAL_H
#define LMP_COMPUTE_DIHEDRAL_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeDihedralLocal : public Compute {
 public:
  ComputeDihedralLocal(class LAMMPS *, int, char **);
  ~ComputeDihedralLocal();
  void init();
  void compute_local();
  double memory_usage();

 private:
  int nvalues,pflag;
  int ncount;

  int nmax;
  double *vector;
  double **array;

  int compute_dihedrals(int);
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

E: Compute dihedral/local used when dihedrals are not allowed

The atom style does not support dihedrals.

E: Invalid keyword in compute dihedral/local command

Self-explanatory.

E: No dihedral style is defined for compute dihedral/local

Self-explanatory.

*/
