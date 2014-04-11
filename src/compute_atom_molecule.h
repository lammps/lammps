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

ComputeStyle(atom/molecule,ComputeAtomMolecule)

#else

#ifndef LMP_COMPUTE_ATOM_MOLECULE_H
#define LMP_COMPUTE_ATOM_MOLECULE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeAtomMolecule : public Compute {
 public:
  ComputeAtomMolecule(class LAMMPS *, int, char **);
  ~ComputeAtomMolecule();
  void init();
  void compute_vector();
  void compute_array();
  double memory_usage();

 private:
  int nvalues,nmolecules;
  tagint idlo,idhi;

  int *which,*argindex,*value2index;
  char **ids;

  int nstride,maxatom;
  double *vone;
  double **aone;
  double *scratch;
  double *peratom;

  void compute_one(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute atom/molecule requires molecular atom style

Self-explanatory.

E: Compute ID for compute atom/molecule does not exist

Self-explanatory.

E: Compute atom/molecule compute does not calculate per-atom values

Self-explanatory.

E: Compute atom/molecule compute does not calculate a per-atom vector

Self-explanatory.

E: Compute atom/molecule compute does not calculate a per-atom array

Self-explanatory.

E: Compute atom/molecule compute array is accessed out-of-range

Self-explanatory.

E: Fix ID for compute atom/molecule does not exist

Self-explanatory.

E: Compute atom/molecule fix does not calculate per-atom values

Self-explanatory.

E: Compute atom/molecule fix does not calculate a per-atom vector

Self-explanatory.

E: Compute atom/molecule fix does not calculate a per-atom array

Self-explanatory.

E: Compute atom/molecule fix array is accessed out-of-range

Self-explanatory.

E: Variable name for compute atom/molecule does not exist

Self-explanatory.

E: Compute atom/molecule variable is not atom-style variable

Self-explanatory.

E: Molecule count changed in compute atom/molecule

Number of molecules must remain constant over time.

E: Fix used in compute atom/molecule not computed at compatible time

The fix must produce per-atom quantities on timesteps that the compute
needs them.

*/
