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

ComputeStyle(inertia/molecule,ComputeInertiaMolecule)

#else

#ifndef LMP_COMPUTE_INERTIA_MOLECULE_H
#define LMP_COMPUTE_INERTIA_MOLECULE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeInertiaMolecule : public Compute {
 public:
  ComputeInertiaMolecule(class LAMMPS *, int, char **);
  ~ComputeInertiaMolecule();
  void init();
  void compute_array();
  double memory_usage();

 private:
  int nmolecules;
  int idlo,idhi;

  double *massproc,*masstotal;
  double **com,**comall;
  double **inertia,**inertiaall;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute inertia/molecule requires molecular atom style

Self-explanatory.

E: Molecule count changed in compute inertia/molecule

Number of molecules must remain constant over time.

*/
