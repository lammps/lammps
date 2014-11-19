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

ComputeStyle(force/molecule,ComputeForceMolecule)

#else

#ifndef LMP_COMPUTE_FORCE_MOLECULE_H
#define LMP_COMPUTE_FORCE_MOLECULE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeForceMolecule : public Compute {
 public:
  ComputeForceMolecule(class LAMMPS *, int, char **);
  ~ComputeForceMolecule();
  void init();
  void compute_array();
  double memory_usage();

 private:
  int nmolecules;
  tagint idlo,idhi;

  double **force,**forceall;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute force/molecule requires molecular atom style

Self-explanatory.

E: Molecule count changed in compute force/molecule

Number of molecules must remain constant over time.

*/
