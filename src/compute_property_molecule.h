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

ComputeStyle(property/molecule,ComputePropertyMolecule)

#else

#ifndef LMP_COMPUTE_PROPERTY_MOLECULE_H
#define LMP_COMPUTE_PROPERTY_MOLECULE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePropertyMolecule : public Compute {
 public:
  ComputePropertyMolecule(class LAMMPS *, int, char **);
  ~ComputePropertyMolecule();
  void init();
  void compute_vector();
  void compute_array();
  double memory_usage();

 private:
  int nvalues,nmolecules;
  tagint idlo,idhi;

  double *buf;

  typedef void (ComputePropertyMolecule::*FnPtrPack)(int);
  FnPtrPack *pack_choice;              // ptrs to pack functions

  void pack_mol(int);
  void pack_count(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute property/molecule requires molecular atom style

Self-explanatory.

E: Invalid keyword in compute property/molecule command

Self-explanatory.

E: Molecule count changed in compute property/molecule

Number of molecules must remain constant over time.

*/
