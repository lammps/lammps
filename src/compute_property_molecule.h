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

#ifndef COMPUTE_PROPERTY_MOLECULE_H
#define COMPUTE_PROPERTY_MOLECULE_H

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
  int idlo,idhi;

  double *buf;

  typedef void (ComputePropertyMolecule::*FnPtrPack)(int);
  FnPtrPack *pack_choice;              // ptrs to pack functions

  void pack_mol(int);
};

}

#endif
