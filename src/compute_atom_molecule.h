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
  int idlo,idhi;

  int *which,*argindex,*flavor,*value2index;
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
