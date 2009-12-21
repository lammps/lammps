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

#ifndef COMPUTE_COM_MOLECULE_H
#define COMPUTE_COM_MOLECULE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeCOMMolecule : public Compute {
 public:
  ComputeCOMMolecule(class LAMMPS *, int, char **);
  ~ComputeCOMMolecule();
  void init();
  void compute_array();
  double memory_usage();

 private:
  int nmolecules;
  int idlo,idhi;

  double *massproc,*masstotal;
  double **com,**comall;
};

}

#endif
