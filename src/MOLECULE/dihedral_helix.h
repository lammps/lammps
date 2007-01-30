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

#ifndef DIHEDRAL_HELIX_H
#define DIHEDRAL_HELIX_H

#include "stdio.h"
#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralHelix : public Dihedral {
 public:
  DihedralHelix(class LAMMPS *);
  ~DihedralHelix();
  void compute(int, int);
  void coeff(int, int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);

 private:
  double *aphi,*bphi,*cphi;

  void allocate();
};

}

#endif
