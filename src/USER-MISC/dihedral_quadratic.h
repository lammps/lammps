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

#ifdef DIHEDRAL_CLASS

DihedralStyle(quadratic,DihedralQuadratic)

#else

#ifndef LMP_DIHEDRAL_QUADRATIC_H
#define LMP_DIHEDRAL_QUADRATIC_H

#include "stdio.h"
#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralQuadratic : public Dihedral {
 public:
  DihedralQuadratic(class LAMMPS *);
  ~DihedralQuadratic();
  void compute(int, int);
  void coeff(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);

 protected:
  double *k,*phi0;
  int *sign,*multiplicity;

  void allocate();
};

}

#endif
#endif
