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

DihedralStyle(opls,DihedralOPLS)

#else

#ifndef LMP_DIHEDRAL_OPLS_H
#define LMP_DIHEDRAL_OPLS_H

#include "stdio.h"
#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralOPLS : public Dihedral {
 public:
  DihedralOPLS(class LAMMPS *);
  virtual ~DihedralOPLS();
  virtual void compute(int, int);
  void coeff(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);

 protected:
  double *k1,*k2,*k3,*k4;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

W: Dihedral problem: %d %ld    

UNDOCUMENTED

E: Incorrect args for dihedral coefficients

Self-explanatory.  Check the input script or data file.

U: Dihedral problem: %d %ld %d %d %d %d

Conformation of the 4 listed dihedral atoms is extreme; you may want
to check your simulation geometry.

*/
