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

#ifdef DIHEDRAL_CLASS

DihedralStyle(multi/harmonic,DihedralMultiHarmonic)

#else

#ifndef LMP_DIHEDRAL_MULTI_HARMONIC_H
#define LMP_DIHEDRAL_MULTI_HARMONIC_H

#include "stdio.h"
#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralMultiHarmonic : public Dihedral {
 public:
  DihedralMultiHarmonic(class LAMMPS *);
  virtual ~DihedralMultiHarmonic();
  virtual void compute(int, int);
  void coeff(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);

 protected:
  double *a1,*a2,*a3,*a4,*a5;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

W: Dihedral problem: %d %ld %ld %ld %ld %ld

Conformation of the 4 listed dihedral atoms is extreme; you may want
to check your simulation geometry.

E: Incorrect args for dihedral coefficients

Self-explanatory.  Check the input script or data file.

*/
