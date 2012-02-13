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

DihedralStyle(charmm,DihedralCharmm)

#else

#ifndef LMP_DIHEDRAL_CHARMM_H
#define LMP_DIHEDRAL_CHARMM_H

#include "stdio.h"
#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralCharmm : public Dihedral {
 public:
  DihedralCharmm(class LAMMPS *);
  virtual ~DihedralCharmm();
  virtual void compute(int, int);
  void coeff(int, char **);
  void init_style();
  void write_restart(FILE *);
  void read_restart(FILE *);

 protected:
  double *k,*weight,*cos_shift,*sin_shift;
  int *multiplicity,*shift;
  double **lj14_1,**lj14_2,**lj14_3,**lj14_4;
  int implicit,weightflag;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

W: Dihedral problem: %d %ld %d %d %d %d

Conformation of the 4 listed dihedral atoms is extreme; you may want
to check your simulation geometry.

E: Incorrect args for dihedral coefficients

Self-explanatory.  Check the input script or data file.

E: Incorrect multiplicity arg for dihedral coefficients

Self-explanatory.  Check the input script or data file.

E: Incorrect weight arg for dihedral coefficients

Self-explanatory.  Check the input script or data file.

E: Dihedral charmm is incompatible with Pair style

Dihedral style charmm must be used with a pair style charmm
in order for the 1-4 epsilon/sigma parameters to be defined.

*/
