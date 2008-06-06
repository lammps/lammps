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

#ifndef DIHEDRAL_CHARMM_H
#define DIHEDRAL_CHARMM_H

#include "stdio.h"
#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralCharmm : public Dihedral {
 public:
  DihedralCharmm(class LAMMPS *);
  ~DihedralCharmm();
  void compute(int, int);
  void coeff(int, int, char **);
  void init_style();
  void write_restart(FILE *);
  void read_restart(FILE *);

 private:
  double *k,*weight,*cos_shift,*sin_shift;
  int *multiplicity,*shift;
  double **lj14_1,**lj14_2,**lj14_3,**lj14_4;
  int implicit,weightflag;

  void allocate();
};

}

#endif
