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

DihedralStyle(fourier,DihedralFourier)

#else

#ifndef LMP_DIHEDRAL_FOURIER_H
#define LMP_DIHEDRAL_FOURIER_H

#include "stdio.h"
#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralFourier : public Dihedral {
 public:
  DihedralFourier(class LAMMPS *);
  virtual ~DihedralFourier();
  virtual void compute(int, int);
  void coeff(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);

 protected:
  double **k,**cos_shift,**sin_shift,**shift;
  int **multiplicity;
  int *nterms;
  int implicit,weightflag;

  void allocate();
};

}

#endif
#endif
