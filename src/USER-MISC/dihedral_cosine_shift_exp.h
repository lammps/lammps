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

DihedralStyle(cosine/shift/exp,DihedralCosineShiftExp)

#else

#ifndef LMP_DIHEDRAL_COSINE_SHIFT_EXP_H
#define LMP_DIHEDRAL_COSINE_SHIFT_EXP_H

#include <cstdio>
#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralCosineShiftExp : public Dihedral {
 public:
  DihedralCosineShiftExp(class LAMMPS *);
  virtual ~DihedralCosineShiftExp();
  virtual void compute(int, int);
  void coeff(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);

 protected:
  bool *doExpansion;
  double *umin,*a,*opt1;
  double *sint;
  double *cost;
  double *theta;

  void allocate();
};

}

#endif
#endif
