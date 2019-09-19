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

DihedralStyle(spherical,DihedralSpherical)

#else

#ifndef LMP_DIHEDRAL_SPHERICAL_H
#define LMP_DIHEDRAL_SPHERICAL_H

#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralSpherical : public Dihedral {
 public:
  DihedralSpherical(class LAMMPS *);
  virtual ~DihedralSpherical();
  virtual void compute(int, int);
  double CalcGeneralizedForces(int, double, double, double,
                               double*, double*, double*);
  void coeff(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);

 protected:
  int    *nterms;
  double **Ccoeff;
  double **phi_mult;
  double **phi_shift;
  double **phi_offset;
  double **theta1_mult;
  double **theta1_shift;
  double **theta1_offset;
  double **theta2_mult;
  double **theta2_shift;
  double **theta2_offset;

  void allocate();
};

}

#endif
#endif
