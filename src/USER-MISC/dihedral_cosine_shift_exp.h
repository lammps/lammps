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

/*
  The torsion angle is defined such that a straight bond has costheta=180.

  U(theta,theta0,umin,a) = -Umin[Exp(-a U)-1]/[[Exp(a)-1]]
  with U = (-1-cos[theta-theta0])/2 = -Cos((theta-theta0)/2)^2

  potential has minimum at theta=theta0 where U() = -Umin
  potential has maximum at theta=theta0+180 where U() = 0

  The spring constant around the minimum is controlled by a and
  is given by  k = a exp(a) Umin/[ 2(Exp[a]-1) ]   for a=0
  the spring constant is k=Umin/2 and the potential reduces to
  the cosineshifted potential.
  
  The potential is implemented such that for a<0.001 a series
  expansion to linear order is used instead of the expression
  above. This ensures a precision of about 1e-5 or better for
  energies and forces, and ensures the potential is well
  behaved for a=0
*/

#ifdef DIHEDRAL_CLASS

DihedralStyle(cosineshiftexp,DihedralCosShiftExp)

#else

#ifndef LMP_DIHEDRAL_COSINESHIFTEDEXP_H
#define LMP_DIHEDRAL_COSINESHIFTEDEXP_H

#include "stdio.h"
#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralCosShiftExp : public Dihedral {
 public:
  DihedralCosShiftExp(class LAMMPS *);
  ~DihedralCosShiftExp();
  void compute(int, int);
  void coeff(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);

 private:
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
