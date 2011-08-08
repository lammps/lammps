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
  The angle is defined from Cos(theta)=  b1.b2 / |b1||b2|  where b1=(R1-R2) and
  b2=(R3-R2) a straight bond has costheta=180.

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

#ifdef ANGLE_CLASS
AngleStyle(cosineshiftexp,AngleCosineShiftExp)
#else

#ifndef LMP_ANGLE_COSINESHIFTEXP_H
#define LMP_ANGLE_COSINESHIFTEXP_H

#include "stdio.h"
#include "angle.h"

namespace LAMMPS_NS {

class AngleCosineShiftExp : public Angle {
 public:
  AngleCosineShiftExp(class LAMMPS *);
  ~AngleCosineShiftExp();
  void compute(int, int);
  void coeff(int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, int, int, int);

 private:
  bool *doExpansion;
  double *umin,*a,*opt1;
  double *theta0;
  double *sint;
  double *cost;

  void allocate();
};

}

#endif
#endif
