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
  The angle is defined from  (R1-R2) . (R3-R2) hence a straight
  bond has costheta=180.

  U(theta,theta0,umin)= -Umin (1+cos[theta-theta0])/2

  potential has minimum at theta=theta0 where U() = -Umin
  potential has maximum at theta=theta0+180 where U() = 0

  At the minimum U(theta,theta0,umin)=-Umin+Umin/4 (theta-theta0)^2+O( ()^4)
  hence the effective spring constant is k=umin/2.

  to match   U(theta,theta0,K)=K(theta-theta0)^2 at the minimum Umin=4*K
*/

#ifdef ANGLE_CLASS

AngleStyle(cosineshift,AngleCosineShift)

#else

#ifndef LMP_ANGLE_COSINESHIFT_H
#define LMP_ANGLE_COSINESHIFT_H

#include "stdio.h"
#include "angle.h"

namespace LAMMPS_NS {

class AngleCosineShift : public Angle {
 public:
  AngleCosineShift(class LAMMPS *);
  ~AngleCosineShift();
  virtual void compute(int, int);
  void coeff(int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, int, int, int);

 private:
  double *k;
  double *a;
  double *theta;
  double *ksint;
  double *kcost;

  void allocate();
};

}

#endif
#endif
