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

#ifdef ANGLE_CLASS

AngleStyle(fourier,AngleFourier)

#else

#ifndef ANGLE_FOURIER_H
#define ANGLE_FOURIER_H

#include "stdio.h"
#include "angle.h"

namespace LAMMPS_NS {

class AngleFourier : public Angle {
 public:
  AngleFourier(class LAMMPS *);
  virtual ~AngleFourier();
  virtual void compute(int, int);
  void coeff(int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);
  virtual double single(int, int, int, int);

 protected:
  double *k,*C0,*C1,*C2;

  void allocate();
};

}

#endif
#endif
