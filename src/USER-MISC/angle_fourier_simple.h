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

#ifdef ANGLE_CLASS

AngleStyle(fourier/simple,AngleFourierSimple)

#else

#ifndef ANGLE_FOURIER_SIMPLE_H
#define ANGLE_FOURIER_SIMPLE_H

#include <cstdio>
#include "angle.h"

namespace LAMMPS_NS {

class AngleFourierSimple : public Angle {
 public:
  AngleFourierSimple(class LAMMPS *);
  virtual ~AngleFourierSimple();
  virtual void compute(int, int);
  void coeff(int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);
  virtual double single(int, int, int, int);

 protected:
  double *k,*C,*N;

  void allocate();
};

}

#endif
#endif
