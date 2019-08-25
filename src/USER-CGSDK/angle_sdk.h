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

AngleStyle(sdk,AngleSDK)

#else

#ifndef LMP_ANGLE_SDK_H
#define LMP_ANGLE_SDK_H

#include <cstdio>
#include "angle.h"

namespace LAMMPS_NS {

class AngleSDK : public Angle {
 public:
  AngleSDK(class LAMMPS *);
  virtual ~AngleSDK();
  virtual void compute(int, int);
  void coeff(int, char **);
  void init_style();
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);
  double single(int, int, int, int);

 protected:
  double *k,*theta0;

  // scaling factor for repulsive 1-3 interaction
  double *repscale;
  // parameters from SDK pair style
  int **lj_type;
  double **lj1,**lj2, **lj3, **lj4;
  double **rminsq,**emin;

  int repflag; // 1 if we have to handle 1-3 repulsion

  void ev_tally13(int, int, int, int, double, double,
                  double, double, double);

  void allocate();
};

}

#endif
#endif
