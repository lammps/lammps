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

#ifndef ANGLE_CHARMM_H
#define ANGLE_CHARMM_H

#include "stdio.h"
#include "angle.h"

namespace LAMMPS_NS {

class AngleCharmm : public Angle {
 public:
  AngleCharmm(class LAMMPS *);
  ~AngleCharmm();
  void compute(int, int);
  void coeff(int, int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, int, int, int);

 private:
  double *k,*theta0,*k_ub,*r_ub;

  void allocate();
};

}

#endif
