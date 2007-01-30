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

#ifndef IMPROPER_CVFF_H
#define IMPROPER_CVFF_H

#include "stdio.h"
#include "improper.h"

namespace LAMMPS_NS {

class ImproperCvff : public Improper {
 public:
  ImproperCvff(class LAMMPS *);
  ~ImproperCvff();
  void compute(int, int);
  void coeff(int, int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);

 private:
  double *k;
  int *sign,*multiplicity;

  void allocate();
};

}

#endif
