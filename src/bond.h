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

#ifndef BOND_H
#define BOND_H

#include "stdio.h"
#include "pointers.h"

namespace LAMMPS_NS {

class Bond : protected Pointers {
 public:
  int allocated;
  int *setflag;
  double energy;
  double eng_vdwl;
  double virial[6];

  Bond(class LAMMPS *);
  virtual ~Bond() {}
  virtual void init();
  virtual void init_style() {}

  virtual void compute(int, int) = 0;
  virtual void settings(int, char **) {}
  virtual void coeff(int, char **) = 0;
  virtual double equilibrium_distance(int) = 0;
  virtual void write_restart(FILE *) = 0;
  virtual void read_restart(FILE *) = 0;
  virtual void single(int, double, int, int, int, double &, double &) = 0;
  virtual int memory_usage() {return 0;}
};

}

#endif
