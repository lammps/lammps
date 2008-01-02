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
  double energy;                  // accumulated energies
  double virial[6];               // accumlated virial
  double *eatom,**vatom;          // accumulated per-atom energy/virial

  Bond(class LAMMPS *);
  virtual ~Bond();
  virtual void init();
  virtual void init_style() {}

  virtual void compute(int, int) = 0;
  virtual void settings(int, char **) {}
  virtual void coeff(int, char **) = 0;
  virtual double equilibrium_distance(int) = 0;
  virtual void write_restart(FILE *) = 0;
  virtual void read_restart(FILE *) = 0;
  virtual double single(int, double, int, int) = 0;
  virtual double memory_usage();

 protected:
  int evflag;
  int eflag_either,eflag_global,eflag_atom;
  int vflag_either,vflag_global,vflag_atom;
  int maxeatom,maxvatom;

  void ev_setup(int, int);
  void ev_tally(int, int, int, int, double, double, double, double, double);
};

}

#endif
