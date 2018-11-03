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

#ifdef FIX_CLASS

FixStyle(plumed,FixPlumed)

#else

#ifndef LMP_FIX_PLUMED_H
#define LMP_FIX_PLUMED_H

#include "fix.h"

// forward declaration
namespace PLMD {
  class Plumed;
}

namespace LAMMPS_NS {

class FixPlumed : public Fix {
 public:
  FixPlumed(class LAMMPS *, int, char **);
  ~FixPlumed();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();

 private:
// pointer to plumed object:
  PLMD::Plumed *p;
// number of atoms local to this process:
  int nlocal;
// array of atom indexes local to this process:
  int *gatindex;
// array of masses for local atoms:
  double *masses;
// array of charges for local atoms:
  double *charges;
// this is something to enable respa
  int nlevels_respa;
// output bias potential
  double bias;
// Compute for the energy
  class Compute *c_pe; 
// Compute for the pressure
  class Compute *c_press; 
// Flag to trigger calculation of the energy and virial
  int plumedNeedsEnergy;
// ID for potential energy and pressure compute
  char  *id_pe, *id_press;
};

};

#endif
#endif
