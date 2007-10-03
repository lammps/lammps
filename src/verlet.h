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

#ifndef VERLET_H
#define VERLET_H

#include "integrate.h"

namespace LAMMPS_NS {

class Verlet : public Integrate {
 public:
  Verlet(class LAMMPS *, int, char **);
  ~Verlet();
  void init();
  void setup();
  void iterate(int);

 private:
  int virial_style;                 // compute virial explicitly or implicitly
  int virial_every;                 // 1 if virial computed every step
  int next_virial;                  // next timestep to compute virial
  int nfix_virial;                  // # of fixes that need virial occasionally
  int *fix_virial_every;            // frequency they require it
  int *next_fix_virial;             // next timestep they need it
  int triclinic;                    // 0 if domain is orthog, 1 if triclinic

  int torqueflag;                   // arrays to zero out every step

  void force_clear(int);
  int fix_virial(int);
};

}

#endif
