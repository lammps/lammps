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

#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "pointers.h"

namespace LAMMPS_NS {

class Integrate : protected Pointers {
 public:
  Integrate(class LAMMPS *lmp, int, char **) : Pointers(lmp) {}
  virtual ~Integrate() {}
  virtual void init() = 0;
  virtual void setup() = 0;
  virtual void iterate(int) = 0;
  virtual void cleanup() {}
  virtual double memory_usage() {return 0.0;}
};

}

#endif
