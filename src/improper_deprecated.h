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

#ifdef IMPROPER_CLASS

ImproperStyle(DEPRECATED,ImproperDeprecated)

#else

#ifndef LMP_IMPROPER_DEPRECATED_H
#define LMP_IMPROPER_DEPRECATED_H

#include "improper.h"

namespace LAMMPS_NS {

class ImproperDeprecated : public Improper {
 public:
  ImproperDeprecated(class LAMMPS *lmp) : Improper(lmp) {}
  virtual ~ImproperDeprecated() {}

  virtual void compute(int, int) {}
  virtual void settings(int, char **);
  virtual void coeff(int, char **) {}
  virtual void write_restart(FILE *) {}
  virtual void read_restart(FILE *) {}
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
