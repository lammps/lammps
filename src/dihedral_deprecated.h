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

#ifdef DIHEDRAL_CLASS

DihedralStyle(DEPRECATED,DihedralDeprecated)

#else

#ifndef LMP_DIHEDRAL_DEPRECATED_H
#define LMP_DIHEDRAL_DEPRECATED_H

#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralDeprecated : public Dihedral {
 public:
  DihedralDeprecated(class LAMMPS *lmp) : Dihedral(lmp) {}
  virtual ~DihedralDeprecated() {}

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
