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

#ifdef COMMAND_CLASS

CommandStyle(DEPRECATED,Deprecated)

#else

#ifndef LMP_DEPRECATED_H
#define LMP_DEPRECATED_H

#include "pointers.h"

namespace LAMMPS_NS {

class Deprecated : protected Pointers {
 public:
  Deprecated(class LAMMPS *lmp) : Pointers(lmp) {};
  void command(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

W: Ignoring unknown or incorrect info command flag

Self-explanatory.  An unknown argument was given to the info command.
Compare your input with the documentation.

E: Unknown name for info package category

Self-explanatory.

E: Unknown name for info newton category

Self-explanatory.

E: Unknown name for info pair category

Self-explanatory.

E: Unknown category for info is_active()

Self-explanatory.

E: Unknown category for info is_available()

Self-explanatory.

E: Unknown category for info is_defined()

Self-explanatory.

*/
