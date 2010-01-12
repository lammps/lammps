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

#ifdef COMMAND_CLASS

CommandStyle(shell,Shell)

#else

#ifndef LMP_SHELL_H
#define LMP_SHELL_H

#include "pointers.h"

namespace LAMMPS_NS {

class Shell : protected Pointers {
 public:
  Shell(class LAMMPS *);
  void command(int, char **);
};

}

#endif
#endif
