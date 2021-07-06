/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_COMMAND_H
#define LMP_COMMAND_H

#include "pointers.h"

namespace LAMMPS_NS {

class Command : protected Pointers {
 public:
  Command(class LAMMPS *lmp) : Pointers(lmp){};
  virtual void command(int, char **) = 0;
};

}    // namespace LAMMPS_NS

#endif
