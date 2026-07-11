/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_COMMAND_H
#define LMP_COMMAND_H

#include "pointers.h"    // IWYU pragma: keep

namespace LAMMPS_NS {

class Command : protected Pointers {
 public:
  Command(class LAMMPS *lmp) : Pointers(lmp), copymode(0) {};
  virtual void command(int, char **) = 0;

 protected:
  int copymode;    // if set, do not deallocate during destruction
                   // required when classes are used as functors by Kokkos
};

}    // namespace LAMMPS_NS

#endif
