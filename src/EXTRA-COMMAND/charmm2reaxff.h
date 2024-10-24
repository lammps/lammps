/* -*- c++ -*-----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(charmm2reaxff,Charmm2Reaxff);
// clang-format on
#else

#ifndef LMP_CHARMM2REAXFF_H
#define LMP_CHARMM2REAXFF_H

#include "command.h"

namespace LAMMPS_NS {

class Charmm2Reaxff : public Command {
 public:
  Charmm2Reaxff(class LAMMPS *lmp) : Command(lmp) {};
  void command(int, char **) override;

};

}    // namespace LAMMPS_NS

#endif
#endif
