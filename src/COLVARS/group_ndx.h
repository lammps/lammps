// -*- c++ -*-

/* ----------------------------------------------------------
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
CommandStyle(group2ndx,Group2Ndx);
// clang-format on
#else

#ifndef LMP_GROUP_NDX_H
#define LMP_GROUP_NDX_H

#include "command.h"

namespace LAMMPS_NS {

class Group2Ndx : public Command {
 public:
  Group2Ndx(class LAMMPS *lmp) : Command(lmp){};
  void command(int, char **) override;

 private:
  void write_group(FILE *, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
