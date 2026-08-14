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
CommandStyle(region2vmd,Region2VMD);
// clang-format on
#else

#ifndef LMP_REGION2VMD_H
#define LMP_REGION2VMD_H

#include "command.h"

namespace LAMMPS_NS {

class Region;

class Region2VMD : public Command {
 public:
  Region2VMD(class LAMMPS *lmp) : Command(lmp) {};
  void command(int, char **) override;

 protected:
  void write_region(FILE *, Region *);
};

}    // namespace LAMMPS_NS

#endif
#endif
