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

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(create_bonds,CreateBonds);
// clang-format on
#else

#ifndef LMP_CREATE_BONDS_H
#define LMP_CREATE_BONDS_H

#include "command.h"

namespace LAMMPS_NS {

class CreateBonds : public Command {
 public:
  CreateBonds(class LAMMPS *);
  void command(int, char **) override;

 private:
  int igroup, group1bit, group2bit;
  int btype, atype, dtype;
  tagint batom1, batom2, aatom1, aatom2, aatom3, datom1, datom2, datom3, datom4;
  double rmin, rmax;

  void many();
  void single_bond();
  void single_angle();
  void single_dihedral();
  void single_improper();
};

}    // namespace LAMMPS_NS

#endif
#endif
