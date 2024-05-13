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

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(write_data,WriteData);
// clang-format on
#else

#ifndef LMP_WRITE_DATA_H
#define LMP_WRITE_DATA_H

#include "command.h"

namespace LAMMPS_NS {

class WriteData : public Command {
 public:
  WriteData(class LAMMPS *);
  void command(int, char **) override;
  void write(const std::string &);

 private:
  int me, nprocs;
  int pairflag;
  int coeffflag;
  int fixflag;
  int triclinic_general;
  int lmapflag;
  FILE *fp;
  bigint nbonds_local, nbonds;
  bigint nangles_local, nangles;
  bigint ndihedrals_local, ndihedrals;
  bigint nimpropers_local, nimpropers;

  void header();
  void type_arrays();
  void force_fields();
  void atoms();
  void velocities();
  void bonds();
  void angles();
  void dihedrals();
  void impropers();
  void bonus(int);
  void fix(class Fix *, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
