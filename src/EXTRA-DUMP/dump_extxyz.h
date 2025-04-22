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

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(extxyz,DumpExtXYZ);
// clang-format on
#else

#ifndef LMP_DUMP_EXTXYZ_H
#define LMP_DUMP_EXTXYZ_H

#include "dump_xyz.h"

namespace LAMMPS_NS {
class DumpExtXYZ : public DumpXYZ {
 public:
  DumpExtXYZ(class LAMMPS *, int, char **);
  ~DumpExtXYZ() override;

 protected:
  int with_vel;
  int with_forces;
  int with_mass;
  int with_pe;
  int with_temp;
  int with_press;
  char *properties_string;

  void update_properties();
  void init_style() override;
  void write_header(bigint) override;
  void pack(tagint *) override;
  int convert_string(int, double *) override;
  int modify_param(int, char **) override;

  void write_lines(int, double *) override;
};
}    // namespace LAMMPS_NS
#endif
#endif
