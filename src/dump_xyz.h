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
DumpStyle(xyz,DumpXYZ);
// clang-format on
#else

#ifndef LMP_DUMP_XYZ_H
#define LMP_DUMP_XYZ_H

#include "dump.h"

namespace LAMMPS_NS {

class DumpXYZ : public Dump {
 public:
  DumpXYZ(class LAMMPS *, int, char **);
  ~DumpXYZ() override;

 protected:
  int ntypes;
  char **typenames;

  void init_style() override;
  void write_header(bigint) override;
  void pack(tagint *) override;
  int convert_string(int, double *) override;
  void write_data(int, double *) override;
  int modify_param(int, char **) override;

  typedef void (DumpXYZ::*FnPtrWrite)(int, double *);
  FnPtrWrite write_choice;    // ptr to write data functions
  void write_string(int, double *);
  void write_lines(int, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif
