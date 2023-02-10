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
// list all deprecated and removed dump styles here
DumpStyle(DEPRECATED,DumpDeprecated);
// clang-format on
#else

#ifndef LMP_DUMP_DEPRECATED_H
#define LMP_DUMP_DEPRECATED_H

#include "dump.h"

namespace LAMMPS_NS {

class DumpDeprecated : public Dump {
 public:
  DumpDeprecated(class LAMMPS *, int, char **);

  void init_style() override {}
  void write_header(bigint) override {}
  void pack(tagint *) override {}
  void write_data(int, double *) override {}
};

}    // namespace LAMMPS_NS

#endif
#endif
