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
DumpStyle(xyz/gz,DumpXYZGZ);
// clang-format on
#else

#ifndef LMP_DUMP_XYZ_GZ_H
#define LMP_DUMP_XYZ_GZ_H

#include "dump_xyz.h"
#include "gz_file_writer.h"

namespace LAMMPS_NS {

class DumpXYZGZ : public DumpXYZ {
 public:
  DumpXYZGZ(class LAMMPS *, int, char **);

 protected:
  GzFileWriter writer;

  void openfile() override;
  void write_header(bigint) override;
  void write_data(int, double *) override;
  void write() override;

  int modify_param(int, char **) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
