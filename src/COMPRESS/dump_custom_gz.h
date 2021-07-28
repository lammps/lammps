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

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(custom/gz,DumpCustomGZ);
// clang-format on
#else

#ifndef LMP_DUMP_CUSTOM_GZ_H
#define LMP_DUMP_CUSTOM_GZ_H

#include "dump_custom.h"
#include "gz_file_writer.h"

namespace LAMMPS_NS {

class DumpCustomGZ : public DumpCustom {
 public:
  DumpCustomGZ(class LAMMPS *, int, char **);
  virtual ~DumpCustomGZ();

 protected:
  GzFileWriter writer;

  virtual void openfile();
  virtual void write_header(bigint);
  virtual void write_data(int, double *);
  virtual void write();

  virtual int modify_param(int, char **);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Dump custom/gz only writes compressed files

The dump custom/gz output file name must have a .gz suffix.

E: Cannot open dump file

Self-explanatory.

*/
