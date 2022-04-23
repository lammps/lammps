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

/* ----------------------------------------------------------------------
   Contributing author: Richard Berger (Temple U)
------------------------------------------------------------------------- */

#ifdef LAMMPS_ZSTD

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(custom/zstd,DumpCustomZstd);
// clang-format on
#else

#ifndef LMP_DUMP_CUSTOM_ZSTD_H
#define LMP_DUMP_CUSTOM_ZSTD_H

#include "dump_custom.h"
#include "zstd_file_writer.h"

namespace LAMMPS_NS {

class DumpCustomZstd : public DumpCustom {
 public:
  DumpCustomZstd(class LAMMPS *, int, char **);

 protected:
  ZstdFileWriter writer;

  void openfile() override;
  void write_header(bigint) override;
  void write_data(int, double *) override;
  void write() override;

  int modify_param(int, char **) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
#endif
