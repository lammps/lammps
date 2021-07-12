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
DumpStyle(cfg/zstd,DumpCFGZstd);
// clang-format on
#else

#ifndef LMP_DUMP_CFG_ZSTD_H
#define LMP_DUMP_CFG_ZSTD_H

#include "dump_cfg.h"
#include "zstd_file_writer.h"

namespace LAMMPS_NS {

class DumpCFGZstd : public DumpCFG {
 public:
  DumpCFGZstd(class LAMMPS *, int, char **);
  virtual ~DumpCFGZstd();

 protected:
  ZstdFileWriter writer;

  virtual void openfile();
  virtual void write_header(bigint);
  virtual void write_data(int, double *);
  virtual void write();

  virtual int modify_param(int, char **);
};

}    // namespace LAMMPS_NS

#endif
#endif
#endif

/* ERROR/WARNING messages:

E: Dump cfg/zstd only writes compressed files

The dump cfg/zstd output file name must have a .zstd suffix.

E: Cannot open dump file

Self-explanatory.

*/
