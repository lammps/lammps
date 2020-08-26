/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(custom/zstd,DumpCustomZstd)

#else

#ifndef LMP_DUMP_CUSTOM_ZSTD_H
#define LMP_DUMP_CUSTOM_ZSTD_H

#include "dump_custom.h"
#include <zstd.h>
#include <stdio.h>

namespace LAMMPS_NS {

class DumpCustomZstd : public DumpCustom {
 public:
  DumpCustomZstd(class LAMMPS *, int, char **);
  virtual ~DumpCustomZstd();

 protected:
  int compression_level;
  int checksum_flag;

  ZSTD_CCtx * cctx;
  FILE * zstdFp;
  char * out_buffer;
  size_t out_buffer_size;

  virtual void openfile();
  virtual void write_header(bigint);
  virtual void write_data(int, double *);
  virtual void write();

  virtual int modify_param(int, char **);

  void zstd_write(const void * buffer, size_t length);
  void zstd_flush();
  void zstd_close();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Dump custom/zstd only writes compressed files

The dump custom/zstd output file name must have a .zst suffix.

E: Cannot open dump file

Self-explanatory.

*/
