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

#ifndef LMP_ZSTD_FILE_WRITER_H
#define LMP_ZSTD_FILE_WRITER_H

#include "file_writer.h"

#include <string>
#include <zstd.h>

#if ZSTD_VERSION_NUMBER < 10400
#error must have at least zstd version 1.4 to compile with -DLAMMPS_ZSTD
#endif


namespace LAMMPS_NS {

class ZstdFileWriter : public FileWriter {
  int compression_level;
  int checksum_flag;

  ZSTD_CCtx *cctx;
  FILE *fp;
  char *out_buffer;
  size_t out_buffer_size;

 public:
  ZstdFileWriter();
  virtual ~ZstdFileWriter();
  virtual void open(const std::string &path, bool append = false) override;
  virtual void close() override;
  virtual void flush() override;
  virtual size_t write(const void *buffer, size_t length) override;
  virtual bool isopen() const override;

  void setCompressionLevel(int level);
  void setChecksum(bool enabled);
};
}    // namespace LAMMPS_NS

#endif
#endif
