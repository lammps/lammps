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

#ifndef LMP_GZ_FILE_WRITER_H
#define LMP_GZ_FILE_WRITER_H

#include "file_writer.h"

#include <string>
#include <zlib.h>

namespace LAMMPS_NS {

class GzFileWriter : public FileWriter {
  int compression_level;

  gzFile gzFp;    // file pointer for the compressed output stream
 public:
  GzFileWriter();
  virtual ~GzFileWriter();
  virtual void open(const std::string &path, bool append = false) override;
  virtual void close() override;
  virtual void flush() override;
  virtual size_t write(const void *buffer, size_t length) override;
  virtual bool isopen() const override;

  void setCompressionLevel(int level);
};
}    // namespace LAMMPS_NS

#endif
