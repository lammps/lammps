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

#ifndef LMP_FILE_WRITER_WRAPPER_H
#define LMP_FILE_WRITER_WRAPPER_H

#include <string>
#include <stdio.h>

#include "file_writer.h"

namespace LAMMPS_NS {

class FileWriterWrapper : public FileWriter {
 public:
  FileWriterWrapper() = delete;
  ~FileWriterWrapper() = default;

  FileWriterWrapper(FILE *file_ptr) { fp = file_ptr; }

  size_t write(const void *buffer, size_t length) final {
    return fwrite(buffer, length, 1, fp) * length;
  }

  void open(const std::string &path, bool append = false) final {
    throw FileWriterException("Cannot call open on FileWriterWrapper");
  }
  void close() final {
    throw FileWriterException("Cannot call close on FileWriterWrapper");
  }
  void flush() final { fflush(fp); }
  [[nodiscard]] bool isopen() const final { return true; }

 private:
  FILE* fp;
};

} // namespace LAMMPS_NS
#endif
