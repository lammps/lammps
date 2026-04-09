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

#ifndef LMP_FILE_WRITER_ROOT_WRAPPER_H
#define LMP_FILE_WRITER_ROOT_WRAPPER_H

#include <string>
#include <stdio.h>

#include "file_writer.h"

namespace LAMMPS_NS {

class FileWriterRootWrapper : public FileWriter {
 public:
  ~FileWriterRootWrapper() = default;
  FileWriterRootWrapper() : FileWriterRootWrapper(-1, nullptr) {}
  FileWriterRootWrapper(int m_rank, FILE *file) : rank(m_rank), fp(file) {}
  
  [[nodiscard]] FILE* get_fp() const override { return fp; }

  size_t write(const void *buffer, size_t length) final {
    if (rank == 0) return fwrite(buffer, length, 1, fp) * length;
    else return 0;
  }

  void open(const std::string &path, bool append = false) final {
    throw FileWriterException("Cannot call open on FileWriterRootWrapper");
  }
  void close() final {
    throw FileWriterException("Cannot call close on FileWriterRootWrapper");
  }
  void flush() final { if (rank == 0) fflush(fp); }
  [[nodiscard]] bool isopen() const final { return true; }

 private:
  int rank;
  FILE* fp;
};

} // namespace LAMMPS_NS
#endif
