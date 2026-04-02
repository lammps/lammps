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

#ifndef LMP_FILE_WRITER_SIZER_H
#define LMP_FILE_WRITER_SIZER_H

#include <string>

#include "lmptype.h"
#include "file_writer.h"

namespace LAMMPS_NS {

class FileWriterSizer : public FileWriter {
 public:
  FileWriterSizer() { len = 0; }
  ~FileWriterSizer() = default;

  bigint size() { return len; }

  size_t write(const void *buffer, size_t length) final {
    len += length;
    return length;
  };
  size_t write_restart_global_size(const Restartable*) final {
    return writev(len);
  }
  size_t write_restart_local_size(const Restartable*) final {
    return writev(len);
  }

  void open(const std::string &path, bool append = false) final { }
  void close() final { len = 0; }
  void flush() final { }
  [[nodiscard]] bool isopen() const final { return true; }
  [[nodiscard]] bool issizer() const final { return true; }

 private:
  bigint len;
};

} // namespace LAMMPS_NS
#endif
