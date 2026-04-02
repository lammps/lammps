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

#ifndef LMP_FILE_WRITER_BUFFER_H
#define LMP_FILE_WRITER_BUFFER_H

#include <string>
#include "file_writer.h"

namespace LAMMPS_NS {

class FileWriterBuffer : public FileWriter {
 public:
  FileWriterBuffer() = delete;
  ~FileWriterBuffer() = default;

  FileWriterBuffer(char* buffer, size_t buffer_length) {
    if (buffer == nullptr && buffer_length > 0) {
      throw FileWriterException("Cannot write to null buffer");
    }
    buf = buffer;
    len = buffer_length;
  }
  FileWriterBuffer(char* b, bigint l) : FileWriterBuffer(b, (size_t)l) { }

  size_t write(const void *buffer, size_t length) final {
    if (length > len) {
      throw FileWriterException("Out of space in FileWriterBuffer");
    }
    memcpy(buf, buffer, length);
    buf += length;
    len -= length;
    return length;
  }

  void open(const std::string &path, bool append = false) final {
    throw FileWriterException("Cannot call open on FileWriterBuffer");
  }
  void close() final {
    buf = nullptr;
    len = 0;
  }
  void flush() final { }
  [[nodiscard]] bool isopen() const final { return buf != nullptr; }

 private:
  char* buf;
  size_t len;
};

} // namespace LAMMPS_NS
#endif
