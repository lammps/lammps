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

/* ----------------------------------------------------------------------
   Contributing author: Matthew Whitlock (Sandia)
------------------------------------------------------------------------- */

#ifndef LMP_BUFFER_READER_H
#define LMP_BUFFER_READER_H

#include <type_traits>
#include <string>
#include <exception>

#include "lmptype.h"
#include "fmt/format.h"

namespace LAMMPS_NS {

class BufferReaderException : public std::exception {
  std::string message;

 public:
  BufferReaderException(const std::string &msg) : message(msg) {}

  [[nodiscard]] const char *what() const noexcept override {
    return message.c_str();
  }
};

class BufferReader {
 public:
  BufferReader() : BufferReader(nullptr, 0) {}
  BufferReader(char *b) : BufferReader(b, MAXBIGINT) {}
  BufferReader(const BufferReader& b, bigint l) : BufferReader(b.buffer, l) {}
  BufferReader(char *b, bigint l) : buffer(b), length(l) {}

  template<typename T>
  T read() {
    static_assert(std::is_trivially_copyable_v<T>);
    if (sizeof(T) > length) {
      throw BufferReaderException(fmt::format(
        "Attempt to read type of size {} with only {} length", sizeof(T), length
      ));
    }
    length -= sizeof(T);
    buffer += sizeof(T);
    return *((T*)(buffer - sizeof(T)));
  }

  template<typename T>
  void read(T& t) { t = read<T>(); }

  // Read a sub-buffer of length len
  BufferReader sub_buf(bigint len) {
    if (len > length) len = length;
    BufferReader ret(buffer, len);
    buffer += len;
    length -= len;
    return ret;
  }

  void skip_bytes(bigint len) {
    if (len > length) len = length;
    buffer += len;
    length -= len;
  }

  char *buffer;
  bigint length;
};

}    // namespace LAMMPS_NS

#endif
