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

#include "lmptype.h"
#include "pointers.h"
#include "error.h"

#define BRERR {FLERR}

namespace LAMMPS_NS {

class BufferReader {
 public:
  BufferReader() : BufferReader(nullptr, 0, nullptr) {}
  BufferReader(char *b, Error* e) : BufferReader(b, MAXBIGINT, e) {}
  BufferReader(char *b, bigint l, Error* e) : error(e), buffer(b), length(l) {}
  BufferReader(std::vector<char>&& v, Error* e)
    : error(e), owned_bytes(std::move(v)), buffer(owned_bytes.data()),
      length(owned_bytes.size()) {}

  struct ErrInfo {
    ErrInfo() = delete;
    ErrInfo(const char* f, int l) : file(f), line(l) {}
    const char* file;
    int line;
  };

  // Read a sub-buffer of length len
  virtual BufferReader sub_buf(const ErrInfo& e, bigint len);

  virtual FILE* get_fp() const { return nullptr; }
  virtual void seek(const ErrInfo& e, bigint p);

  void seek_offset(const ErrInfo& e, bigint o) { this->seek(e, pos() + o); }

  bigint size() const { return length; }
  bigint pos() const { return position; }
  bigint remaining() const { return length - position; }
  explicit operator bool() const { return get_fp() || buffer; }

 protected:
  BufferReader(Error* e) : error(e), length(0), buffer(nullptr) {}

  // Reads directly into destination
  virtual void read_raw_buf(const ErrInfo& e, char* dst, bigint count);

  // Returns an ephemeral buffer that can only be used immediately
  virtual char* get_raw_buf(const ErrInfo& e, bigint count);

  template<typename T>
  T* get_buf(const ErrInfo& e, bigint count) {
    static_assert(std::is_trivially_copyable_v<T>);
    return (T*)(this->get_raw_buf(e, count*sizeof(T)));
  }

 public:
  template<typename T>
  void read_buf(const ErrInfo& e, T* t, bigint count) {
    if (count < 0) error->one(e.file, e.line, "Illegal buffer size {}!", count);
    if constexpr (std::is_trivially_copyable_v<T>) {
      this->read_raw_buf(e, (char*)t, count*sizeof(T));
    } else {
      for (int i = 0; i < count; i++) t[i] = read<T>(e);
    }
  }

  // Caller responsible for deallocating returned array
  template<typename T>
  T* read_buf(const ErrInfo& e, bigint count) {
    if (count < 0) error->one(e.file, e.line, "Illegal buffer size {}!", count);
    if (count == 0) return nullptr;
    T* t = new T[count];
    read_buf(e, t, count);
    return t;
  }

  template<typename T>
  void read(const ErrInfo& e, T& t) { read_buf(e, &t, 1); }

  // This is the function to override for custom type handling, all others
  // eventually just invoke this for non-trivially-copyable types
  template<typename T>
  T read(const ErrInfo& e) { return *get_buf<T>(e, 1); }

 protected:
  Error* error;
  std::vector<char> owned_bytes;

  char *buffer;
  bigint length;
  bigint position = 0;
};

template<> std::string BufferReader::read<std::string>(const ErrInfo& e);
// Caller responsible for deallocating array
template<> char* BufferReader::read<char*>(const ErrInfo& e);

}    // namespace LAMMPS_NS

#endif
