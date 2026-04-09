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

#include "buffer_reader.h"
#include "error.h"
#include "fmt/format.h"
#include <cstring>


namespace LAMMPS_NS {

void BufferReader::read_raw_buf(const ErrInfo& e, char* dst, bigint count) {
  memcpy(dst, this->get_raw_buf(e, count), count);
}

char* BufferReader::get_raw_buf(const ErrInfo& e, bigint count) {
  if (count < 0) error->one(e.file, e.line, "Attempt to read {} bytes", count);
  if (count > remaining()) error->one(e.file, e.line,
    "Attempt to read {} bytes with only {} remaining", count, remaining()
  );
  char* ret = buffer+position;
  position += count;
  return ret;
}

BufferReader BufferReader::sub_buf(const ErrInfo& e, bigint len) {
  if (len < 0) error->one(e.file, e.line, "Illegal sub_buf size {}!", len);
  if (len > remaining()) len = remaining();
  return BufferReader(get_buf<char>(e, len), len, error);
}

void BufferReader::seek(const ErrInfo& e, bigint p) {
  if (p < 0 || p > length) error->one(
    e.file, e.line, "Illegal seek to {} on {} sized buffer", p, this->size()
  );
  position = p;
}

template<> std::string BufferReader::read<std::string>(const ErrInfo& e) {
  int size = read<int>(e);
  if (size < 0) error->one(e.file, e.line, "Illegal string size {}!", size);
  return std::string(get_buf<char>(e, size), size-1);
}

template<> char* BufferReader::read<char*>(const ErrInfo& e) {
  int size = read<int>(e);
  if (size < 0) error->one(e.file, e.line, "Illegal string size {}!", size);
  return size == 0 ? nullptr : read_buf<char>(e, size);
}

}
