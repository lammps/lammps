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

#include "buffer_reader_file.h"
#include "error.h"
#include "utils.h"
#include "platform.h"
#include "fmt/format.h"
#include <cstring>

namespace LAMMPS_NS {

BufferReaderFile::BufferReaderFile(
  FILE* f, const std::string& fn, Error* e
) : BufferReader(e), fname(fn), fp(f) {
  if (fp) {
    bigint curpos = platform::ftell(fp);
    platform::fseek(fp, platform::END_OF_FILE);
    bigint endpos = platform::ftell(fp);
    platform::fseek(fp, curpos);

    length = endpos - curpos;
  }
}


void BufferReaderFile::read_raw_buf(const ErrInfo& e, char* dst, bigint count) {
  if (!count) return;
  if (count < 0) error->one(e.file, e.line, "Attempt to read {} bytes", count);
  if (count > remaining()) error->one(e.file, e.line,
    "Attempt to read {} bytes with only {} remaining", count, remaining()
  );
  position += count;
  utils::sfread(e.file, e.line, dst, 1, count, fp, fname_cstr, error);
}

BufferReader BufferReaderFile::sub_buf(const ErrInfo& e, bigint len) {
  if (len > remaining()) len = remaining();
  std::vector<char> bytes(len);
  this->read_raw_buf(e, bytes.data(), len);
  return BufferReader(std::move(bytes), error);
}

void BufferReaderFile::seek(const ErrInfo& e, bigint p) {
  if (p < 0 || p > this->size()) error->one(
    e.file, e.line, "Illegal seek to {} on {} sized file", p, this->size()
  );
  position = p;
  platform::fseek(fp, position);
}

char* BufferReaderFile::get_raw_buf(const ErrInfo& e, bigint count) {
  owned_bytes.resize(count);
  this->read_raw_buf(e, owned_bytes.data(), count);
  return owned_bytes.data();
}

}
