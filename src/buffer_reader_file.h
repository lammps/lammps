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

#ifndef LMP_BUFFER_READER_FILE_H
#define LMP_BUFFER_READER_FILE_H

#include "buffer_reader.h"
#include "lmptype.h"

namespace LAMMPS_NS {

class BufferReaderFile : public BufferReader {
 public:
  ~BufferReaderFile() = default;
  BufferReaderFile() : BufferReader() {}
  BufferReaderFile(FILE* f, Error* e)
    : BufferReaderFile(f, "", e) { }
  BufferReaderFile(FILE* f, const std::string& fn, Error* e);

  BufferReader sub_buf(const ErrInfo& e, bigint len) override;
  FILE* get_fp() const override { return fp; }
  void seek(const ErrInfo& e, bigint len) override;

  const std::string fname;

 protected:
  void read_raw_buf(const ErrInfo& e, char* b, bigint count) override;
  char* get_raw_buf(const ErrInfo& e, bigint count) override;

 private:
  FILE* fp;
  const char* fname_cstr = fname.empty() ? nullptr : fname.c_str();
};

}    // namespace LAMMPS_NS

#endif
