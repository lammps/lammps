/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributed by Timothy Sirk
------------------------------------------------------------------------- */

#ifdef READER_CLASS
// clang-format off
ReaderStyle(native,ReaderNative);
// clang-format on
#else

#ifndef LMP_READER_NATIVE_H
#define LMP_READER_NATIVE_H

#include "reader.h"

#include <map>
#include <string>

namespace LAMMPS_NS {

class ReaderNative : public Reader {
 public:
  ReaderNative(class LAMMPS *);
  ~ReaderNative() override;

  int read_time(bigint &) override;
  void skip() override;
  bigint read_header(double[3][3], int &, int &, int, int, int *, char **, int, int, int &, int &,
                     int &, int &) override;
  void read_atoms(int, int, double **) override;

 private:
  int revision;

  std::string magic_string;
  std::string unit_style;
  int *fieldindex;

  char *line;         // line read from dump file
  double *databuf;    // buffer for binary data
  int nwords;         // # of per-atom columns in dump file

  int size_one;       // number of double for one atom
  size_t maxbuf;      // maximum buffer size
  int nchunk;         // number of chunks in the binary file
  int ichunk;         // index of current reading chunk
  int natom_chunk;    // number of atoms in the current chunks
  int iatom_chunk;    // index of current atom in the current chunk

  int find_label(const std::string &label, const std::map<std::string, int> &labels);
  void read_lines(int);

  void read_buf(void *, size_t, size_t);
  void read_double_chunk(size_t);
  void skip_buf(size_t);
  void skip_reading_magic_str();
  bool is_known_magic_str() const;
  std::string read_binary_str(size_t);
};

}    // namespace LAMMPS_NS

#endif
#endif
