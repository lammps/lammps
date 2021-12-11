/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributed by Lixin Sun
------------------------------------------------------------------------- */

#ifdef READER_CLASS
// clang-format off
ReaderStyle(native_bin,ReaderNative);
// clang-format on
#else

#ifndef LMP_READER_NATIVE_BIN_H
#define LMP_READER_NATIVE_BIN_H

#include "reader_native.h"

#include <map>

namespace LAMMPS_NS {

class ReaderNativeBin : public ReaderNative {
 public:
  ReaderNativeBin(class LAMMPS *);
  ~ReaderNativeBin();

  int read_time(bigint &);
  void skip();
  bigint read_header(double[3][3], int &, int &, int, int, int *, char **, int, int, int &, int &,
                     int &, int &);
  void read_atoms(int, int, double **);

 private:
  char *line;    // line read from dump file

  int nwords;         // # of per-atom columns in dump file
  int *fieldindex;    //

  int endian;
  int revision;
  char *magic_string;
  char *unit_style;
  int size_one;
  int maxbuf = 0;
  double *buf = nullptr;

  int find_label(const std::string &label, const std::map<std::string, int> &labels);
  void read_lines(int);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Dump file is incorrectly formatted

Self-explanatory.

E: Unexpected end of dump file

A read operation from the file failed.

*/
