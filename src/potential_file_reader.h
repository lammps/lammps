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
   Contributing authors: Richard Berger (Temple U)
------------------------------------------------------------------------- */

#ifndef LMP_POTENTIAL_FILE_READER_H
#define LMP_POTENTIAL_FILE_READER_H

#include "pointers.h"     // IWYU pragma: export
#include "tokenizer.h"    // IWYU pragma: export

namespace LAMMPS_NS {
class TextFileReader;

class PotentialFileReader : protected Pointers {
 protected:
  TextFileReader *reader;
  std::string filename;
  std::string filetype;
  int unit_convert;

  TextFileReader *open_potential(const std::string &path);

 public:
  PotentialFileReader(class LAMMPS *lmp, const std::string &filename,
                      const std::string &potential_name, const int auto_convert = 0);
  PotentialFileReader(class LAMMPS *lmp, const std::string &filename,
                      const std::string &potential_name, const std::string &name_suffix,
                      const int auto_convert = 0);
  ~PotentialFileReader() override;

  void ignore_comments(bool value);

  void rewind();
  void skip_line();
  char *next_line(int nparams = 0);
  void next_dvector(double *list, int n);
  ValueTokenizer next_values(int nparams,
                             const std::string &separators = TOKENIZER_DEFAULT_SEPARATORS);

  // convenience functions
  double next_double();
  int next_int();
  tagint next_tagint();
  bigint next_bigint();
  std::string next_string();

  // unit conversion info
  int get_unit_convert() const { return unit_convert; }
};

}    // namespace LAMMPS_NS

#endif
