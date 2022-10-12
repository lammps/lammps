/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Richard Berger (Temple U)
------------------------------------------------------------------------- */

#ifndef LMP_TEXT_FILE_READER_H
#define LMP_TEXT_FILE_READER_H

#include "tokenizer.h"    // IWYU pragma: export

#include <cstdio>

namespace LAMMPS_NS {
class TextFileReader {
  std::string filetype;
  bool closefp;
  int bufsize;
  char *line;
  FILE *fp;

 public:
  bool ignore_comments;    //!< Controls whether comments are ignored

  TextFileReader(const std::string &filename, const std::string &filetype);
  TextFileReader(FILE *fp, std::string filetype);
  TextFileReader() = delete;
  TextFileReader(const TextFileReader &) = delete;
  TextFileReader(const TextFileReader &&) = delete;
  TextFileReader &operator=(const TextFileReader &) = delete;
  virtual ~TextFileReader();

  void set_bufsize(int);
  void rewind();
  void skip_line();
  char *next_line(int nparams = 0);

  void next_dvector(double *list, int n);
  ValueTokenizer next_values(int nparams,
                             const std::string &separators = TOKENIZER_DEFAULT_SEPARATORS);
};

class FileReaderException : public std::exception {
  std::string message;

 public:
  FileReaderException(const std::string &msg) : message(msg) {}

  const char *what() const noexcept override { return message.c_str(); }
};

class EOFException : public FileReaderException {
 public:
  EOFException(const std::string &msg) : FileReaderException(msg) {}
};

}    // namespace LAMMPS_NS

#endif
