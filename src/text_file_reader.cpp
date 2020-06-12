/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
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

#include "lammps.h"
#include "force.h"
#include "error.h"
#include "comm.h"
#include "utils.h"
#include "text_file_reader.h"
#include "tokenizer.h"
#include "fmt/format.h"

#include <cstring>

using namespace LAMMPS_NS;

TextFileReader::TextFileReader(const std::string &filename, const std::string &filetype)
  : filename(filename), filetype(filetype), ignore_comments(true)
{
  fp = fopen(filename.c_str(), "r");

  if (fp == nullptr) {
    throw FileReaderException(fmt::format("cannot open {} file {}", filetype, filename));
  }
}

TextFileReader::~TextFileReader() {
  fclose(fp);
}

void TextFileReader::skip_line() {
  char *ptr = fgets(line, MAXLINE, fp);
  if (ptr == nullptr) {
    // EOF
    throw EOFException(fmt::format("Missing line in {} file!", filetype));
  }
}

char *TextFileReader::next_line(int nparams) {
  // concatenate lines until have nparams words
  int n = 0;
  int nwords = 0;

  char *ptr = fgets(line, MAXLINE, fp);

  if (ptr == nullptr) {
    // EOF
    return nullptr;
  }

  // strip comment
  if (ignore_comments && (ptr = strchr(line, '#'))) *ptr = '\0';

  nwords = utils::count_words(line);

  if (nwords > 0) {
    n = strlen(line);
  }

  while(nwords == 0 || nwords < nparams) {
    char *ptr = fgets(&line[n], MAXLINE - n, fp);

    if (ptr == nullptr) {
      // EOF
      if (nwords > 0 && nwords < nparams) {
        throw EOFException(fmt::format("Incorrect format in {} file! {}/{} parameters", filetype, nwords, nparams));
      }
      return nullptr;
    }


    // strip comment
    if (ignore_comments && (ptr = strchr(line, '#'))) *ptr = '\0';

    nwords += utils::count_words(&line[n]);

    // skip line if blank
    if (nwords > 0) {
      n = strlen(line);
    }
  }

  return line;
}

void TextFileReader::next_dvector(double * list, int n) {
  int i = 0;
  while (i < n) {
    char *ptr = next_line();

    if (ptr == nullptr) {
      // EOF
      if (i < n) {
        throw FileReaderException(fmt::format("Incorrect format in {} file! {}/{} values", filetype, i, n));
      }
    }

    ValueTokenizer values(line);
    while(values.has_next()) {
      list[i++] = values.next_double();
    }
  }
}

ValueTokenizer TextFileReader::next_values(int nparams, const std::string & separators) {
  return ValueTokenizer(next_line(nparams), separators);
}
