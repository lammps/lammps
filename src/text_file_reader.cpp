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

#include "text_file_reader.h"

#include "fmt/format.h"
#include "tokenizer.h"
#include "utils.h"

#include <cstring>
#include <utility>

using namespace LAMMPS_NS;

/** Class for reading and parsing text files
 *
 * The value of the class member variable *ignore_comments* controls
 * whether any text following the pound sign (#) should be ignored (true)
 * or not (false). Default: true, i.e. ignore.
\verbatim embed:rst

*See also*
   :cpp:class:`TextFileReader`

\endverbatim
 *
 * \param  filename  Name of file to be read
 * \param  filetype  Description of file type for error messages */

TextFileReader::TextFileReader(const std::string &filename, const std::string &filetype) :
    filetype(filetype), closefp(true), line(nullptr), ignore_comments(true)
{
  set_bufsize(1024);
  fp = fopen(filename.c_str(), "r");

  if (fp == nullptr) {
    throw FileReaderException(
        fmt::format("cannot open {} file {}: {}", filetype, filename, utils::getsyserror()));
  }
}

/**
 * \overload
 *
\verbatim embed:rst

This function is useful in combination with :cpp:func:`utils::open_potential`.

.. note::

   The FILE pointer is not closed in the destructor, but will be advanced
   when reading from it.

\endverbatim
 *
 * \param  fp        File descriptor of the already opened file
 * \param  filetype  Description of file type for error messages */

TextFileReader::TextFileReader(FILE *fp, std::string filetype) :
    filetype(std::move(filetype)), closefp(false), line(nullptr), fp(fp), ignore_comments(true)
{
  set_bufsize(1024);
  if (fp == nullptr) throw FileReaderException("Invalid file descriptor");
}

/** Closes the file */

TextFileReader::~TextFileReader()
{
  if (closefp) fclose(fp);
  delete[] line;
}

/** adjust line buffer size */

void TextFileReader::set_bufsize(int newsize)
{
  if (newsize < 100) {
    throw FileReaderException(
        fmt::format("line buffer size {} for {} file too small, must be > 100", newsize, filetype));
  }
  delete[] line;
  bufsize = newsize;
  line = new char[bufsize];
}

/** Reset file to the beginning */

void TextFileReader::rewind()
{
  ::rewind(fp);
}

/** Read the next line and ignore it */

void TextFileReader::skip_line()
{
  char *ptr = fgets(line, bufsize, fp);
  if (ptr == nullptr) {
    // EOF
    throw EOFException(fmt::format("Missing line in {} file!", filetype));
  }
}

/** Read the next line(s) until *nparams* words have been read.
 *
 * This reads a line and counts the words in it, if the number
 * is less than the requested number, it will read the next
 * line, as well.  Output will be a string with all read lines
 * combined.  The purpose is to somewhat replicate the reading
 * behavior of formatted files in Fortran.
 *
 * If the *ignore_comments* class member has the value *true*,
 * then any text read in is truncated at the first '#' character.
 *
 * \param   nparams  Number of words that must be read. Default: 0
 * \return           String with the concatenated text */

char *TextFileReader::next_line(int nparams)
{
  // concatenate lines until have nparams words
  int n = 0;
  int nwords = 0;

  char *ptr = fgets(line, bufsize, fp);

  if (ptr == nullptr) {
    // EOF
    return nullptr;
  }

  // strip comment
  if (ignore_comments && (ptr = strchr(line, '#'))) *ptr = '\0';

  nwords = utils::count_words(line);
  if (nwords > 0) n = strlen(line);

  while (nwords == 0 || nwords < nparams) {
    ptr = fgets(&line[n], bufsize - n, fp);

    if (ptr == nullptr) {
      // EOF
      if (nwords > 0 && nwords < nparams) {
        throw EOFException(fmt::format("Incorrect format in {} file! {}/{} parameters", filetype,
                                       nwords, nparams));
      }
      return nullptr;
    }

    // strip comment
    if (ignore_comments && (ptr = strchr(line, '#'))) *ptr = '\0';

    nwords += utils::count_words(&line[n]);

    // skip line if blank
    if (nwords > 0) { n = strlen(line); }
  }

  return line;
}

/** Read lines until *n* doubles have been read and stored in array *list*
 *
 * This reads lines from the file using the next_line() function,
 * and splits them into floating-point numbers using the
 * ValueTokenizer class and stores the number is the provided list.
 *
 * \param  list  Pointer to array with suitable storage for *n* doubles
 * \param  n     Number of doubles to be read */

void TextFileReader::next_dvector(double *list, int n)
{
  int i = 0;
  while (i < n) {
    char *ptr = next_line();

    if (ptr == nullptr) {
      // EOF
      if (i < n) {
        throw FileReaderException(
            fmt::format("Incorrect format in {} file! {}/{} values", filetype, i, n));
      }
    }

    ValueTokenizer values(line);
    while (values.has_next() && i < n) { list[i++] = values.next_double(); }
  }
}

/** Read text until *nparams* words are read and passed to a tokenizer object for custom parsing.
 *
 * This reads lines from the file using the next_line() function,
 * and splits them into floating-point numbers using the
 * ValueTokenizer class and stores the number is the provided list.
 *
 * \param   nparams     Number of words to be read
 * \param   separators  String with list of separators.
 * \return              ValueTokenizer object for read in text */

ValueTokenizer TextFileReader::next_values(int nparams, const std::string &separators)
{
  char *ptr = next_line(nparams);
  if (ptr == nullptr) throw EOFException(fmt::format("Missing line in {} file!", filetype));
  return {line, separators};
}
