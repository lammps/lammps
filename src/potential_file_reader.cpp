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

#include "potential_file_reader.h"

#include "comm.h"
#include "error.h"
#include "text_file_reader.h"
#include "tokenizer.h"
#include "update.h"

using namespace LAMMPS_NS;

/** Class for reading and parsing LAMMPS potential files
 *
 * The value of the class member variable *ignore_comments* controls
 * whether any text following the pound sign (#) should be ignored (true)
 * or not (false). Default: true, i.e. ignore.
\verbatim embed:rst

*See also*
   :cpp:class:`TextFileReader`

\endverbatim
 *
 * \param  lmp             Pointer to LAMMPS instance
 * \param  filename        Name of file to be read
 * \param  potential_name  Name of potential style for error messages
 * \param  name_suffix     Suffix added to potential name in error messages
 * \param  auto_convert    Bitmask of supported unit conversions
 */

PotentialFileReader::PotentialFileReader(LAMMPS *lmp, const std::string &filename,
                                         const std::string &potential_name,
                                         const std::string &name_suffix, const int auto_convert) :
    Pointers(lmp),
    reader(nullptr), filename(filename), filetype(potential_name + name_suffix),
    unit_convert(auto_convert)
{
  if (comm->me != 0) { error->one(FLERR, "FileReader should only be called by proc 0!"); }

  try {
    reader = open_potential(filename);
    if (!reader) {
      error->one(FLERR, "cannot open {} potential file {}: {}", potential_name, filename,
                 utils::getsyserror());
    }
  } catch (FileReaderException &e) {
    error->one(FLERR, e.what());
  }
}

/*
 * \param  lmp             Pointer to LAMMPS instance
 * \param  filename        Name of file to be read
 * \param  potential_name  Name of potential style for error messages
 * \param  auto_convert    Bitmask of supported unit conversions
 */
PotentialFileReader::PotentialFileReader(LAMMPS *lmp, const std::string &filename,
                                         const std::string &potential_name,
                                         const int auto_convert) :
    PotentialFileReader(lmp, filename, potential_name, " potential", auto_convert)
{
}

/** Closes the file */

PotentialFileReader::~PotentialFileReader()
{
  delete reader;
}

/** Set comment (= text after '#') handling preference for the file to be read
 *
 * \param   value   Comment text is ignored if true, or not if false */
void PotentialFileReader::ignore_comments(bool value)
{
  reader->ignore_comments = value;
}

/** Read a line but ignore its content */

void PotentialFileReader::skip_line()
{
  try {
    reader->skip_line();
  } catch (FileReaderException &e) {
    error->one(FLERR, e.what());
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
 * \param   nparams  Number of words that must be read. Default: 0
 * \return           String with the concatenated text */

char *PotentialFileReader::next_line(int nparams)
{
  try {
    return reader->next_line(nparams);
  } catch (FileReaderException &e) {
    error->one(FLERR, e.what());
  }
  return nullptr;
}

/** Read lines until *n* doubles have been read and stored in array *list*
 *
 * This reads lines from the file using the next_line() function,
 * and splits them into floating-point numbers using the
 * ValueTokenizer class and stores the number is the provided list.
 *
 * \param  list  Pointer to array with suitable storage for *n* doubles
 * \param  n     Number of doubles to be read */

void PotentialFileReader::next_dvector(double *list, int n)
{
  try {
    return reader->next_dvector(list, n);
  } catch (FileReaderException &e) {
    error->one(FLERR, e.what());
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

ValueTokenizer PotentialFileReader::next_values(int nparams, const std::string &separators)
{
  try {
    return reader->next_values(nparams, separators);
  } catch (FileReaderException &e) {
    error->one(FLERR, e.what());
  }
  return ValueTokenizer("");
}

/** Read next line and convert first word to a double
 *
 * \return  Value of first word in line as double */

double PotentialFileReader::next_double()
{
  try {
    char *line = reader->next_line(1);
    return ValueTokenizer(line).next_double();
  } catch (FileReaderException &e) {
    error->one(FLERR, e.what());
  }
  return 0.0;
}

/** Read next line and convert first word to an int
 *
 * \return  Value of first word in line as int */

int PotentialFileReader::next_int()
{
  try {
    char *line = reader->next_line(1);
    return ValueTokenizer(line).next_int();
  } catch (FileReaderException &e) {
    error->one(FLERR, e.what());
  }
  return 0;
}

/** Read next line and convert first word to a tagint
 *
 * \return  Value of first word in line as tagint */

tagint PotentialFileReader::next_tagint()
{
  try {
    char *line = reader->next_line(1);
    return ValueTokenizer(line).next_tagint();
  } catch (FileReaderException &e) {
    error->one(FLERR, e.what());
  }
  return 0;
}

/** Read next line and convert first word to a bigint
 *
 * \return  Value of first word in line as bigint */

bigint PotentialFileReader::next_bigint()
{
  try {
    char *line = reader->next_line(1);
    return ValueTokenizer(line).next_bigint();
  } catch (FileReaderException &e) {
    error->one(FLERR, e.what());
  }
  return 0;
}

/** Read next line and return first word
 *
 * \return  First word of read in line */

std::string PotentialFileReader::next_string()
{
  try {
    char *line = reader->next_line(1);
    return ValueTokenizer(line).next_string();
  } catch (FileReaderException &e) {
    error->one(FLERR, e.what());
  }
  return "";
}

/** Look up and open the potential file
 *
\verbatim embed:rst

*See also*
   :cpp:func:`utils::open_potential`,
   :cpp:class:`TextFileReader`

\endverbatim
 * \param   path  Path of the potential file to open
 * \return        Pointer to TextFileReader object created */

TextFileReader *PotentialFileReader::open_potential(const std::string &path)
{
  std::string filepath = utils::get_potential_file_path(path);

  if (!filepath.empty()) {
    std::string unit_style = lmp->update->unit_style;
    std::string date = utils::get_potential_date(filepath, filetype);
    std::string units = utils::get_potential_units(filepath, filetype);

    if (!date.empty())
      utils::logmesg(lmp, "Reading {} file {} with DATE: {}\n", filetype, filename, date);

    if (units.empty()) {
      unit_convert = utils::NOCONVERT;
    } else {
      if (units == unit_style) {
        unit_convert = utils::NOCONVERT;
      } else {
        if ((units == "metal") && (unit_style == "real") && (unit_convert & utils::METAL2REAL)) {
          unit_convert = utils::METAL2REAL;
        } else if ((units == "real") && (unit_style == "metal") &&
                   (unit_convert & utils::REAL2METAL)) {
          unit_convert = utils::REAL2METAL;
        } else {
          lmp->error->one(FLERR, "{} file {} requires {} units but {} units are in use", filetype,
                          filename, units, unit_style);
        }
      }
    }
    if (unit_convert != utils::NOCONVERT)
      lmp->error->warning(FLERR, "Converting {} in {} units to {} units", filetype, units,
                          unit_style);
    return new TextFileReader(filepath, filetype);
  }
  return nullptr;
}
