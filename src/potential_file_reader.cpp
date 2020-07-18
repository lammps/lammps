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
#include "potential_file_reader.h"
#include "update.h"
#include "utils.h"
#include "tokenizer.h"
#include "fmt/format.h"

#include <cstring>

using namespace LAMMPS_NS;

PotentialFileReader::PotentialFileReader(LAMMPS *lmp,
                                         const std::string &filename,
                                         const std::string &potential_name,
                                         const int auto_convert) :
  Pointers(lmp),
  reader(nullptr),
  filename(filename),
  filetype(potential_name + " potential"),
  unit_convert(auto_convert)
{
  if (comm->me != 0) {
    error->one(FLERR, "FileReader should only be called by proc 0!");
  }

  try {
    reader = open_potential(filename);
    if(!reader) {
      error->one(FLERR, fmt::format("cannot open {} potential file {}", potential_name, filename));
    }
  } catch (FileReaderException & e) {
    error->one(FLERR, e.what());
  }
}

PotentialFileReader::~PotentialFileReader() {
  delete reader;
}

void PotentialFileReader::ignore_comments(bool value) {
  reader->ignore_comments = value;
}

void PotentialFileReader::skip_line() {
  try {
    reader->skip_line();
  } catch (FileReaderException & e) {
    error->one(FLERR, e.what());
  }
}

char *PotentialFileReader::next_line(int nparams) {
  try {
    return reader->next_line(nparams);
  } catch (FileReaderException & e) {
    error->one(FLERR, e.what());
  }
  return nullptr;
}

void PotentialFileReader::next_dvector(double * list, int n) {
  try {
    return reader->next_dvector(list, n);
  } catch (FileReaderException & e) {
    error->one(FLERR, e.what());
  }
}

ValueTokenizer PotentialFileReader::next_values(int nparams, const std::string & separators) {
  try {
    return reader->next_values(nparams, separators);
  } catch (FileReaderException & e) {
    error->one(FLERR, e.what());
  }
  return ValueTokenizer("");
}

double PotentialFileReader::next_double() {
  try {
    char * line = reader->next_line(1);
    return ValueTokenizer(line).next_double();
  } catch (FileReaderException & e) {
    error->one(FLERR, e.what());
  }
  return 0.0;
}

int PotentialFileReader::next_int() {
  try {
    char * line = reader->next_line(1);
    return ValueTokenizer(line).next_int();
  } catch (FileReaderException & e) {
    error->one(FLERR, e.what());
  }
  return 0;
}

tagint PotentialFileReader::next_tagint() {
  try {
    char * line = reader->next_line(1);
    return ValueTokenizer(line).next_tagint();
  } catch (FileReaderException & e) {
    error->one(FLERR, e.what());
  }
  return 0;
}

bigint PotentialFileReader::next_bigint() {
  try {
    char * line = reader->next_line(1);
    return ValueTokenizer(line).next_bigint();
  } catch (FileReaderException & e) {
    error->one(FLERR, e.what());
  }
  return 0;
}

std::string PotentialFileReader::next_string() {
  try {
    char * line = reader->next_line(1);
    return ValueTokenizer(line).next_string();
  } catch (FileReaderException & e) {
    error->one(FLERR, e.what());
  }
  return "";
}

TextFileReader *PotentialFileReader::open_potential(const std::string &path) {
  std::string filepath = utils::get_potential_file_path(path);

  if (!filepath.empty()) {
    std::string unit_style = lmp->update->unit_style;
    std::string date       = utils::get_potential_date(filepath, filetype);
    std::string units      = utils::get_potential_units(filepath, filetype);

    if (!date.empty())
      utils::logmesg(lmp, fmt::format("Reading {} file {} with DATE: {}\n",
                                      filetype, filename, date));

    if (units.empty()) {
      unit_convert = utils::NOCONVERT;
    } else {
      if (units == unit_style) {
        unit_convert = utils::NOCONVERT;
      } else {
        if ((units == "metal") && (unit_style == "real") && (unit_convert & utils::METAL2REAL)) {
          unit_convert = utils::METAL2REAL;
        } else if ((units == "real") && (unit_style == "metal") && (unit_convert & utils::REAL2METAL)) {
          unit_convert = utils::REAL2METAL;
        } else {
          lmp->error->one(FLERR, fmt::format("{} file {} requires {} units "
                                             "but {} units are in use", filetype,
                                             filename, units, unit_style));
        }
      }
    }
    if (unit_convert != utils::NOCONVERT)
      lmp->error->warning(FLERR, fmt::format("Converting {} in {} units to {} "
                                             "units", filetype, units, unit_style));
    return new TextFileReader(filepath, filetype);
  }
  return nullptr;
}
