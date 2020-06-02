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
#include "utils.h"
#include "tokenizer.h"
#include "fmt/format.h"

#include <cstring>

using namespace LAMMPS_NS;

PotentialFileReader::PotentialFileReader(LAMMPS *lmp,
                                         const std::string &filename,
                                         const std::string &potential_name) :
  Pointers(lmp),
  reader(nullptr),
  filename(filename),
  filetype(potential_name + " potential")
{
  if (comm->me != 0) {
    error->one(FLERR, "FileReader should only be called by proc 0!");
  }

  try {
    reader = open_potential(filename);
  } catch (FileReaderException & e) {
    error->one(FLERR, e.what());
  }
}

PotentialFileReader::~PotentialFileReader() {
  delete reader;
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

void PotentialFileReader::next_dvector(int n, double * list) {
  try {
    return reader->next_dvector(n, list);
  } catch (FileReaderException & e) {
    error->one(FLERR, e.what());
  }
}

/* ----------------------------------------------------------------------
   open a potential file as specified by name
   if fails, search in dir specified by env variable LAMMPS_POTENTIALS
------------------------------------------------------------------------- */

TextFileReader * PotentialFileReader::open_potential(const std::string& path) {
  // attempt to open file directly
  // if successful, return filename
  std::string filepath = path;
  std::string filename = utils::path_basename(path);
  std::string date;

  if(utils::file_is_readable(filepath)) {
    date = get_potential_date(filepath);
  } else {
    // try the environment variable directory
    const char *path = getenv("LAMMPS_POTENTIALS");

    if (path != nullptr){
      std::string pot = utils::path_basename(filepath);
      filepath = utils::path_join(path, pot);

      if (utils::file_is_readable(filepath)) {
        date = get_potential_date(filepath);
      } else {
        return nullptr;
      }
    } else {
      return nullptr;
    }
  }

  if(!date.empty()) {
    utils::logmesg(lmp, fmt::format("Reading potential file {} with DATE: {}", filename, date));
  }

  return new TextFileReader(filepath, filetype);
}

/* ----------------------------------------------------------------------
   read first line of potential file
   if has DATE field, print following word
------------------------------------------------------------------------- */

std::string PotentialFileReader::get_potential_date(const std::string & path) {
  TextFileReader reader(path, filetype);
  reader.ignore_comments = false;
  char * line = nullptr;

  while (line = reader.next_line()) {
    ValueTokenizer values(line);
    while (values.has_next()) {
      std::string word = values.next_string();
      if (word == "DATE:") {
        if (values.has_next()) {
          std::string date = values.next_string();
          return date;
        }
      }
    }
  }
  return "";
}
