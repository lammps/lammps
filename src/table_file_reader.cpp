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

#include "table_file_reader.h"

#include "text_file_reader.h"

using namespace LAMMPS_NS;

TableFileReader::TableFileReader(LAMMPS *lmp, const std::string &filename, const std::string &type,
                                 const int auto_convert) :
    PotentialFileReader(lmp, filename, type + " table", auto_convert)
{
}

char *TableFileReader::find_section_start(const std::string &keyword)
{
  char *line = nullptr;
  while ((line = reader->next_line())) {
    ValueTokenizer values(line);
    std::string word = values.next_string();

    if (word == keyword) {
      // matching keyword
      return line;
    }
  }
  return nullptr;
}
