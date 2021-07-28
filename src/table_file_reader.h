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

#ifndef LMP_TABLE_FILE_READER_H
#define LMP_TABLE_FILE_READER_H

#include "potential_file_reader.h"

namespace LAMMPS_NS {
class TableFileReader : public PotentialFileReader {
 public:
  TableFileReader(class LAMMPS *lmp, const std::string &filename, const std::string &type,
                  const int auto_convert = 0);
  virtual ~TableFileReader();

  char *find_section_start(const std::string &keyword);
};

}    // namespace LAMMPS_NS

#endif
