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

#ifndef LMP_POTENTIAL_FILE_READER_H
#define LMP_POTENTIAL_FILE_READER_H

#include <string>

#include "pointers.h"

namespace LAMMPS_NS
{
  class PotentialFileReader : protected Pointers {
    std::string potential_name;
    static const int MAXLINE = 1024;
    char line[MAXLINE];
    FILE *fp;

  public:
    PotentialFileReader(class LAMMPS *lmp, const std::string &filename, const std::string &potential_name);
    ~PotentialFileReader();

    char *next_line(int nparams);
  };

} // namespace LAMMPS_NS

#endif
