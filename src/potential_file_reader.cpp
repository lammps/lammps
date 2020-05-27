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

#include <cstring>

using namespace LAMMPS_NS;

PotentialFileReader::PotentialFileReader(LAMMPS *lmp,
                                         const std::string &filename,
                                         const std::string &potential_name) : Pointers(lmp) {
  if (comm->me != 0) {
    error->one(FLERR, "PotentialFileReader should only be called by proc 0!");
  }

  fp = force->open_potential(filename.c_str());

  if (fp == NULL) {
    char str[128];
    snprintf(str, 128, "cannot open %s potential file %s", potential_name.c_str(), filename.c_str());
    error->one(FLERR, str);
  }
}

PotentialFileReader::~PotentialFileReader() {
  fclose(fp);
}

char *PotentialFileReader::next_line(int nparams) {
  // concatenate lines until have nparams words
  int n = 0;
  int nwords = 0;

  do {
    char *ptr = fgets(&line[n], MAXLINE - n, fp);

    if (ptr == nullptr) {
      // EOF
      if (nwords > 0 && nwords < nparams) {
        char str[128];
        snprintf(str, 128, "Incorrect format in %s potential file", potential_name.c_str());
        error->one(FLERR, str);
      }

      return nullptr;
    }

    // strip comment
    if ((ptr = strchr(line, '#'))) *ptr = '\0';

    nwords = utils::count_words(line);

    // skip line if blank
    if (nwords == 0)
      continue;

    n = strlen(line) + 1;
  } while (nwords < nparams);

  return line;
}
