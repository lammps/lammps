// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "reader.h"

#include "error.h"

#include <cstring>

using namespace LAMMPS_NS;

// only proc 0 calls methods of this class, except for constructor/destructor

/* ---------------------------------------------------------------------- */

Reader::Reader(LAMMPS *lmp) : Pointers(lmp)
{
  fp = nullptr;
}

/* ----------------------------------------------------------------------
   try to open given file
   generic version for ASCII files that may be compressed
------------------------------------------------------------------------- */

void Reader::open_file(const std::string &file)
{
  if (fp != nullptr) close_file();

  if (platform::has_compress_extension(file)) {
    compressed = 1;
    fp = platform::compressed_read(file);
    if (!fp) error->one(FLERR,"Cannot open compressed file for reading");
  } else {
    compressed = 0;

    // leleslx: if there is .bin ending
    std::size_t dot = file.find_last_of('.');
    char reading_mode[3] = "r"; 
    if (dot != std::string::npos) {
      const std::string ext = file.substr(dot + 1);
      if ("bin" == ext){
        sprintf(reading_mode, "rb");
      }
    }

    fp = fopen(file.c_str(), reading_mode);
  }

  if (!fp) error->one(FLERR,"Cannot open file {}: {}", file, utils::getsyserror());
}

/* ----------------------------------------------------------------------
   close current file if open
   generic version for ASCII files that may be compressed
------------------------------------------------------------------------- */

void Reader::close_file()
{
  if (fp == nullptr) return;
  if (compressed) platform::pclose(fp);
  else fclose(fp);
  fp = nullptr;
}

/* ----------------------------------------------------------------------
   detect unused arguments
------------------------------------------------------------------------- */

void Reader::settings(int narg, char** /*args*/)
{
  if (narg > 0)
    error->all(FLERR,"Illegal read_dump command");
}
