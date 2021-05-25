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

void Reader::open_file(const char *file)
{
  if (fp != nullptr) close_file();

  if (utils::strmatch(file,"\\.gz$")) {
    compressed = 1;

#ifdef LAMMPS_GZIP
    auto gunzip = fmt::format("gzip -c -d {}",file);

#ifdef _WIN32
    fp = _popen(gunzip.c_str(),"rb");
#else
    fp = popen(gunzip.c_str(),"r");
#endif

#else
    error->one(FLERR,"Cannot open gzipped file without gzip support");
#endif
  } else {
    compressed = 0;
    fp = fopen(file,"r");
  }

  if (fp == nullptr)
    error->one(FLERR,"Cannot open file {}: {}",
                                 file, utils::getsyserror());
}

/* ----------------------------------------------------------------------
   close current file if open
   generic version for ASCII files that may be compressed
------------------------------------------------------------------------- */

void Reader::close_file()
{
  if (fp == nullptr) return;
  if (compressed) pclose(fp);
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
