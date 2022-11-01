/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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
  binary = false;
  compressed = false;
}

// avoid resource leak
Reader::~Reader()
{
  if (fp != nullptr) close_file();
}

/* ----------------------------------------------------------------------
   try to open given file
   generic version for ASCII files with optional compression or for native binary dumps
------------------------------------------------------------------------- */

void Reader::open_file(const std::string &file)
{
  if (fp != nullptr) close_file();

  if (platform::has_compress_extension(file)) {
    compressed = true;
    fp = platform::compressed_read(file);
    if (!fp) error->one(FLERR, "Cannot open compressed file for reading");
  } else {
    compressed = false;
    if (utils::strmatch(file, "\\.bin$")) {
      binary = true;
      fp = fopen(file.c_str(), "rb");
    } else {
      fp = fopen(file.c_str(), "r");
      binary = false;
    }
  }

  if (!fp) error->one(FLERR, "Cannot open file {}: {}", file, utils::getsyserror());
}

/* ----------------------------------------------------------------------
   close current file if open
   generic version for ASCII files that may be compressed
------------------------------------------------------------------------- */

void Reader::close_file()
{
  if (fp == nullptr) return;
  if (compressed)
    platform::pclose(fp);
  else
    fclose(fp);
  fp = nullptr;
}

/* ----------------------------------------------------------------------
   detect unused arguments
------------------------------------------------------------------------- */

void Reader::settings(int narg, char ** /*args*/)
{
  if (narg > 0) error->all(FLERR, "Illegal read_dump command");
}
