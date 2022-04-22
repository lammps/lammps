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
   Contributing author: Richard Berger (Temple U)
------------------------------------------------------------------------- */

#include "gz_file_writer.h"
#include "fmt/format.h"
#include <cstdio>

using namespace LAMMPS_NS;

GzFileWriter::GzFileWriter() : FileWriter(), compression_level(Z_BEST_COMPRESSION), gzFp(nullptr) {}

/* ---------------------------------------------------------------------- */

GzFileWriter::~GzFileWriter()
{
  GzFileWriter::close();
}

/* ---------------------------------------------------------------------- */

void GzFileWriter::open(const std::string &path, bool append)
{
  if (isopen()) return;

  std::string mode;
  if (append) {
    mode = fmt::format("ab{}", mode, compression_level);
  } else {
    mode = fmt::format("wb{}", mode, compression_level);
  }

  gzFp = gzopen(path.c_str(), mode.c_str());

  if (gzFp == nullptr) throw FileWriterException(fmt::format("Could not open file '{}'", path));
}

/* ---------------------------------------------------------------------- */

size_t GzFileWriter::write(const void *buffer, size_t length)
{
  if (!isopen()) return 0;

  return gzwrite(gzFp, buffer, length);
}

/* ---------------------------------------------------------------------- */

void GzFileWriter::flush()
{
  if (!isopen()) return;

  gzflush(gzFp, Z_SYNC_FLUSH);
}

/* ---------------------------------------------------------------------- */

void GzFileWriter::close()
{
  if (!GzFileWriter::isopen()) return;

  gzclose(gzFp);
  gzFp = nullptr;
}

/* ---------------------------------------------------------------------- */

bool GzFileWriter::isopen() const
{
  return gzFp;
}

/* ---------------------------------------------------------------------- */

void GzFileWriter::setCompressionLevel(int level)
{
  if (isopen())
    throw FileWriterException("Compression level can not be changed while file is open");

  const int min_level = Z_DEFAULT_COMPRESSION;
  const int max_level = Z_BEST_COMPRESSION;

  if (level < min_level || level > max_level)
    throw FileWriterException(
        fmt::format("Compression level must in the range of [{}, {}]", min_level, max_level));

  compression_level = level;
}
