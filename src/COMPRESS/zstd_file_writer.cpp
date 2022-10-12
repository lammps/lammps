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

#ifdef LAMMPS_ZSTD

#include "zstd_file_writer.h"
#include "fmt/format.h"
#include <cstdio>

using namespace LAMMPS_NS;

ZstdFileWriter::ZstdFileWriter() :
    compression_level(0), checksum_flag(1), cctx(nullptr), fp(nullptr)
{
  out_buffer_size = ZSTD_CStreamOutSize();
  out_buffer = new char[out_buffer_size];
}

/* ---------------------------------------------------------------------- */

ZstdFileWriter::~ZstdFileWriter()
{
  ZstdFileWriter::close();

  delete[] out_buffer;
  out_buffer = nullptr;
  out_buffer_size = 0;
}

/* ---------------------------------------------------------------------- */

void ZstdFileWriter::open(const std::string &path, bool append)
{
  if (isopen()) return;

  if (append) {
    fp = fopen(path.c_str(), "ab");
  } else {
    fp = fopen(path.c_str(), "wb");
  }

  if (!fp) { throw FileWriterException(fmt::format("Could not open file '{}'", path)); }

  cctx = ZSTD_createCCtx();

  if (!cctx) {
    fclose(fp);
    fp = nullptr;
    throw FileWriterException("Could not create Zstd context");
  }

  ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, compression_level);
  ZSTD_CCtx_setParameter(cctx, ZSTD_c_checksumFlag, checksum_flag);
}

/* ---------------------------------------------------------------------- */

size_t ZstdFileWriter::write(const void *buffer, size_t length)
{
  if (!isopen()) return 0;

  ZSTD_inBuffer input = {buffer, length, 0};
  ZSTD_EndDirective mode = ZSTD_e_continue;

  do {
    ZSTD_outBuffer output = {out_buffer, out_buffer_size, 0};
    ZSTD_compressStream2(cctx, &output, &input, mode);
    fwrite(out_buffer, sizeof(char), output.pos, fp);
  } while (input.pos < input.size);

  return length;
}

/* ---------------------------------------------------------------------- */

void ZstdFileWriter::flush()
{
  if (!isopen()) return;

  size_t remaining;
  ZSTD_inBuffer input = {nullptr, 0, 0};
  ZSTD_EndDirective mode = ZSTD_e_flush;

  do {
    ZSTD_outBuffer output = {out_buffer, out_buffer_size, 0};
    remaining = ZSTD_compressStream2(cctx, &output, &input, mode);
    fwrite(out_buffer, sizeof(char), output.pos, fp);
  } while (remaining);

  fflush(fp);
}

/* ---------------------------------------------------------------------- */

void ZstdFileWriter::close()
{
  if (!ZstdFileWriter::isopen()) return;

  size_t remaining;
  ZSTD_inBuffer input = {nullptr, 0, 0};
  ZSTD_EndDirective mode = ZSTD_e_end;

  do {
    ZSTD_outBuffer output = {out_buffer, out_buffer_size, 0};
    remaining = ZSTD_compressStream2(cctx, &output, &input, mode);
    fwrite(out_buffer, sizeof(char), output.pos, fp);
  } while (remaining);

  ZSTD_freeCCtx(cctx);
  cctx = nullptr;
  fclose(fp);
  fp = nullptr;
}

/* ---------------------------------------------------------------------- */

bool ZstdFileWriter::isopen() const
{
  return fp && cctx;
}

/* ---------------------------------------------------------------------- */

void ZstdFileWriter::setCompressionLevel(int level)
{
  if (isopen())
    throw FileWriterException("Compression level can not be changed while file is open");

  const int min_level = ZSTD_minCLevel();
  const int max_level = ZSTD_maxCLevel();

  if (level < min_level || level > max_level)
    throw FileWriterException(
        fmt::format("Compression level must in the range of [{}, {}]", min_level, max_level));

  compression_level = level;
}

/* ---------------------------------------------------------------------- */

void ZstdFileWriter::setChecksum(bool enabled)
{
  if (isopen()) throw FileWriterException("Checksum flag can not be changed while file is open");
  checksum_flag = enabled ? 1 : 0;
}

#endif
