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
   Contributing author: Matthew Whitlock (Sandia)
------------------------------------------------------------------------- */

#include "buffer_reader_root_file.h"
#include "error.h"
#include "utils.h"
#include "platform.h"
#include "fmt/format.h"
#include <cstring>


namespace LAMMPS_NS {

void BufferReaderRootFile::read_raw_buf_root_only(
  const ErrInfo& e, char* dst, bigint count
) {
  if (count < 0) error->all(e.file, e.line, "Attempt to read {} bytes", count);
  if (count > remaining()) error->all(e.file, e.line,
    "Attempt to read {} bytes with only {} remaining", count, remaining()
  );

  if (m_rank == 0) BufferReaderFile::read_raw_buf(e, dst, count);
  else position += count;
}

void BufferReaderRootFile::read_raw_buf(const ErrInfo& e, char* dst, bigint count) {
  read_raw_buf_root_only(e, dst, count);
  MPI_Bcast(dst, count, MPI_BYTE, 0, m_comm);
}

void BufferReaderRootFile::seek(const ErrInfo& e, bigint p) {
  if (p < 0 || p > this->size()) error->all(
    e.file, e.line, "Illegal seek to {} on {} sized file", p, this->size()
  );
  if (m_rank == 0) BufferReaderFile::seek(e, p);
  else position = p;
}

}
