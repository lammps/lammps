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

#ifndef LMP_BUFFER_READER_ROOT_FILE_H
#define LMP_BUFFER_READER_ROOT_FILE_H

#include <mpi.h>

#include "buffer_reader_file.h"
#include "lmptype.h"

namespace LAMMPS_NS {

class BufferReaderRootFile : public BufferReaderFile {
 public:
  ~BufferReaderRootFile() = default;
  BufferReaderRootFile() : BufferReaderFile() {}
  BufferReaderRootFile(MPI_Comm mpi_comm, FILE* f, Error* e)
    : BufferReaderRootFile(mpi_comm, f, "", e) {}
  BufferReaderRootFile(MPI_Comm comm, FILE* f, const std::string& fn, Error* e)
    : BufferReaderFile(rank(comm) == 0 ? f : nullptr, fn, e), m_comm(comm) {
    MPI_Bcast(&length, 1, MPI_LMP_BIGINT, 0, m_comm);
  }

  void read_raw_buf(const ErrInfo& e, char* b, bigint count) override;
  void seek(const ErrInfo& e, bigint p) override;

  // Quick helper to maintain global state but not read data on non-root ranks
  void read_raw_buf_root_only(const ErrInfo& e, char* b, bigint count);

 private:
  MPI_Comm m_comm = MPI_COMM_SELF;
  int m_rank = rank(m_comm);

  int rank(MPI_Comm c) {
    int r;
    MPI_Comm_rank(c, &r);
    return r;
  }
};

}    // namespace LAMMPS_NS

#endif
