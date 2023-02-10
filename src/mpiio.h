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

#ifndef LMP_MPIIO_H
#define LMP_MPIIO_H

// true interface to MPIIO package
// used when MPIIO package is installed

#ifdef LMP_MPIIO

#if defined(MPI_STUBS)
#error "The MPIIO package cannot be compiled in serial with MPI STUBS"
#endif

#include "restart_mpiio.h"    // IWYU pragma: export

#else

// dummy interface to MPIIO package
// needed for compiling when MPIIO package is not installed

namespace LAMMPS_NS {

class RestartMPIIO {
 public:
  int mpiio_exists;

  RestartMPIIO(class LAMMPS *) { mpiio_exists = 0; }
  ~RestartMPIIO() {}
  void openForRead(const char *) {}
  void openForWrite(const char *) {}
  void write(MPI_Offset, int, double *) {}
  void read(MPI_Offset, long, double *) {}
  void close() {}
};

}    // namespace LAMMPS_NS

#endif
#endif
