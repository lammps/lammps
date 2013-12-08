/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_MPIIO_H
#define LMP_MPIIO_H

// true interface to MPIIO package
// used when MPIIO pacakge is installed

#ifdef LMP_MPIIO

#include "restart_mpiio.h"

#else

// dummy interface to MPIIO package
// needed for compiling when MPIIO package is not installed

namespace LAMMPS_NS {

class RestartMPIIO {
 public:
  int mpiio_exists;

  RestartMPIIO(class LAMMPS *) {mpiio_exists = 0;}
  ~RestartMPIIO() {}
  void open() {}
  void write() {}
  void read() {}
  void close() {}
};

}

#endif
#endif
