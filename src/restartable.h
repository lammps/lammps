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

#ifndef LMP_RESTARTABLE_H
#define LMP_RESTARTABLE_H

#include <cstddef>

#include "lmptype.h"
#include "buffer_reader.h"

namespace LAMMPS_NS {

class FileWriter;

class Restartable {
 public:
  virtual void write_restart_global(FileWriter*) const = 0;
  virtual void write_restart_local(FileWriter*) const {}

  virtual void read_restart_global(BufferReader) = 0;
  virtual void read_restart_local(BufferReader) {}
};

} // namespace LAMMPS_NS

#endif
