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

#include "restartable.h"
#include "file_writer.h"
#include "file_writer_sizer.h"

using namespace LAMMPS_NS;

size_t FileWriter::write_restart_global_size(const Restartable* r) {
  FileWriterSizer sizer;
  r->write_restart_global(&sizer);
  return this->writev(sizer.size());
}

size_t FileWriter::write_restart_local_size(const Restartable* r) {
  FileWriterSizer sizer;
  r->write_restart_local(&sizer);
  return this->writev(sizer.size());
}
