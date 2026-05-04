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

#include <string>
#include <cstring>

using namespace LAMMPS_NS;

bool FileWriter::issizer() const { return false; }

size_t FileWriter::write_restart_global_size(const Restartable* r) {
  if (!r->restartable_global) {
    return writev<bigint>(0);
  } else {
    FileWriterSizer sizer;
    r->write_restart_global(sizer);
    return this->writev(sizer.size());
  }
}

size_t FileWriter::write_restart_local_size(const Restartable* r) {
  if (!r->restartable_local) {
    return writev<bigint>(0);
  } else {
    FileWriterSizer sizer;
    r->write_restart_local(sizer);
    return this->writev(sizer.size());
  }
}

template<> size_t FileWriter::writev<std::string>(const std::string& str) {
  size_t ret = writev((int)(str.size()+1));
  return ret + writev(str.c_str(), str.size()+1);
}

template<> size_t FileWriter::writev<char*>(char* const& str) {
  int n = strlen(str) + 1;
  size_t ret = writev(n);
  return ret + writev(str, n);
}
