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

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(read_restart,ReadRestart);
// clang-format on
#else

#ifndef LMP_READ_RESTART_H
#define LMP_READ_RESTART_H

#include "command.h"
#include "safe_pointers.h"
#include "buffer_reader.h"
#include <memory>

namespace LAMMPS_NS {

class ReadRestart : public Command {
 public:
  ReadRestart(class LAMMPS *);
  void command(int, char **) override;

 private:
  int me, nprocs;
  SafeFilePtr fp;

  enum ReadMode {
    READ_RESTART,
    READ_RESTART_GLOBAL,
    READ_RESTART_LOCAL,
  } mode;

  int multiproc;         // 0 = restart file is a single file
                         // 1 = restart file is parallel (multiple files)
  int multiproc_file;    // # of parallel files in restart
  int nprocs_file;       // total # of procs that wrote restart file
  int revision;          // revision number of the restart file format

  int xperiodic, yperiodic, zperiodic;

  std::string file_search(const std::string &);
  void header(BufferReader&);
  void header_flag(BufferReader&, int);
  void type_arrays(BufferReader&);
  void force_fields(BufferReader&);

  void magic_string(BufferReader&);
  void endian(BufferReader&);
  void format_revision(BufferReader&);
  void check_eof_magic(BufferReader&);
  void file_layout(BufferReader&);

  char *read_string(BufferReader&);
  void read_int_vec(BufferReader&, int, int *);
  void read_double_vec(BufferReader&, int, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif
