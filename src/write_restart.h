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
CommandStyle(write_restart,WriteRestart);
// clang-format on
#else

#ifndef LMP_WRITE_RESTART_H
#define LMP_WRITE_RESTART_H

#include "command.h"
#include "safe_pointers.h"
#include "restartable.h"
#include "file_writer.h"

namespace LAMMPS_NS {

class WriteRestart : public Command {
 public:
  WriteRestart(class LAMMPS *);
  void command(int, char **) override;
  void multiproc_options(int, int, char **);
  void write(const std::string &);

  void write_restart_global(FileWriter&);
  void write_restart_local(FileWriter&);

 protected:
  int me, nprocs;
  SafeFilePtr fp;

  enum WriteMode {
    WRITE_RESTART,
    WRITE_RESTART_GLOBAL,
    WRITE_RESTART_LOCAL
  } mode;

  bigint natoms;    // natoms (sum of nlocal) to write into file
  int noinit;

  int multiproc;        // 0 = proc 0 writes for all
                        // else # of procs writing files
  int nclusterprocs;    // # of procs in my cluster that write to one file
  int filewriter;       // 1 if this proc writes a file, else 0
  int fileproc;         // ID of proc in my cluster who writes to file
  int icluster;         // which cluster I am in

  void global(FileWriter&);

  void header(FileWriter&);
  void type_arrays(FileWriter&);
  void force_fields(FileWriter&);
  void file_layout(FileWriter&, int, int);

  bool is_restartable(const Restartable *);
  void write_restartable(FileWriter&, Restartable *);

  void magic_string(FileWriter&);

  template<typename T>
  void write_val(FileWriter& fw, const int& tag, const T& val) {
    fw.writev(tag);
    fw.writev(val);
  }
  template<typename T>
  void write_vec(FileWriter& fw, const int& tag, const int& len, const T *vec) {
    write_val(fw, tag, len);
    fw.writev(vec, len);
  }
};
}    // namespace LAMMPS_NS
#endif
#endif
