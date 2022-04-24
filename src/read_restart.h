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

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(read_restart,ReadRestart);
// clang-format on
#else

#ifndef LMP_READ_RESTART_H
#define LMP_READ_RESTART_H

#include "command.h"

namespace LAMMPS_NS {

class ReadRestart : public Command {
 public:
  ReadRestart(class LAMMPS *);
  void command(int, char **) override;

 private:
  int me, nprocs;
  FILE *fp;

  int multiproc;         // 0 = restart file is a single file
                         // 1 = restart file is parallel (multiple files)
  int multiproc_file;    // # of parallel files in restart
  int nprocs_file;       // total # of procs that wrote restart file
  int revision;          // revision number of the restart file format

  // MPI-IO values

  int mpiioflag;                // 1 for MPIIO output, else 0
  class RestartMPIIO *mpiio;    // MPIIO for restart file input
  bigint assignedChunkSize;
  MPI_Offset assignedChunkOffset, headerOffset;

  std::string file_search(const std::string &);
  void header();
  void type_arrays();
  void force_fields();

  void magic_string();
  void endian();
  void format_revision();
  void check_eof_magic();
  void file_layout();

  int read_int();
  bigint read_bigint();
  double read_double();
  char *read_string();
  void read_int_vec(int, int *);
  void read_double_vec(int, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif
