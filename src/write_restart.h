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

namespace LAMMPS_NS {

class WriteRestart : public Command {
 public:
  WriteRestart(class LAMMPS *);
  void command(int, char **) override;
  void multiproc_options(int, int, int, char **);
  void write(const std::string &);

 private:
  int me, nprocs;
  FILE *fp;
  bigint natoms;    // natoms (sum of nlocal) to write into file
  int noinit;

  int multiproc;        // 0 = proc 0 writes for all
                        // else # of procs writing files
  int nclusterprocs;    // # of procs in my cluster that write to one file
  int filewriter;       // 1 if this proc writes a file, else 0
  int fileproc;         // ID of proc in my cluster who writes to file
  int icluster;         // which cluster I am in

  // MPI-IO values

  int mpiioflag;                // 1 for MPIIO output, else 0
  class RestartMPIIO *mpiio;    // MPIIO for restart file output
  MPI_Offset headerOffset;

  void header();
  void type_arrays();
  void force_fields();
  void file_layout(int);

  void magic_string();
  void endian();
  void version_numeric();

  void write_int(int, int);
  void write_bigint(int, bigint);
  void write_double(int, double);
  void write_string(int, const std::string &);
  void write_int_vec(int, int, int *);
  void write_double_vec(int, int, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif
