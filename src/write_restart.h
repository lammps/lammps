/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(write_restart,WriteRestart)

#else

#ifndef LMP_WRITE_RESTART_H
#define LMP_WRITE_RESTART_H

#include "pointers.h"

namespace LAMMPS_NS {

class WriteRestart : protected Pointers {
 public:
  WriteRestart(class LAMMPS *);
  void command(int, char **);
  void multiproc_options(int, int, int, char **);
  void write(char *);

 private:
  int me,nprocs;
  FILE *fp;
  bigint natoms;         // natoms (sum of nlocal) to write into file
  int noinit;

  int multiproc;             // 0 = proc 0 writes for all
                             // else # of procs writing files
  int nclusterprocs;         // # of procs in my cluster that write to one file
  int filewriter;            // 1 if this proc writes a file, else 0
  int fileproc;              // ID of proc in my cluster who writes to file
  int icluster;              // which cluster I am in

  // MPI-IO values

  int mpiioflag;               // 1 for MPIIO output, else 0
  class RestartMPIIO *mpiio;   // MPIIO for restart file output
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
  void write_string(int, const char *);
  void write_int_vec(int, int, int *);
  void write_double_vec(int, int, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Write_restart command before simulation box is defined

The write_restart command cannot be used before a read_data,
read_restart, or create_box command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Restart file MPI-IO output not allowed with % in filename

This is because a % signifies one file per processor and MPI-IO
creates one large file for all processors.

E: Writing to MPI-IO filename when MPIIO package is not installed

Self-explanatory.

E: Cannot use write_restart fileper without % in restart file name

Self-explanatory.

E: Cannot use write_restart nfile without % in restart file name

Self-explanatory.

E: Atom count is inconsistent, cannot write restart file

Sum of atoms across processors does not equal initial total count.
This is probably because you have lost some atoms.

E: Cannot open restart file %s

Self-explanatory.

*/
