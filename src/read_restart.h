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

CommandStyle(read_restart,ReadRestart)

#else

#ifndef LMP_READ_RESTART_H
#define LMP_READ_RESTART_H

#include "pointers.h"

namespace LAMMPS_NS {

class ReadRestart : protected Pointers {
 public:
  ReadRestart(class LAMMPS *);
  void command(int, char **);

 private:
  int me,nprocs;
  FILE *fp;

  int multiproc;             // 0 = restart file is a single file
                             // 1 = restart file is parallel (multiple files)
  int multiproc_file;        // # of parallel files in restart
  int nprocs_file;           // total # of procs that wrote restart file
  int revision;              // revision number of the restart file format

  // MPI-IO values

  int mpiioflag;               // 1 for MPIIO output, else 0
  class RestartMPIIO *mpiio;   // MPIIO for restart file input
  bigint assignedChunkSize;
  MPI_Offset assignedChunkOffset,headerOffset;

  void file_search(char *, char *);
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

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot read_restart after simulation box is defined

The read_restart command cannot be used after a read_data,
read_restart, or create_box command.

E: Read restart MPI-IO input not allowed with % in filename

This is because a % signifies one file per processor and MPI-IO
creates one large file for all processors.

E: Reading from MPI-IO filename when MPIIO package is not installed

Self-explanatory.

E: Cannot open restart file %s

Self-explanatory.

E: Invalid flag in peratom section of restart file

The format of this section of the file is not correct.

E: Did not assign all restart atoms correctly

Atoms read in from the restart file were not assigned correctly to
processors.  This is likely due to some atom coordinates being outside
a non-periodic simulation box.  Normally this should not happen.  You
may wish to use the "remap" option on the read_restart command to see
if this helps.

E: Cannot open dir to search for restart file

Using a "*" in the name of the restart file will open the current
directory to search for matching file names.

E: Found no restart file matching pattern

When using a "*" in the restart file name, no matching file was found.

E: Restart file incompatible with current version

This is probably because you are trying to read a file created with a
version of LAMMPS that is too old compared to the current version.
Use your older version of LAMMPS and convert the restart file
to a data file.

E: Smallint setting in lmptype.h is not compatible

Smallint stored in restart file is not consistent with LAMMPS version
you are running.

E: Imageint setting in lmptype.h is not compatible

Format of imageint stored in restart file is not consistent with
LAMMPS version you are running.  See the settings in src/lmptype.h

E: Tagint setting in lmptype.h is not compatible

Format of tagint stored in restart file is not consistent with LAMMPS
version you are running.  See the settings in src/lmptype.h

E: Bigint setting in lmptype.h is not compatible

Format of bigint stored in restart file is not consistent with LAMMPS
version you are running.  See the settings in src/lmptype.h

E: Cannot run 2d simulation with non-periodic Z dimension

Use the boundary command to make the z dimension periodic in order to
run a 2d simulation.

W: Restart file used different # of processors

The restart file was written out by a LAMMPS simulation running on a
different number of processors.  Due to round-off, the trajectories of
your restarted simulation may diverge a little more quickly than if
you ran on the same # of processors.

W: Restart file used different 3d processor grid

The restart file was written out by a LAMMPS simulation running on a
different 3d grid of processors.  Due to round-off, the trajectories
of your restarted simulation may diverge a little more quickly than if
you ran on the same # of processors.

W: Restart file used different newton pair setting, using input script value

The input script value will override the setting in the restart file.

W: Restart file used different newton bond setting, using restart file value

The restart file value will override the setting in the input script.

W: Restart file used different boundary settings, using restart file values

Your input script cannot change these restart file settings.

E: Illegal or unset periodicity in restart

This error should not normally occur unless the restart file is invalid.

E: Invalid flag in header section of restart file

Unrecognized entry in restart file.

E: Invalid flag in type arrays section of restart file

Unrecognized entry in restart file.

E: Invalid flag in force field section of restart file

Unrecognized entry in restart file.

E: Restart file is not a multi-proc file

The file is inconsistent with the filename you specified for it.

E: Restart file is a multi-proc file

The file is inconsistent with the filename you specified for it.

E: Restart file is a MPI-IO file

The file is inconsistent with the filename you specified for it.

E: Restart file is not a MPI-IO file

The file is inconsistent with the filename you specified for it.

E: Invalid LAMMPS restart file

The file does not appear to be a LAMMPS restart file since
it doesn't contain the correct magic string at the beginning.

E: Restart file byte ordering is swapped

The file was written on a machine with different byte-ordering than
the machine you are reading it on.  Convert it to a text data file
instead, on the machine you wrote it on.

E: Restart file byte ordering is not recognized

The file does not appear to be a LAMMPS restart file since it doesn't
contain a recognized byte-orderomg flag at the beginning.

E: Illegal size string or corrupt restart

This error should not normally occur unless the restart file is invalid.

E: Illegal size integer vector read requested

This error should not normally occur unless the restart file is invalid.

E: Illegal size double vector read requested

This error should not normally occur unless the restart file is invalid.

*/
