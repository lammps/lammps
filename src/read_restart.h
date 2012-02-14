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

#ifdef COMMAND_CLASS

CommandStyle(read_restart,ReadRestart)

#else

#ifndef LMP_READ_RESTART_H
#define LMP_READ_RESTART_H

#include "stdio.h"
#include "pointers.h"

namespace LAMMPS_NS {

class ReadRestart : protected Pointers {
 public:
  ReadRestart(class LAMMPS *);
  void command(int, char **);

 private:
  int me,nprocs,nprocs_file;
  FILE *fp;
  int nfix_restart_global,nfix_restart_peratom;
  int swapflag;

  void file_search(char *, char *);
  void header();
  void type_arrays();
  void force_fields();

  void nread_int(int *, int, FILE *);
  void nread_double(double *, int, FILE *);
  void nread_char(char *, int, FILE *);
  int read_int();
  double read_double();
  char *read_char();
  bigint read_bigint();
  int autodetect(FILE **, char *);
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

E: Cannot open restart file %s

Self-explanatory.

E: Did not assign all atoms correctly

Atoms read in from a data file were not assigned correctly to
processors.  This is likely due to some atom coordinates being
outside a non-periodic simulation box.

E: Cannot open dir to search for restart file

Using a "*" in the name of the restart file will open the current
directory to search for matching file names.

E: Found no restart file matching pattern

When using a "*" in the restart file name, no matching file was found.

W: Restart file version does not match LAMMPS version

This may cause problems when reading the restart file.

E: Smallint setting in lmptype.h is not compatible

Smallint stored in restart file is not consistent with LAMMPS version
you are running.

E: Tagint setting in lmptype.h is not compatible

Smallint stored in restart file is not consistent with LAMMPS version
you are running.

E: Bigint setting in lmptype.h is not compatible

Bigint stored in restart file is not consistent with LAMMPS version
you are running.

E: Cannot run 2d simulation with nonperiodic Z dimension

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

E: Invalid flag in header section of restart file

Unrecognized entry in restart file.

E: Invalid flag in type arrays section of restart file

Unrecognized entry in restart file.

E: Invalid flag in force field section of restart file

Unrecognized entry in restart file.

*/
