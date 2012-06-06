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

CommandStyle(write_restart,WriteRestart)

#else

#ifndef LMP_WRITE_RESTART_H
#define LMP_WRITE_RESTART_H

#include "stdio.h"
#include "pointers.h"

namespace LAMMPS_NS {

class WriteRestart : protected Pointers {
 public:
  WriteRestart(class LAMMPS *);
  void command(int, char **);
  void write(char *);

 private:
  int me,nprocs;
  FILE *fp;
  bigint natoms;         // natoms (sum of nlocal) to write into file

  void header();
  void type_arrays();
  void force_fields();

  void write_int(int, int);
  void write_double(int, double);
  void write_char(int, char *);
  void write_bigint(int, bigint);
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

E: Atom count is inconsistent, cannot write restart file

Sum of atoms across processors does not equal initial total count.
This is probably because you have lost some atoms.

E: Cannot open restart file %s

Self-explanatory.

*/
