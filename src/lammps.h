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

#ifndef LMP_LAMMPS_H
#define LMP_LAMMPS_H

#include "stdio.h"

namespace LAMMPS_NS {

class LAMMPS {
 public:
                                 // ptrs to fundamental LAMMPS classes
  class Memory *memory;          // memory allocation functions
  class Error *error;            // error handling
  class Universe *universe;      // universe of processors
  class Input *input;            // input script processing
                                 // ptrs to top-level LAMMPS-specific classes
  class Atom *atom;              // atom-based quantities
  class Update *update;          // integrators/minimizers
  class Neighbor *neighbor;      // neighbor lists
  class Comm *comm;              // inter-processor communication
  class Domain *domain;          // simulation box
  class Force *force;            // inter-particle forces
  class Modify *modify;          // fixes and computes
  class Group *group;            // groups of atoms
  class Output *output;          // thermo/dump/restart
  class Timer *timer;            // CPU timing info

  MPI_Comm world;                // MPI communicator
  FILE *infile;                  // infile
  FILE *screen;                  // screen output
  FILE *logfile;                 // logfile

  char *suffix;                  // suffix to add to input script style names
  int suffix_enable;             // 1 if suffix enabled, 0 if disabled
  class Cuda *cuda;              // CUDA accelerator class

  LAMMPS(int, char **, MPI_Comm);
  ~LAMMPS();
  void create();
  void post_create();
  void init();
  void destroy();

  void print_styles();
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid command-line argument

One or more command-line arguments is invalid.  Check the syntax of
the command you are using to launch LAMMPS.

E: Cannot use -reorder after -partition

Self-explanatory.  See doc page discussion of command-line switches.

E: Processor partitions are inconsistent

The total number of processors in all partitions must match the number
of processors LAMMPS is running on.

E: Must use -in switch with multiple partitions

A multi-partition simulation cannot read the input script from stdin.
The -in command-line option must be used to specify a file.

E: Can only use -pscreen with multiple partitions

Self-explanatory.  See doc page discussion of command-line switches.

E: Can only use -plog with multiple partitions

Self-explanatory.  See doc page discussion of command-line switches.

E: Cannot open universe screen file

For a multi-partition run, the master screen file cannot be opened.
Check that the directory you are running in allows for files to be
created.

E: Cannot open log.lammps

The default LAMMPS log file cannot be opened.  Check that the
directory you are running in allows for files to be created.

E: Cannot open universe log file

For a multi-partition run, the master log file cannot be opened.
Check that the directory you are running in allows for files to be
created.

E: Cannot open input script %s

Self-explanatory.

E: Cannot open screen file

The screen file specified as a command-line argument cannot be
opened.  Check that the directory you are running in allows for files
to be created.

E: Cannot open logfile

The LAMMPS log file named in a command-line argument cannot be opened.
Check that the path and name are correct.

E: Smallint setting in lmptype.h is invalid

It has to be the size of an integer.

E: Tagint setting in lmptype.h is invalid

Tagint must be as large or larger than smallint.

E: Bigint setting in lmptype.h is invalid

Size of bigint is less than size of tagint.

E: MPI_LMP_TAGINT and tagint in lmptype.h are not compatible

The size of the MPI datatype does not match the size of a tagint.

E: MPI_LMP_BIGINT and bigint in lmptype.h are not compatible

The size of the MPI datatype does not match the size of a bigint.

E: Small, tag, big integers are not sized correctly

See description of these 3 data types in src/lmptype.h.

E: 64-bit atom IDs are not yet supported

See description of this data type in src/lmptype.h.

E: Cannot use -cuda on without USER-CUDA installed

The USER-CUDA package must be installed via "make yes-user-cuda"
before LAMMPS is built.

*/
