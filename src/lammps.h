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

#ifndef LMP_LAMMPS_H
#define LMP_LAMMPS_H


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

  double initclock;              // wall clock at instantiation

  char *suffix,*suffix2;         // suffixes to add to input script style names
  int suffix_enable;             // 1 if suffixes are enabled, 0 if disabled
  char *exename;                 // pointer to argv[0]
  char ***packargs;              // arguments for cmdline package commands
  int num_package;               // number of cmdline package commands
  int cite_enable;               // 1 if generating log.cite, 0 if disabled

  int clientserver;              // 0 = neither, 1 = client, 2 = server
  void *cslib;                   // client/server messaging via CSlib
  MPI_Comm cscomm;               // MPI comm for client+server in mpi/one mode

  class KokkosLMP *kokkos;       // KOKKOS accelerator class
  class AtomKokkos *atomKK;      // KOKKOS version of Atom class
  class MemoryKokkos *memoryKK;  // KOKKOS version of Memory class

  class Python * python;         // Python interface

  class CiteMe *citeme;          // citation info

  const char *match_style(const char *style, const char *name);
  static const char * installed_packages[];
  static bool is_installed_pkg(const char *pkg);

  static const bool has_git_info;
  static const char git_commit[];
  static const char git_branch[];
  static const char git_descriptor[];

  LAMMPS(int, char **, MPI_Comm);
  ~LAMMPS();
  void create();
  void post_create();
  void init();
  void destroy();
  void print_config(FILE *);    // print compile time settings

 private:
  struct package_styles_lists *pkg_lists;
  void init_pkg_lists();
  void help();
  LAMMPS() {};                   // prohibit using the default constructor
  LAMMPS(const LAMMPS &) {};     // prohibit using the copy constructor
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid command-line argument

One or more command-line arguments is invalid.  Check the syntax of
the command you are using to launch LAMMPS.

E: Cannot use -reorder after -partition

Self-explanatory.  See doc page discussion of command-line switches.

E: Processor partitions do not match number of allocated processors

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

E: Cannot open log.lammps for writing

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

E: Imageint setting in lmptype.h is invalid

Imageint must be as large or larger than smallint.

E: Tagint setting in lmptype.h is invalid

Tagint must be as large or larger than smallint.

E: Bigint setting in lmptype.h is invalid

Size of bigint is less than size of tagint.

E: MPI_LMP_TAGINT and tagint in lmptype.h are not compatible

The size of the MPI datatype does not match the size of a tagint.

E: MPI_LMP_BIGINT and bigint in lmptype.h are not compatible

The size of the MPI datatype does not match the size of a bigint.

E: Small to big integers are not sized correctly

This error occurs when the sizes of smallint, imageint, tagint, bigint,
as defined in src/lmptype.h are not what is expected.  Contact
the developers if this occurs.

E: Cannot use -kokkos on without KOKKOS installed

Self-explanatory.

E: Using suffix gpu without GPU package installed

Self-explanatory.

E: Using suffix intel without USER-INTEL package installed

Self-explanatory.

E: Using suffix kk without KOKKOS package enabled

Self-explanatory.

E: Using suffix omp without USER-OMP package installed

Self-explanatory.

E: Too many -pk arguments in command line

The string formed by concatenating the arguments is too long.  Use a
package command in the input script instead.

U: Cannot use -cuda on and -kokkos on together

This is not allowed since both packages can use GPUs.

*/
