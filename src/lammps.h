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

#ifndef LMP_LAMMPS_H
#define LMP_LAMMPS_H

#include <cstdio>
#include <mpi.h>

namespace LAMMPS_NS {

class LAMMPS {
 public:
  // ptrs to fundamental LAMMPS classes
  class Memory *memory;            // memory allocation functions
  class Error *error;              // error handling
  class Universe *universe;        // universe of processors
  class Input *input;              // input script processing
                                   // ptrs to top-level LAMMPS-specific classes
  class Atom *atom;                // atom-based quantities
  class Update *update;            // integrators/minimizers
  class Neighbor *neighbor;        // neighbor lists
  class Comm *comm;                // inter-processor communication
  class Domain *domain;            // simulation box
  class Force *force;              // inter-particle forces
  class Modify *modify;            // fixes and computes
  class Group *group;              // groups of atoms
  class Output *output;            // thermo/dump/restart
  class Timer *timer;              // CPU timing info
                                   //
  class KokkosLMP *kokkos;         // KOKKOS accelerator class
  class AtomKokkos *atomKK;        // KOKKOS version of Atom class
  class MemoryKokkos *memoryKK;    // KOKKOS version of Memory class
  class Python *python;            // Python interface
  class CiteMe *citeme;            // handle citation info

  const char *version;    // LAMMPS version string = date
  int num_ver;            // numeric version id derived from *version*
                          // that is constructed so that will be greater
                          // for newer versions in numeric or string
                          // value comparisons
                          //
  MPI_Comm world;         // MPI communicator
  FILE *infile;           // infile
  FILE *screen;           // screen output
  FILE *logfile;          // logfile
                          //
  double initclock;       // wall clock at instantiation
  int skiprunflag;        // 1 inserts timer command to skip run and minimize loops

  char *suffix, *suffix2, *suffixp;    // suffixes to add to input script style names
  int suffix_enable;                   // 1 if suffixes are enabled, 0 if disabled
  char *exename;                       // pointer to argv[0]

  char ***packargs;    // arguments for cmdline package commands
  int num_package;     // number of cmdline package commands

  MPI_Comm external_comm;    // MPI comm encompassing external programs
                             // when multiple programs launched by mpirun
                             // set by -mpicolor command line arg

  void *mdicomm;    // for use with MDI code coupling library

  const char *match_style(const char *style, const char *name);
  static const char *installed_packages[];
  static bool is_installed_pkg(const char *pkg);

  static bool has_git_info();
  static const char *git_commit();
  static const char *git_branch();
  static const char *git_descriptor();

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
  /// Default constructor. Declared private to prohibit its use
  LAMMPS(){};
  /// Copy constructor. Declared private to prohibit its use
  LAMMPS(const LAMMPS &){};
};

}    // namespace LAMMPS_NS

#endif
