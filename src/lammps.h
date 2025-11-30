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

#ifndef LMP_LAMMPS_H
#define LMP_LAMMPS_H

#include <cstdio>
#include <mpi.h>
#include <string>
#include <vector>

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
  class CiteMe *citeme;            ///< handle citation info

  const char *version;       // LAMMPS version string = date
  int num_ver;               // numeric version id derived from *version*
                             // that is constructed so that will be greater
                             // for newer versions in numeric or string
                             // value comparisons
  int restart_ver;           // -1 or numeric version id of LAMMPS version in restart
                             // file, in case LAMMPS was initialized from a restart
                             //
  MPI_Comm world;            // MPI communicator
  FILE *infile;              // file pointer of input
  FILE *screen;              // screen output
  FILE *logfile;             // logfile output
                             //
  double initclock;          // wall clock at instantiation
  int skiprunflag;           // 1 inserts timer command to skip run and minimize loops
                             //
  char *suffix, *suffix2;    // suffixes to add to input script style names
  int suffix_enable;         // 1 if suffixes are enabled, 0 if disabled
  int pair_only_flag;        // 1 if only force field pair styles are accelerated, 0 if all
  char *exename;             // pointer to argv[0]
                             //
  char ***packargs;          // arguments for cmdline package commands
  int num_package;           // number of cmdline package commands
                             //
  MPI_Comm external_comm;    // MPI comm encompassing external programs
                             // when multiple programs launched by mpirun
                             // set by -mpicolor command-line arg
                             //
  void *mdicomm;             // for use with MDI code coupling library

  // disable most auto-generated default members
  LAMMPS() = delete;
  LAMMPS(const LAMMPS &) = delete;
  LAMMPS(LAMMPS &&) = delete;
  LAMMPS &operator=(const LAMMPS &) = delete;
  LAMMPS &operator=(LAMMPS &&) = delete;

  // constructor using default C-language style argument vector
  LAMMPS(int, char **, MPI_Comm);

  // overloaded constructor with argument C++ string list argument vector
  using argv = std::vector<std::string>;
  LAMMPS(argv &args, MPI_Comm);

  ~LAMMPS() noexcept(false);

  void create();
  void post_create();
  void destroy();
  void init();

  // get name of package that a style is part of
  const char *get_style_pkg(const char *style, const char *name) const;
  static const char *installed_packages[];
  const char *non_pair_suffix() const;

  static bool has_git_info();
  static const char *git_commit();
  static const char *git_branch();
  static const char *git_descriptor();

  /// Print out LAMMPS compile time settings
  void print_config(FILE *) const;

 private:
  struct package_styles_lists *pkg_lists;
  void init_pkg_lists();
  /// Print detailed help message. Used with '-h' or '-help' command line flag.
  void help() const;
};

}    // namespace LAMMPS_NS

#endif
