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

#ifndef LMP_UNIVERSE_H
#define LMP_UNIVERSE_H

#include "pointers.h"

namespace LAMMPS_NS {

class Universe : protected Pointers {
 public:
  MPI_Comm uworld;    // communicator for entire universe
  int me, nprocs;     // my place in universe

  FILE *uscreen;     // universe screen output
  FILE *ulogfile;    // universe logfile

  int existflag;           // 1 if universe exists due to -partition flag
  int nworlds;             // # of worlds in universe
  int iworld;              // which world I am in
  int *procs_per_world;    // # of procs in each world
  int *root_proc;          // root proc in each world

  MPI_Comm uorig;    // original communicator passed to LAMMPS instance
  int *uni2orig;     // proc I in universe uworld is
                     // proc uni2orig[I] in original communicator

  Universe(class LAMMPS *, MPI_Comm);
  ~Universe();
  void reorder(char *, char *);
  void add_world(char *);
  int consistent();
};

}    // namespace LAMMPS_NS

#endif

/* ERROR/WARNING messages:

E: Invalid -reorder N value

Self-explanatory.

E: Nprocs not a multiple of N for -reorder

Self-explanatory.

E: Cannot open -reorder file

Self-explanatory.

E: Unexpected end of -reorder file

Self-explanatory.

E: Invalid entry in -reorder file

Self-explanatory.

E: Invalid command-line argument

One or more command-line arguments is invalid.  Check the syntax of
the command you are using to launch LAMMPS.

E: Invalid partition string '%s'

UNDOCUMENTED

*/
