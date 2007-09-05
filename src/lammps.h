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

#ifndef LAMMPS_H
#define LAMMPS_H

#include "mpi.h"
#include "stdio.h"

namespace LAMMPS_NS {

class LAMMPS {
 public:
                                 // ptrs to fundamental LAMMPS classes
  class Memory *memory;             // memory allocation functions
  class Error *error;               // error handling
  class Universe *universe;         // universe of processors
  class Input *input;               // input script processing
                                 // ptrs to top-level LAMMPS-specific classes
  class Atom *atom;                 // atom-based quantities
  class Update *update;             // integrators/minimizers
  class Neighbor *neighbor;         // neighbor lists
  class Comm *comm;                 // inter-processor communication
  class Domain *domain;             // simulation box
  class Force *force;               // inter-particle forces
  class Modify *modify;             // fixes and computes
  class Group *group;               // groups of atoms
  class Output *output;             // thermo/dump/restart
  class Timer *timer;               // CPU timing info

  MPI_Comm world;                // MPI communicator
  FILE *infile;                  // infile
  FILE *screen;                  // screen output
  FILE *logfile;                 // logfile

  LAMMPS(int, char **, MPI_Comm);
  ~LAMMPS();
  void create();
  void init();
  void destroy();
};

}

#endif
