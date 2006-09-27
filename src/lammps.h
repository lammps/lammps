/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

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

class System;
class Universe;
class Input;
class Memory;
class Error;

class Atom;
class Update;
class Neighbor;
class Comm;
class Domain;
class Force;
class Modify;
class Group;
class Output;
class Timer;

class LAMMPS {
 public:
  static System *sys;             // simulation system
  static Universe *universe;      // universe of processors
  static Input *input;            // input script processing
  static Memory *memory;          // memory allocation functions
  static Error *error;            // error handling

  static Atom *atom;              // atom-based quantities
  static Update *update;          // integrators/minimizers
  static Neighbor *neighbor;      // neighbor lists
  static Comm *comm;              // inter-processor communication
  static Domain *domain;          // simulation box
  static Force *force;            // inter-particle forces
  static Modify *modify;          // fixes
  static Group *group;            // groups of atoms
  static Output *output;          // thermo/dump/restart
  static Timer *timer;            // CPU timing info

  static MPI_Comm world;          // communicator for my world of procs
  static FILE *infile;            // infile for my world
  static FILE *screen;            // screen output for my world
  static FILE *logfile;           // logfile for my world

  void open(int, char **, MPI_Comm);
  void close();
};

#endif
