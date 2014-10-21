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

// Pointers class contains ptrs to master copy of
//   fundamental LAMMPS class ptrs stored in lammps.h
// every LAMMPS class inherits from Pointers to access lammps.h ptrs
// these variables are auto-initialized by Pointer class constructor
// *& variables are really pointers to the pointers in lammps.h
// & enables them to be accessed directly in any class, e.g. atom->x

#ifndef LMP_POINTERS_H
#define LMP_POINTERS_H

#include "lmptype.h"
#include "mpi.h"
#include "lammps.h"

namespace LAMMPS_NS {

// universal defines inside namespace

#define FLERR __FILE__,__LINE__

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

class Pointers {
 public:
  Pointers(LAMMPS *ptr) :
    lmp(ptr),
    memory(ptr->memory),
    error(ptr->error),
    universe(ptr->universe),
    input(ptr->input),
    atom(ptr->atom),
    update(ptr->update),
    neighbor(ptr->neighbor),
    comm(ptr->comm),
    domain(ptr->domain),
    force(ptr->force),
    modify(ptr->modify),
    group(ptr->group),
    output(ptr->output),
    timer(ptr->timer),
    world(ptr->world),
    infile(ptr->infile),
    screen(ptr->screen),
    logfile(ptr->logfile) {}
  virtual ~Pointers() {}

 protected:
  LAMMPS *lmp;
  Memory *&memory;
  Error *&error;
  Universe *&universe;
  Input *&input;

  Atom *&atom;
  Update *&update;
  Neighbor *&neighbor;
  Comm *&comm;
  Domain *&domain;
  Force *&force;
  Modify *&modify;
  Group *&group;
  Output *&output;
  Timer *&timer;

  MPI_Comm &world;
  FILE *&infile;
  FILE *&screen;
  FILE *&logfile;
};

}

#endif
