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


#ifndef LMP_POINTERS_H
#define LMP_POINTERS_H

#include "lmptype.h"   // IWYU pragma: export
#include <mpi.h>       // IWYU pragma: export
#include <cstddef>     // IWYU pragme: export
#include <cstdio>      // IWYU pragma: export
#include "lammps.h"    // IWYU pragma: export

namespace LAMMPS_NS {

// universal defines inside namespace

#define FLERR __FILE__,__LINE__

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

// enum used for KOKKOS host/device flags

enum ExecutionSpace{Host,Device};

// global forward declarations

template <class T> class MyPoolChunk;
template <class T> class MyPage;

/** \brief Base class for most fundamental and top-level LAMMPS classes
 *
 * The Pointers class contains references to the pointers of
 * the constituent class instances in the LAMMPS class but no
 * other functionality outside of initializing those references.
 * This way it is possible to provide access to members of
 * the LAMMPS class without having to derive from it and thus
 * avoid that the LAMMPS constructor is re-run when the derived
 * class is instantiated.
 */
class Pointers {
 public:
  /** Constructor. Get references to class members of the LAMMPS class.
   *
   * \param ptr pointer to the LAMMPS class instance
   */
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
    logfile(ptr->logfile),
    atomKK(ptr->atomKK),
    memoryKK(ptr->memoryKK),
    python(ptr->python) {}
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

  class AtomKokkos *&atomKK;
  class MemoryKokkos *&memoryKK;
  class Python *&python;
};

}

#endif
