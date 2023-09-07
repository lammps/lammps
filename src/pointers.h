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
// clang-format off

// Pointers class contains ptrs to master copy of
//   fundamental LAMMPS class ptrs stored in lammps.h
// every LAMMPS class inherits from Pointers to access lammps.h ptrs
// these variables are auto-initialized by Pointer class constructor
// *& variables are really pointers to the pointers in lammps.h
// & enables them to be accessed directly in any class, e.g. atom->x

#ifndef LMP_POINTERS_H
#define LMP_POINTERS_H

#include "lmptype.h"    // IWYU pragma: export

#include <mpi.h>        // IWYU pragma: export
#include <cstddef>      // IWYU pragme: export
#include <cstdio>       // IWYU pragma: export
#include <string>       // IWYU pragma: export
#include <vector>       // IWYU pragma: export

#include "fmt/format.h" // IWYU pragma: export
#include "lammps.h"     // IWYU pragma: export
#include "platform.h"   // IWYU pragma: export
#include "utils.h"      // IWYU pragma: export

namespace LAMMPS_NS {

// universal defines inside namespace

#define FLERR __FILE__,__LINE__

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

// enum used for KOKKOS host/device flags

enum ExecutionSpace{ Host, Device };

// global forward declarations

template <class T> class MyPoolChunk;
template <class T> class MyPage;

/** \class LAMMPS_NS::Pointers
 * \brief Base class for LAMMPS features
 *
 * The Pointers class contains references to many of the pointers
 * and members of the LAMMPS_NS::LAMMPS class. Derived classes thus
 * gain access to the constituent class instances in the LAMMPS
 * composite class and thus to the core functionality of LAMMPS.
 *
 * This kind of construct is needed, since the LAMMPS constructor
 * should only be run once per LAMMPS instance and thus classes
 * cannot be derived from LAMMPS itself. The Pointers class
 * constructor instead only initializes C++ references to component
 * pointer in the LAMMPS class. */

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
    logfile(ptr->logfile),
    atomKK(ptr->atomKK),
    memoryKK(ptr->memoryKK),
    python(ptr->python) {}
  virtual ~Pointers() noexcept(false) {}

  // remove other default members

  Pointers() = delete;
  Pointers(const Pointers &) = default;
  Pointers(Pointers &&) = delete;
  Pointers & operator=(const Pointers&) = delete;
  Pointers & operator=(Pointers&&) = delete;

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
