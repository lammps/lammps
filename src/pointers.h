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
 * the constituent class instances in the LAMMPS class.
 * Since most classes in LAMMPS are either directly or indirectly
 * derived from the Pointers class, they have access to
 * all pointers in the LAMMPS class.  These references
 * are initialized by the Pointers constructor and thus
 * all operations and data that are exported as public in those
 * base classes can be accessed and executed by all derived classes.
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
  LAMMPS *lmp;                  //!< pointer to LAMMPS class
  Memory *&memory;              //!< reference to Memory class instance in LAMMPS class
  Error *&error;                //!< reference to Error class instance in LAMMPS class
  Universe *&universe;          //!< reference to Univers class instance in LAMMPS class
  Input *&input;                //!< reference to Input class instance in LAMMPS class

  Atom *&atom;                  //!< reference to Atom class instance in LAMMPS class
  Update *&update;              //!< reference to Update class instance in LAMMPS class
  Neighbor *&neighbor;          //!< reference to Neighbor class instance in LAMMPS class
  Comm *&comm;                  //!< reference to Comm class instance in LAMMPS class
  Domain *&domain;              //!< reference to Domain class instance in LAMMPS class
  Force *&force;                //!< reference to Force class instance in LAMMPS class
  Modify *&modify;              //!< reference to Modify class instance in LAMMPS class
  Group *&group;                //!< reference to Group class instance in LAMMPS class
  Output *&output;              //!< reference to Output class instance in LAMMPS class
  Timer *&timer;                //!< reference to Timer class instance in LAMMPS class

  MPI_Comm &world;              //!< reference to MPI "world" communicator in LAMMPS class
  FILE *&infile;                //!< reference to input file pointer in LAMMPS class
  FILE *&screen;                //!< reference to screen output file pointer in LAMMPS class
  FILE *&logfile;               //!< reference to log file output file pointer in LAMMPS class

  class AtomKokkos *&atomKK;    //!< reference to AtomKokkos class instance in LAMMPS class
  class MemoryKokkos *&memoryKK; //!< reference to AemoryKokkos class instance in LAMMPS class
  class Python *&python;         //!< reference to Python class instance in LAMMPS class
};

}

#endif
