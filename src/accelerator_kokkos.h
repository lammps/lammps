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

#ifndef LMP_ACCELERATOR_KOKKOS_H
#define LMP_ACCELERATOR_KOKKOS_H

// true interface to KOKKOS
// used when KOKKOS is installed

#ifdef LMP_KOKKOS

#include "atom_kokkos.h"          // IWYU pragma: export
#include "comm_kokkos.h"          // IWYU pragma: export
#include "comm_tiled_kokkos.h"    // IWYU pragma: export
#include "domain_kokkos.h"        // IWYU pragma: export
#include "kokkos.h"               // IWYU pragma: export
#include "memory_kokkos.h"        // IWYU pragma: export
#include "modify_kokkos.h"        // IWYU pragma: export
#include "neighbor_kokkos.h"      // IWYU pragma: export

#define LAMMPS_INLINE KOKKOS_INLINE_FUNCTION

#else

// dummy interface to KOKKOS
// needed for compiling when KOKKOS is not installed

#include "atom.h"
#include "comm_brick.h"
#include "comm_tiled.h"
#include "domain.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"

#define LAMMPS_INLINE inline

namespace LAMMPS_NS {

class KokkosLMP {
 public:
  int kokkos_exists;
  int nthreads;
  int ngpus;

  KokkosLMP(class LAMMPS *, int, char **) { kokkos_exists = 0; }
  ~KokkosLMP() {}
  static void finalize() {}
  void accelerator(int, char **) {}
  int neigh_list_kokkos(int) { return 0; }
  int neigh_count(int) { return 0; }
};

class AtomKokkos : public Atom {
 public:
  tagint **k_special;
  AtomKokkos(class LAMMPS *lmp) : Atom(lmp) {}
  void sync(const ExecutionSpace /*space*/, unsigned int /*mask*/) {}
  void modified(const ExecutionSpace /*space*/, unsigned int /*mask*/) {}
};

class CommKokkos : public CommBrick {
 public:
  CommKokkos(class LAMMPS *lmp) : CommBrick(lmp) {}
};

class CommTiledKokkos : public CommTiled {
 public:
  CommTiledKokkos(class LAMMPS *lmp) : CommTiled(lmp) {}
  CommTiledKokkos(class LAMMPS *lmp, Comm *oldcomm) : CommTiled(lmp, oldcomm) {}
};

class DomainKokkos : public Domain {
 public:
  DomainKokkos(class LAMMPS *lmp) : Domain(lmp) {}
};

class NeighborKokkos : public Neighbor {
 public:
  NeighborKokkos(class LAMMPS *lmp) : Neighbor(lmp) {}
};

class MemoryKokkos : public Memory {
 public:
  MemoryKokkos(class LAMMPS *lmp) : Memory(lmp) {}
  void grow_kokkos(tagint **, tagint **, int, int, const char *) {}
};

class ModifyKokkos : public Modify {
 public:
  ModifyKokkos(class LAMMPS *lmp) : Modify(lmp) {}
};

class DAT {
 public:
  typedef double tdual_xfloat_1d;
  typedef double tdual_FFT_SCALAR_1d;
  typedef int tdual_int_1d;
  typedef int tdual_int_2d;
};

}    // namespace LAMMPS_NS

#endif
#endif
