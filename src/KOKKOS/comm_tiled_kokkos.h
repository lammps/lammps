// clang-format off
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

#ifndef LMP_COMM_TILED_KOKKOS_H
#define LMP_COMM_TILED_KOKKOS_H

#include "comm_tiled.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class CommTiledKokkos : public CommTiled {
 public:
  CommTiledKokkos(class LAMMPS *);
  CommTiledKokkos(class LAMMPS *, class Comm *);

  using CommTiled::forward_comm;
  using CommTiled::reverse_comm;
  void forward_comm(int dummy = 0) override;    // forward comm of atom coords
  void reverse_comm() override;                 // reverse comm of forces
  void exchange() override;                     // move atoms to new procs
  void borders() override;                      // setup list of atoms to comm

  void forward_comm(class Pair *) override;    // forward comm from a Pair
  void reverse_comm(class Pair *) override;    // reverse comm from a Pair
  void forward_comm(class Fix *, int size=0) override;
                                                   // forward comm from a Fix
  void reverse_comm(class Fix *, int size=0) override;
                                                   // reverse comm from a Fix
  void reverse_comm_variable(class Fix *) override;
                                     // variable size reverse comm from a Fix
  void forward_comm(class Compute *) override;  // forward from a Compute
  void reverse_comm(class Compute *) override;  // reverse from a Compute
  void forward_comm(class Dump *) override;    // forward comm from a Dump
  void reverse_comm(class Dump *) override;    // reverse comm from a Dump

  void forward_comm_array(int, double **) override;          // forward comm of array
};
}
#endif

