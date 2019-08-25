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

#ifndef LMP_COMM_TILED_KOKKOS_H
#define LMP_COMM_TILED_KOKKOS_H

#include "comm_tiled.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class CommTiledKokkos : public CommTiled {
 public:
  CommTiledKokkos(class LAMMPS *);
  CommTiledKokkos(class LAMMPS *, class Comm *);
  virtual ~CommTiledKokkos();

  void forward_comm(int dummy = 0);    // forward comm of atom coords
  void reverse_comm();                 // reverse comm of forces
  void exchange();                     // move atoms to new procs
  void borders();                      // setup list of atoms to comm

  void forward_comm_pair(class Pair *);    // forward comm from a Pair
  void reverse_comm_pair(class Pair *);    // reverse comm from a Pair
  void forward_comm_fix(class Fix *, int size=0);
                                                   // forward comm from a Fix
  void reverse_comm_fix(class Fix *, int size=0);
                                                   // reverse comm from a Fix
  void reverse_comm_fix_variable(class Fix *);
                                     // variable size reverse comm from a Fix
  void forward_comm_compute(class Compute *);  // forward from a Compute
  void reverse_comm_compute(class Compute *);  // reverse from a Compute
  void forward_comm_dump(class Dump *);    // forward comm from a Dump
  void reverse_comm_dump(class Dump *);    // reverse comm from a Dump

  void forward_comm_array(int, double **);          // forward comm of array
  int exchange_variable(int, double *, double *&);  // exchange on neigh stencil

 private:

};

}

#endif

/* ERROR/WARNING messages:

*/
