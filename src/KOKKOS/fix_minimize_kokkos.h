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

#ifdef FIX_CLASS

FixStyle(MINIMIZE/kk,FixMinimizeKokkos)
FixStyle(MINIMIZE/kk/device,FixMinimizeKokkos)
FixStyle(MINIMIZE/kk/host,FixMinimizeKokkos)

#else

#ifndef LMP_FIX_MINIMIZE_KOKKOS_H
#define LMP_FIX_MINIMIZE_KOKKOS_H

#include "fix_minimize.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class FixMinimizeKokkos : public FixMinimize {
  friend class MinLineSearchKokkos;

 public:
  FixMinimizeKokkos(class LAMMPS *, int, char **);
  virtual ~FixMinimizeKokkos();
  void init() {}

  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

  void add_vector_kokkos();
  DAT::t_float_1d request_vector_kokkos(int);
  void reset_coords();

  DAT::tdual_float_2d k_vectors;
  DAT::t_float_2d d_vectors;
  HAT::t_float_2d h_vectors;
};

}

#endif
#endif
/* ERROR/WARNING messages:

*/
