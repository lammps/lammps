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

#ifndef KOKKOS_BASE_H
#define KOKKOS_BASE_H

#include "kokkos_type.h"

namespace LAMMPS_NS {

class KokkosBase {
 public:
  KokkosBase() {}

  //Kspace
  virtual void pack_forward_kspace_kokkos(int, DAT::tdual_FFT_SCALAR_1d &, int, DAT::tdual_int_2d &, int) {};
  virtual void unpack_forward_kspace_kokkos(int, DAT::tdual_FFT_SCALAR_1d &, int, DAT::tdual_int_2d &, int) {};
  virtual void pack_reverse_kspace_kokkos(int, DAT::tdual_FFT_SCALAR_1d &, int, DAT::tdual_int_2d &, int) {};
  virtual void unpack_reverse_kspace_kokkos(int, DAT::tdual_FFT_SCALAR_1d &, int, DAT::tdual_int_2d &, int) {};

  // Pair
  virtual int pack_forward_comm_kokkos(int, DAT::tdual_int_2d,
                                       int, DAT::tdual_xfloat_1d &,
                                       int, int *) {return 0;};
  virtual void unpack_forward_comm_kokkos(int, int, DAT::tdual_xfloat_1d &) {}

  // Region
  virtual void match_all_kokkos(int, DAT::tdual_int_1d) {}
};

}

#endif

/* ERROR/WARNING messages:

*/
