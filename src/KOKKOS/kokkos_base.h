/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/*
Copyright (C) 2020 Advanced Micro Devices, Inc. All Rights Reserved.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#ifndef KOKKOS_BASE_H
#define KOKKOS_BASE_H

#include "kokkos_type.h"

namespace LAMMPS_NS {

class KokkosBase {
 public:
  KokkosBase() {}

  // Pair
  virtual int pack_forward_comm_kokkos(int, DAT::tdual_int_2d,
                                       int, DAT::tdual_xfloat_1d &,
                                       int, int *) {return 0;};
  virtual void unpack_forward_comm_kokkos(int, int, DAT::tdual_xfloat_1d &) {}

  // Fix
  virtual int pack_forward_comm_fix_kokkos(int, DAT::tdual_int_2d,
                                           int, DAT::tdual_xfloat_1d &,
                                           int, int *) {return 0;};
  virtual void unpack_forward_comm_fix_kokkos(int, int, DAT::tdual_xfloat_1d &) {}


  // Region
  virtual void match_all_kokkos(int, DAT::tdual_int_1d) {}
};

}

#endif

/* ERROR/WARNING messages:

*/
