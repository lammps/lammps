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

#ifndef LMP_NEIGH_LIST_KOKKOS_H
#define LMP_NEIGH_LIST_KOKKOS_H

#include "pointers.h"

#include "neigh_list.h"         // IWYU pragma: export
#include "kokkos_type.h"

namespace LAMMPS_NS {

class AtomNeighbors
{
 public:
  const int num_neighs;

  KOKKOS_INLINE_FUNCTION
  AtomNeighbors(int* const & firstneigh, const int & _num_neighs,
                const int & stride):
  num_neighs(_num_neighs), _firstneigh(firstneigh), _stride(stride) {};
  KOKKOS_INLINE_FUNCTION
  int& operator()(const int &i) const {
    return _firstneigh[(bigint) i*_stride];
  }

 private:
  int* const _firstneigh;
  const int _stride;
};

class AtomNeighborsConst
{
 public:
  const int* const _firstneigh;
  const int num_neighs;

  KOKKOS_INLINE_FUNCTION
  AtomNeighborsConst(const int* const & firstneigh, const int & _num_neighs,
                     const int & stride):
  _firstneigh(firstneigh), num_neighs(_num_neighs), _stride(stride) {};
  KOKKOS_INLINE_FUNCTION
  const int& operator()(const int &i) const {
    return _firstneigh[(bigint) i*_stride];
  }

 private:
  //const int* const _firstneigh;
  const int _stride;
};

template<class DeviceType>
class NeighListKokkos: public NeighList {
  int _stride;

public:
  int maxneighs;

  void grow(int nmax);
  typename ArrayTypes<DeviceType>::t_neighbors_2d d_neighbors;
  typename ArrayTypes<DeviceType>::t_neighbors_2d_lr d_neighbors_transpose;
  DAT::tdual_int_1d k_ilist;   // local indices of I atoms
  typename ArrayTypes<DeviceType>::t_int_1d d_ilist;
  typename ArrayTypes<DeviceType>::t_int_1d d_numneigh;

  NeighListKokkos(class LAMMPS *lmp);

  KOKKOS_INLINE_FUNCTION
  AtomNeighbors get_neighbors(const int &i) const {
    return AtomNeighbors(&d_neighbors(i,0),d_numneigh(i),
                         &d_neighbors(i,1)-&d_neighbors(i,0));
  }

  KOKKOS_INLINE_FUNCTION
  AtomNeighbors get_neighbors_transpose(const int &i) const {
    return AtomNeighbors(&d_neighbors_transpose(i,0),d_numneigh(i),
                         &d_neighbors_transpose(i,1)-&d_neighbors_transpose(i,0));
  }

  KOKKOS_INLINE_FUNCTION
  static AtomNeighborsConst static_neighbors_const(int i,
           typename ArrayTypes<DeviceType>::t_neighbors_2d_const const& d_neighbors,
           typename ArrayTypes<DeviceType>::t_int_1d_const const& d_numneigh) {
    return AtomNeighborsConst(&d_neighbors(i,0),d_numneigh(i),
                              &d_neighbors(i,1)-&d_neighbors(i,0));
  }

  KOKKOS_INLINE_FUNCTION
  AtomNeighborsConst get_neighbors_const(const int &i) const {
    return AtomNeighborsConst(&d_neighbors(i,0),d_numneigh(i),
                              &d_neighbors(i,1)-&d_neighbors(i,0));
  }

  KOKKOS_INLINE_FUNCTION
  int& num_neighs(const int & i) const {
    return d_numneigh(i);
  }
 private:
  int maxatoms;
};

}

#endif
