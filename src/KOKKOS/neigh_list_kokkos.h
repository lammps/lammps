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

#ifndef LMP_NEIGH_LIST_KOKKOS_H
#define LMP_NEIGH_LIST_KOKKOS_H

#include "pointers.h"
#include "neigh_list.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

enum{FULL,HALFTHREAD,HALF,N2,FULLCLUSTER};

class AtomNeighbors
{
 public:
  const int num_neighs;

  KOKKOS_INLINE_FUNCTION
  AtomNeighbors(int* const & firstneigh, const int & _num_neighs, 
                const int & stride):
  _firstneigh(firstneigh), _stride(stride), num_neighs(_num_neighs) {};
  KOKKOS_INLINE_FUNCTION
  int& operator()(const int &i) const {
    return _firstneigh[i*_stride];
  }

 private:
  int* const _firstneigh;
  const int _stride;
};

class AtomNeighborsConst
{
 public:
  const int* const _firstneigh;
  const int numneigh;

  KOKKOS_INLINE_FUNCTION
  AtomNeighborsConst(int* const & firstneigh, const int & _numneigh, 
                     const int & stride):
  _firstneigh(firstneigh), _stride(stride), numneigh(_numneigh) {};
  KOKKOS_INLINE_FUNCTION
  const int& operator()(const int &i) const {
    return _firstneigh[i*_stride];
  }

 private:
  //const int* const _firstneigh;
  const int _stride;
};

template<class Device>
class NeighListKokkos: public NeighList {
  int _stride;

public:
  int maxneighs;

  void clean_copy();
  void grow(int nmax);
  typename ArrayTypes<Device>::t_neighbors_2d d_neighbors;
  typename ArrayTypes<Device>::t_int_1d d_ilist;   // local indices of I atoms
  typename ArrayTypes<Device>::t_int_1d d_numneigh; // # of J neighs for each I
  typename ArrayTypes<Device>::t_int_1d d_stencil;  // # of J neighs for each I
  typename ArrayTypes<LMPHostType>::t_int_1d h_stencil; // # of J neighs per I

  NeighListKokkos(class LAMMPS *lmp):
  NeighList(lmp) {_stride = 1; maxneighs = 16;};
  ~NeighListKokkos() {stencil = NULL; numneigh = NULL; ilist = NULL;};

  KOKKOS_INLINE_FUNCTION
  AtomNeighbors get_neighbors(const int &i) const {
    return AtomNeighbors(&d_neighbors(i,0),d_numneigh(i),
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
  void stencil_allocate(int smax, int style);
};

}

#endif
