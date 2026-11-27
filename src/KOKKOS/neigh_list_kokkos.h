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

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  AtomNeighbors(int* const & firstneigh, const int & _num_neighs,
                const int & stride):
  num_neighs(_num_neighs), _firstneigh(firstneigh), _stride(stride) {};
// NOLINTNEXTLINE
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

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  AtomNeighborsConst(const int* const & firstneigh, const int & _num_neighs,
                     const int & stride):
  _firstneigh(firstneigh), num_neighs(_num_neighs), _stride(stride) {};
// NOLINTNEXTLINE
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
  typedef ArrayTypes<DeviceType> AT;

public:
  int maxneighs;

  void grow(int nmax);
  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_neighbors_2d_lr d_neighbors_transpose;
  DAT::tdual_int_1d k_ilist;   // local indices of I atoms
  typename AT::t_int_1d d_ilist;
  typename AT::t_int_1d d_numneigh;

  // Dynamic pruning (enabled by 'package kokkos neigh/prune yes'): a tight
  // inner list re-pruned from d_neighbors/d_numneigh (the master) between
  // reneighbor steps. Same (nmax, maxneighs) shape; d_ilist is shared.
  typename AT::t_neighbors_2d d_inner_neighbors;
  typename AT::t_int_1d d_inner_numneigh;
  bigint prune_last_step = -1;    // step the inner list was last pruned
  bigint prune_master_step = -1;  // master lastcall the inner list mirrors

  // Cluster-pair neighbor list (enabled by 'package kokkos neigh/cluster yes')
  int max_jclusters;
  int cluster_hash_sh = 0;         // shared-mem hash slots for the cluster build (power of 2)
  bigint cluster_built_step = -1;  // timestep of the flat-list build the cluster list mirrors
  int cluster_newton_built = -1;   // mode of the last cluster build (1 = cluster-level Newton
                                   // over a full-style flat list, 0 = legacy per-side tiles)
  typename AT::t_int_1d d_cluster_numneigh;
  typename AT::t_int_2d d_cluster_jlist;
  // d_cluster_excl(ci, 2*cj_idx)   = bits 0..31: sbmask(2 bits) for pairs  0..15
  // d_cluster_excl(ci, 2*cj_idx+1) = bits 0..31: sbmask(2 bits) for pairs 16..31
  typename AT::t_int_2d d_cluster_excl;
  // d_cluster_pres(ci, cj_idx): bit p set iff tile pair p = ki*CJ+kj is an entry
  // of the flat list. Encodes half-list single storage, deleted special pairs,
  // and ghost ownership exactly; pairs without their bit are never computed.
  typename AT::t_int_2d d_cluster_pres;
  typename AT::t_int_1d d_cluster_scratch;  // [0]=overflow flag (1=jlist,2=hash), [1]=new max j-clusters, [2]=hash-full sentinel
  // Dedicated j-side cluster ordering over all atoms including ghosts (which
  // comm leaves spatially unsorted). Built per reneighbor; d_xcl/d_typecl are
  // regathered every step so the force kernels read packed, coalesced tiles.
  int cluster_nall = 0;
  typename AT::t_int_1d d_cl2atom;        // cluster slot -> atom index
  typename AT::t_int_1d d_atom2cl;        // atom index -> cluster slot
  typename AT::t_kkfloat_1d_3_lr d_xcl;   // coords in cluster order
  typename AT::t_int_1d d_typecl;         // types in cluster order
  void grow_cluster_order(int nall);
  void grow_clusters(int num_iclusters, int max_jc);
  [[noreturn]] void cluster_fatal(const std::string &file, int line, const std::string &msg);

  NeighListKokkos(class LAMMPS *lmp);

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  AtomNeighbors get_neighbors(const int &i) const {
    return AtomNeighbors(&d_neighbors(i,0),d_numneigh(i),
                         &d_neighbors(i,1)-&d_neighbors(i,0));
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  AtomNeighbors get_neighbors_transpose(const int &i) const {
    return AtomNeighbors(&d_neighbors_transpose(i,0),d_numneigh(i),
                         &d_neighbors_transpose(i,1)-&d_neighbors_transpose(i,0));
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  static AtomNeighborsConst static_neighbors_const(int i,
           typename AT::t_neighbors_2d_const const& d_neighbors,
           typename AT::t_int_1d_const const& d_numneigh) {
    return AtomNeighborsConst(&d_neighbors(i,0),d_numneigh(i),
                              &d_neighbors(i,1)-&d_neighbors(i,0));
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  AtomNeighborsConst get_neighbors_const(const int &i) const {
    return AtomNeighborsConst(&d_neighbors(i,0),d_numneigh(i),
                              &d_neighbors(i,1)-&d_neighbors(i,0));
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int& num_neighs(const int & i) const {
    return d_numneigh(i);
  }
 private:
  int maxatoms;
};

}

#endif
