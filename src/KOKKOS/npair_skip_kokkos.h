/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef NPAIR_CLASS
// clang-format off
typedef NPairSkipKokkos<LMPDeviceType> NPairKokkosSkipDevice;
NPairStyle(skip/kk/device,
           NPairKokkosSkipDevice,
           NP_SKIP | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_KOKKOS_DEVICE);

typedef NPairSkipKokkos<LMPDeviceType> NPairKokkosSkipGhostDevice;
NPairStyle(skip/ghost/kk/device,
           NPairKokkosSkipGhostDevice,
           NP_SKIP | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_GHOST | NP_KOKKOS_DEVICE);

typedef NPairSkipKokkos<LMPHostType> NPairKokkosSkipHost;
NPairStyle(skip/kk/host,
           NPairKokkosSkipHost,
           NP_SKIP | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_KOKKOS_HOST);

typedef NPairSkipKokkos<LMPHostType> NPairKokkosSkipGhostHost;
NPairStyle(skip/ghost/kk/host,
           NPairKokkosSkipGhostHost,
           NP_SKIP | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_GHOST | NP_KOKKOS_HOST);
// clang-format on
#else

// clang-format off
#ifndef LMP_NPAIR_SKIP_KOKKOS_H
#define LMP_NPAIR_SKIP_KOKKOS_H

#include "npair.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagNPairSkipCompute{};
struct TagNPairSkipCountLocal{};

template<class DeviceType>
class NPairSkipKokkos : public NPair {
 public:
  typedef DeviceType device_type;
  typedef int value_type;
  typedef ArrayTypes<DeviceType> AT;

  NPairSkipKokkos(class LAMMPS *);
  ~NPairSkipKokkos() {}
  void build(class NeighList *);

  KOKKOS_INLINE_FUNCTION
  void operator()(TagNPairSkipCompute, const int&, int&, const bool&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagNPairSkipCountLocal, const int&, int&) const;

 private:
  int nlocal,num_skip;

  typename AT::t_int_1d_randomread type;

  typename AT::t_int_scalar d_inum;

  typename AT::t_neighbors_2d_const d_neighbors_skip;
  typename AT::t_int_1d_const d_ilist_skip;
  typename AT::t_int_1d_const d_numneigh_skip;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d d_ilist;
  typename AT::t_int_1d d_numneigh;

  DAT::tdual_int_1d k_iskip;
  DAT::tdual_int_2d k_ijskip;
  typename AT::t_int_1d d_iskip;
  typename AT::t_int_2d d_ijskip;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
