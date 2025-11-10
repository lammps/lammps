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

#ifdef NPAIR_CLASS
// clang-format off
using NPairKokkosSkipDevice = NPairSkipKokkos<LMPDeviceType,0>;
NPairStyle(skip/kk/device,
           NPairKokkosSkipDevice,
           NP_SKIP | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_KOKKOS_DEVICE);

using NPairKokkosSkipGhostDevice = NPairSkipKokkos<LMPDeviceType,0>;
NPairStyle(skip/ghost/kk/device,
           NPairKokkosSkipGhostDevice,
           NP_SKIP | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_GHOST | NP_KOKKOS_DEVICE);

using NPairKokkosSkipHost = NPairSkipKokkos<LMPHostType,0>;
NPairStyle(skip/kk/host,
           NPairKokkosSkipHost,
           NP_SKIP | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_KOKKOS_HOST);

using NPairKokkosSkipGhostHost = NPairSkipKokkos<LMPHostType,0>;
NPairStyle(skip/ghost/kk/host,
           NPairKokkosSkipGhostHost,
           NP_SKIP | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_GHOST | NP_KOKKOS_HOST);

using NPairKokkosSkipTrimDevice = NPairSkipKokkos<LMPDeviceType,1>;
NPairStyle(skip/trim/kk/device,
           NPairKokkosSkipTrimDevice,
           NP_SKIP | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_TRIM |NP_KOKKOS_DEVICE);

using NPairKokkosSkipTrimGhostDevice = NPairSkipKokkos<LMPDeviceType,1>;
NPairStyle(skip/trim/ghost/kk/device,
           NPairKokkosSkipTrimGhostDevice,
           NP_SKIP | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_TRIM | NP_GHOST | NP_KOKKOS_DEVICE);

using NPairKokkosSkipTrimHost = NPairSkipKokkos<LMPHostType,1>;
NPairStyle(skip/trim/kk/host,
           NPairKokkosSkipTrimHost,
           NP_SKIP | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_TRIM | NP_KOKKOS_HOST);

using NPairKokkosSkipTrimGhostHost = NPairSkipKokkos<LMPHostType,1>;
NPairStyle(skip/trim/ghost/kk/host,
           NPairKokkosSkipTrimGhostHost,
           NP_SKIP | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_TRIM | NP_GHOST | NP_KOKKOS_HOST);

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

template<class DeviceType, int TRIM>
class NPairSkipKokkos : public NPair {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef int value_type;

  NPairSkipKokkos(class LAMMPS *);
  void build(class NeighList *) override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagNPairSkipCompute, const int&, int&, const bool&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagNPairSkipCountLocal, const int&, int&) const;

 private:
  int nlocal,num_skip,cutsq_custom;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
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

