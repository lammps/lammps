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

// Trim off

// Newton, no triclinic

typedef NPairHalffullKokkos<LMPDeviceType,1,0,0> NPairKokkosHalffullNewtonDevice;
NPairStyle(halffull/newton/kk/device,
           NPairKokkosHalffullNewtonDevice,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,1,0,0> NPairKokkosHalffullNewtonHost;
NPairStyle(halffull/newton/kk/host,
           NPairKokkosHalffullNewtonHost,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_KOKKOS_HOST);

typedef NPairHalffullKokkos<LMPDeviceType,1,0,0> NPairKokkosHalffullNewtonDevice;
NPairStyle(halffull/newton/skip/kk/device,
           NPairKokkosHalffullNewtonDevice,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_SKIP | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,1,0,0> NPairKokkosHalffullNewtonHost;
NPairStyle(halffull/newton/skip/kk/host,
           NPairKokkosHalffullNewtonHost,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_SKIP | NP_KOKKOS_HOST);

// Newton, triclinic

typedef NPairHalffullKokkos<LMPDeviceType,1,1,0> NPairKokkosHalffullNewtonTriDevice;
NPairStyle(halffull/newton/tri/kk/device,
           NPairKokkosHalffullNewtonTriDevice,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,1,1,0> NPairKokkosHalffullNewtonTriHost;
NPairStyle(halffull/newton/tri/kk/host,
           NPairKokkosHalffullNewtonTriHost,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_KOKKOS_HOST);

typedef NPairHalffullKokkos<LMPDeviceType,1,1,0> NPairKokkosHalffullNewtonTriDevice;
NPairStyle(halffull/newton/tri/skip/kk/device,
           NPairKokkosHalffullNewtonTriDevice,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,1,1,0> NPairKokkosHalffullNewtonTriHost;
NPairStyle(halffull/newton/tri/skip/kk/host,
           NPairKokkosHalffullNewtonTriHost,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_KOKKOS_HOST);

// Newtoff (can be triclinic but template param always set to 0)

typedef NPairHalffullKokkos<LMPDeviceType,0,0,0> NPairKokkosHalffullNewtoffDevice;
NPairStyle(halffull/newtoff/kk/device,
           NPairKokkosHalffullNewtoffDevice,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,0,0,0> NPairKokkosHalffullNewtoffHost;
NPairStyle(halffull/newtoff/kk/host,
           NPairKokkosHalffullNewtoffHost,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_KOKKOS_HOST);

typedef NPairHalffullKokkos<LMPDeviceType,0,0,0> NPairKokkosHalffullNewtoffDevice;
NPairStyle(halffull/newtoff/skip/kk/device,
           NPairKokkosHalffullNewtoffDevice,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,0,0,0> NPairKokkosHalffullNewtoffHost;
NPairStyle(halffull/newtoff/skip/kk/host,
           NPairKokkosHalffullNewtoffHost,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_KOKKOS_HOST);

//************ Ghost **************

// Newton, no triclinic

typedef NPairHalffullKokkos<LMPDeviceType,1,0,0> NPairKokkosHalffullNewtonDevice;
NPairStyle(halffull/newton/ghost/kk/device,
           NPairKokkosHalffullNewtonDevice,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_GHOST | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,1,0,0> NPairKokkosHalffullNewtonHost;
NPairStyle(halffull/newton/ghost/kk/host,
           NPairKokkosHalffullNewtonHost,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_GHOST | NP_KOKKOS_HOST);

typedef NPairHalffullKokkos<LMPDeviceType,1,0,0> NPairKokkosHalffullNewtonDevice;
NPairStyle(halffull/newton/skip/ghost/kk/device,
           NPairKokkosHalffullNewtonDevice,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_GHOST | NP_SKIP | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,1,0,0> NPairKokkosHalffullNewtonHost;
NPairStyle(halffull/newton/skip/ghost/kk/host,
           NPairKokkosHalffullNewtonHost,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_GHOST | NP_SKIP | NP_KOKKOS_HOST);

// Newton, triclinic

typedef NPairHalffullKokkos<LMPDeviceType,1,1,0> NPairKokkosHalffullNewtonTriDevice;
NPairStyle(halffull/newton/tri/ghost/kk/device,
           NPairKokkosHalffullNewtonTriDevice,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,1,1,0> NPairKokkosHalffullNewtonTriHost;
NPairStyle(halffull/newton/tri/ghost/kk/host,
           NPairKokkosHalffullNewtonTriHost,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_KOKKOS_HOST);

typedef NPairHalffullKokkos<LMPDeviceType,1,1,0> NPairKokkosHalffullNewtonTriDevice;
NPairStyle(halffull/newton/tri/skip/ghost/kk/device,
           NPairKokkosHalffullNewtonTriDevice,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_SKIP | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,1,1,0> NPairKokkosHalffullNewtonTriHost;
NPairStyle(halffull/newton/tri/skip/ghost/kk/host,
           NPairKokkosHalffullNewtonTriHost,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_SKIP | NP_KOKKOS_HOST);

// Newtoff (can be triclinic but template param always set to 0)

typedef NPairHalffullKokkos<LMPDeviceType,0,0,0> NPairKokkosHalffullNewtoffDevice;
NPairStyle(halffull/newtoff/ghost/kk/device,
           NPairKokkosHalffullNewtoffDevice,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,0,0,0> NPairKokkosHalffullNewtoffHost;
NPairStyle(halffull/newtoff/ghost/kk/host,
           NPairKokkosHalffullNewtoffHost,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_KOKKOS_HOST);

typedef NPairHalffullKokkos<LMPDeviceType,0,0,0> NPairKokkosHalffullNewtoffDevice;
NPairStyle(halffull/newtoff/skip/ghost/kk/device,
           NPairKokkosHalffullNewtoffDevice,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_SKIP | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,0,0,0> NPairKokkosHalffullNewtoffHost;
NPairStyle(halffull/newtoff/skip/ghost/kk/host,
           NPairKokkosHalffullNewtoffHost,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_SKIP | NP_KOKKOS_HOST);

//************ Trim **************

// Newton, no triclinic

typedef NPairHalffullKokkos<LMPDeviceType,1,0,1> NPairKokkosHalffullNewtonTrimDevice;
NPairStyle(halffull/newton/trim/kk/device,
           NPairKokkosHalffullNewtonTrimDevice,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRIM | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,1,0,1> NPairKokkosHalffullNewtonTrimHost;
NPairStyle(halffull/newton/trim/kk/host,
           NPairKokkosHalffullNewtonTrimHost,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRIM | NP_KOKKOS_HOST);

typedef NPairHalffullKokkos<LMPDeviceType,1,0,1> NPairKokkosHalffullNewtonTrimDevice;
NPairStyle(halffull/newton/trim/skip/kk/device,
           NPairKokkosHalffullNewtonTrimDevice,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_SKIP | NP_TRIM | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,1,0,1> NPairKokkosHalffullNewtonTrimHost;
NPairStyle(halffull/newton/trim/skip/kk/host,
           NPairKokkosHalffullNewtonTrimHost,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_SKIP | NP_TRIM | NP_KOKKOS_HOST);

// Newton, triclinic

typedef NPairHalffullKokkos<LMPDeviceType,1,1,1> NPairKokkosHalffullNewtonTriTrimDevice;
NPairStyle(halffull/newton/tri/trim/kk/device,
           NPairKokkosHalffullNewtonTriTrimDevice,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_TRIM | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,1,1,1> NPairKokkosHalffullNewtonTriTrimHost;
NPairStyle(halffull/newton/tri/trim/kk/host,
           NPairKokkosHalffullNewtonTriTrimHost,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_TRIM | NP_KOKKOS_HOST);

typedef NPairHalffullKokkos<LMPDeviceType,1,1,1> NPairKokkosHalffullNewtonTriTrimDevice;
NPairStyle(halffull/newton/tri/trim/skip/kk/device,
           NPairKokkosHalffullNewtonTrimDevice,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_TRIM | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,1,1,1> NPairKokkosHalffullNewtonTriTrimHost;
NPairStyle(halffull/newton/tri/trim/skip/kk/host,
           NPairKokkosHalffullNewtonTriTrimHost,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_TRIM | NP_KOKKOS_HOST);

// Newtoff (can be triclinic but template param always set to 0)

typedef NPairHalffullKokkos<LMPDeviceType,0,0,1> NPairKokkosHalffullNewtoffTrimDevice;
NPairStyle(halffull/newtoff/trim/kk/device,
           NPairKokkosHalffullNewtoffTrimDevice,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_TRIM | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,0,0,1> NPairKokkosHalffullNewtoffTrimHost;
NPairStyle(halffull/newtoff/trim/kk/host,
           NPairKokkosHalffullNewtoffTrimHost,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_TRIM | NP_KOKKOS_HOST);

typedef NPairHalffullKokkos<LMPDeviceType,0,0,1> NPairKokkosHalffullNewtoffTrimDevice;
NPairStyle(halffull/newtoff/trim/skip/kk/device,
           NPairKokkosHalffullNewtoffTrimDevice,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_TRIM | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,0,0,1> NPairKokkosHalffullNewtoffTrimHost;
NPairStyle(halffull/newtoff/trim/skip/kk/host,
           NPairKokkosHalffullNewtoffTrimHost,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_SKIP |  NP_TRIM | NP_KOKKOS_HOST);

//************ Ghost **************

// Newton, no triclinic

typedef NPairHalffullKokkos<LMPDeviceType,1,0,1> NPairKokkosHalffullNewtonTrimDevice;
NPairStyle(halffull/newton/tri/trim/ghost/kk/device,
           NPairKokkosHalffullNewtonTrimDevice,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_GHOST | NP_TRIM | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,1,0,1> NPairKokkosHalffullNewtonTrimHost;
NPairStyle(halffull/newton/trim/ghost/kk/host,
           NPairKokkosHalffullNewtonTrimHost,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_GHOST | NP_TRIM | NP_KOKKOS_HOST);

typedef NPairHalffullKokkos<LMPDeviceType,1,0,1> NPairKokkosHalffullNewtonTrimDevice;
NPairStyle(halffull/newton/trim/skip/ghost/kk/device,
           NPairKokkosHalffullNewtonTrimDevice,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_GHOST | NP_SKIP | NP_TRIM | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,1,0,1> NPairKokkosHalffullNewtonTrimHost;
NPairStyle(halffull/newton/trim/skip/ghost/kk/host,
           NPairKokkosHalffullNewtonTrimHost,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_GHOST | NP_SKIP | NP_TRIM | NP_KOKKOS_HOST);

// Newton, triclinic

typedef NPairHalffullKokkos<LMPDeviceType,1,1,1> NPairKokkosHalffullNewtonTriTrimDevice;
NPairStyle(halffull/newton/tri/trim/ghost/kk/device,
           NPairKokkosHalffullNewtonTriTrimDevice,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_TRIM | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,1,1,1> NPairKokkosHalffullNewtonTriTrimHost;
NPairStyle(halffull/newton/tri/trim/ghost/kk/host,
           NPairKokkosHalffullNewtonTriTrimHost,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_TRIM | NP_KOKKOS_HOST);

typedef NPairHalffullKokkos<LMPDeviceType,1,1,1> NPairKokkosHalffullNewtonTriTrimDevice;
NPairStyle(halffull/newton/tri/trim/skip/ghost/kk/device,
           NPairKokkosHalffullNewtonTriTrimDevice,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_SKIP | NP_TRIM | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,1,1,1> NPairKokkosHalffullNewtonTriTrimHost;
NPairStyle(halffull/newton/tri/trim/skip/ghost/kk/host,
           NPairKokkosHalffullNewtonTriTrimHost,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_SKIP | NP_TRIM | NP_KOKKOS_HOST);

// Newtoff (can be triclinic but template param always set to 0)

typedef NPairHalffullKokkos<LMPDeviceType,0,0,1> NPairKokkosHalffullNewtoffTrimDevice;
NPairStyle(halffull/newtoff/trim/ghost/kk/device,
           NPairKokkosHalffullNewtoffTrimDevice,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_TRIM | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,0,0,1> NPairKokkosHalffullNewtoffTrimHost;
NPairStyle(halffull/newtoff/trim/ghost/kk/host,
           NPairKokkosHalffullNewtoffTrimHost,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_TRIM | NP_KOKKOS_HOST);

typedef NPairHalffullKokkos<LMPDeviceType,0,0,1> NPairKokkosHalffullNewtoffTrimDevice;
NPairStyle(halffull/newtoff/trim/skip/ghost/kk/device,
           NPairKokkosHalffullNewtoffTrimDevice,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_SKIP | NP_TRIM | NP_KOKKOS_DEVICE);

typedef NPairHalffullKokkos<LMPHostType,0,0,1> NPairKokkosHalffullNewtoffTrimHost;
NPairStyle(halffull/newtoff/trim/skip/ghost/kk/host,
           NPairKokkosHalffullNewtoffTrimHost,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_SKIP | NP_TRIM | NP_KOKKOS_HOST);

// clang-format on
#else

// clang-format off
#ifndef LMP_NPAIR_HALFFULL_KOKKOS_H
#define LMP_NPAIR_HALFFULL_KOKKOS_H

#include "npair.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagNPairHalffullCompute{};

template<class DeviceType, int NEWTON, int TRI, int TRIM>
class NPairHalffullKokkos : public NPair {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  NPairHalffullKokkos(class LAMMPS *);
  void build(class NeighList *) override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagNPairHalffullCompute, const int&) const;

 private:
  int nlocal,triclinic;
  double cutsq_custom,delta;

  typename AT::t_x_array_randomread x;

  typename AT::t_neighbors_2d_const d_neighbors_full;
  typename AT::t_int_1d_const d_ilist_full;
  typename AT::t_int_1d_const d_numneigh_full;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d d_ilist;
  typename AT::t_int_1d d_numneigh;
};

}

#endif
#endif

