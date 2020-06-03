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

#ifdef NPAIR_CLASS

// Newton

typedef NPairHalffullKokkos<Device,1> NPairKokkosHalffullNewtonDevice;
NPairStyle(halffull/newton/kk/device,
           NPairKokkosHalffullNewtonDevice,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_KOKKOS_DEVICE)

typedef NPairHalffullKokkos<Host,1> NPairKokkosHalffullNewtonHost;
NPairStyle(halffull/newton/kk/host,
           NPairKokkosHalffullNewtonHost,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_KOKKOS_HOST)

typedef NPairHalffullKokkos<Device,1> NPairKokkosHalffullNewtonDevice;
NPairStyle(halffull/newton/skip/kk/device,
           NPairKokkosHalffullNewtonDevice,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_KOKKOS_DEVICE)

typedef NPairHalffullKokkos<Host,1> NPairKokkosHalffullNewtonHost;
NPairStyle(halffull/newton/skip/kk/host,
           NPairKokkosHalffullNewtonHost,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_KOKKOS_HOST)

// Newtoff

typedef NPairHalffullKokkos<Device,0> NPairKokkosHalffullNewtoffDevice;
NPairStyle(halffull/newtoff/kk/device,
           NPairKokkosHalffullNewtoffDevice,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_KOKKOS_DEVICE)

typedef NPairHalffullKokkos<Host,0> NPairKokkosHalffullNewtoffHost;
NPairStyle(halffull/newtoff/kk/host,
           NPairKokkosHalffullNewtoffHost,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_KOKKOS_HOST)

typedef NPairHalffullKokkos<Device,0> NPairKokkosHalffullNewtoffDevice;
NPairStyle(halffull/newtoff/skip/kk/device,
           NPairKokkosHalffullNewtoffDevice,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_KOKKOS_DEVICE)

typedef NPairHalffullKokkos<Host,0> NPairKokkosHalffullNewtoffHost;
NPairStyle(halffull/newtoff/skip/kk/host,
           NPairKokkosHalffullNewtoffHost,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_KOKKOS_HOST)

//************ Ghost **************

// Newton

typedef NPairHalffullKokkos<Device,1> NPairKokkosHalffullNewtonGhostDevice;
NPairStyle(halffull/newton/ghost/kk/device,
           NPairKokkosHalffullNewtonGhostDevice,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_KOKKOS_DEVICE)

typedef NPairHalffullKokkos<Host,1> NPairKokkosHalffullNewtonHost;
NPairStyle(halffull/newton/ghost/kk/host,
           NPairKokkosHalffullNewtonHost,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_KOKKOS_HOST)

typedef NPairHalffullKokkos<Device,1> NPairKokkosHalffullNewtonGhostDevice;
NPairStyle(halffull/newton/skip/ghost/kk/device,
           NPairKokkosHalffullNewtonGhostDevice,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_SKIP | NP_KOKKOS_DEVICE)

typedef NPairHalffullKokkos<Host,1> NPairKokkosHalffullNewtonHost;
NPairStyle(halffull/newton/skip/ghost/kk/host,
           NPairKokkosHalffullNewtonHost,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_SKIP | NP_KOKKOS_HOST)

// Newtoff

typedef NPairHalffullKokkos<Device,0> NPairKokkosHalffullNewtoffGhostDevice;
NPairStyle(halffull/newtoff/ghost/kk/device,
           NPairKokkosHalffullNewtoffGhostDevice,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_KOKKOS_DEVICE)

typedef NPairHalffullKokkos<Host,0> NPairKokkosHalffullNewtoffHost;
NPairStyle(halffull/newtoff/ghost/kk/host,
           NPairKokkosHalffullNewtoffHost,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_KOKKOS_HOST)

typedef NPairHalffullKokkos<Device,0> NPairKokkosHalffullNewtoffGhostDevice;
NPairStyle(halffull/newtoff/skip/ghost/kk/device,
           NPairKokkosHalffullNewtoffGhostDevice,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_SKIP | NP_KOKKOS_DEVICE)

typedef NPairHalffullKokkos<Host,0> NPairKokkosHalffullNewtoffHost;
NPairStyle(halffull/newtoff/skip/ghost/kk/host,
           NPairKokkosHalffullNewtoffHost,
           NP_HALF_FULL | NP_NEWTOFF | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_SKIP | NP_KOKKOS_HOST)

#else

#ifndef LMP_NPAIR_HALFFULL_KOKKOS_H
#define LMP_NPAIR_HALFFULL_KOKKOS_H

#include "npair.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagNPairHalffullCompute{};

template<ExecutionSpace Space, int NEWTON>
class NPairHalffullKokkos : public NPair {
 public:
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;

  NPairHalffullKokkos(class LAMMPS *);
  ~NPairHalffullKokkos() {}
  void build(class NeighList *);

  KOKKOS_INLINE_FUNCTION
  void operator()(TagNPairHalffullCompute, const int&) const;

 private:
  int nlocal;

  typename AT::t_float_1d_3_randomread x;

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

/* ERROR/WARNING messages:

*/
