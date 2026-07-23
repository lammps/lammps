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

/* ----------------------------------------------------------------------
   lj/cut/coul/long2/kk -- experimental hand-written GPU variant of
   lj/cut/coul/long/kk.  Physics is identical; the difference is the
   force-only compute path, which bypasses the generic pair_kokkos.h
   dispatch (PairComputeFunctor / ScatterView / CoulLongTable templating)
   in favor of one lean kernel: a single fused loop with analytical Ewald
   (no coul table lookup), register-resident i-force accumulation, and
   direct atomic j-scatter.  Energy/virial (thermo) steps delegate to the
   proven base-class kernel, so only the hot force-only path is new.

   Scope: GPU, single pair style, half list + newton on (-> HALFTHREAD),
   ntypes <= MAX_TYPES_STACKPARAMS.  Any other configuration transparently
   falls back to PairLJCutCoulLongKokkos::compute().  See kokkos_neigh.md.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/coul/long2/kk,PairLJCutCoulLong2Kokkos<LMPDeviceType>);
PairStyle(lj/cut/coul/long2/kk/device,PairLJCutCoulLong2Kokkos<LMPDeviceType>);
PairStyle(lj/cut/coul/long2/kk/host,PairLJCutCoulLong2Kokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LJ_CUT_COUL_LONG2_KOKKOS_H
#define LMP_PAIR_LJ_CUT_COUL_LONG2_KOKKOS_H

#include "pair_lj_cut_coul_long_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairLJCutCoulLong2Kokkos : public PairLJCutCoulLongKokkos<DeviceType> {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  PairLJCutCoulLong2Kokkos(class LAMMPS *);

  void compute(int, int) override;
};

}

#endif
#endif
