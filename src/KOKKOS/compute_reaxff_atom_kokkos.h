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
   Contributing author: Richard Berger (LANL)
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(reaxff/atom/kk,ComputeReaxFFAtomKokkos<LMPDeviceType>);
ComputeStyle(reaxff/atom/kk/device,ComputeReaxFFAtomKokkos<LMPDeviceType>);
ComputeStyle(reaxff/atom/kk/host,ComputeReaxFFAtomKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_COMPUTE_REAXFF_BONDS_KOKKOS_H
#define LMP_COMPUTE_REAXFF_BONDS_KOKKOS_H

#include "compute_reaxff_atom.h"
#include "pair_reaxff_kokkos.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class ComputeReaxFFAtomKokkos : public ComputeReaxFFAtom {
 public:
  using device_type = DeviceType;
  using AT = ArrayTypes<DeviceType>;

  ComputeReaxFFAtomKokkos(class LAMMPS *, int, char **);
  ~ComputeReaxFFAtomKokkos() override;
  void init() override;
  void compute_local() override;
  void compute_peratom() override;
  void compute_bonds() override;
  double memory_usage() override;

 private:
  int nbuf;
  double *buf;
  typename AT::tdual_float_1d k_buf;

  auto device_pair() {
    return static_cast<PairReaxFFKokkos<LMPDeviceType>*>(reaxff);
  }

  auto host_pair() {
    return static_cast<PairReaxFFKokkos<LMPHostType>*>(reaxff);
  }
};

}    // namespace LAMMPS_NS

#endif
#endif
