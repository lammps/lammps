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
ComputeStyle(reaxff/bonds/local/kk,ComputeReaxFFBondsLocalKokkos<LMPDeviceType>);
ComputeStyle(reaxff/bonds/local/kk/device,ComputeReaxFFBondsLocalKokkos<LMPDeviceType>);
ComputeStyle(reaxff/bonds/local/kk/host,ComputeReaxFFBondsLocalKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_COMPUTE_REAXFF_BONDS_LOCAL_KOKKOS_H
#define LMP_COMPUTE_REAXFF_BONDS_LOCAL_KOKKOS_H

#include "compute_reaxff_bonds_local.h"
#include "pair_reaxff_kokkos.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class ComputeReaxFFBondsLocalKokkos : public Compute {
 public:
  using device_type = DeviceType;
  using AT = ArrayTypes<DeviceType>;

  ComputeReaxFFBondsLocalKokkos(class LAMMPS *, int, char **);
  ~ComputeReaxFFBondsLocalKokkos() override;
  void init() override;
  void compute_local() override;
  double memory_usage() override;

 private:
  int nlocal;
  int nvalues;
  int prev_nvalues;

  double **alocal;
  typename AT::tdual_float_2d k_alocal;
  PairReaxFF *reaxff;

  auto device_pair() {
    return dynamic_cast<PairReaxFFKokkos<LMPDeviceType>*>(reaxff);
  }

  auto host_pair() {
    return dynamic_cast<PairReaxFFKokkos<LMPHostType>*>(reaxff);
  }
};

}    // namespace LAMMPS_NS

#endif
#endif
