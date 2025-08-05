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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(mliap/property/atom/kk,ComputeMLIAPPropertyAtomKokkos<LMPDeviceType>);
// clang-format on
#else

#ifndef LMP_COMPUTE_MLIAP_PROPERTY_ATOM_KOKKOS_H
#define LMP_COMPUTE_MLIAP_PROPERTY_ATOM_KOKKOS_H

#include "compute.h"
#include "mliap_data_kokkos.h"
#include "pair_mliap_kokkos.h"

namespace LAMMPS_NS {

template <class DeviceType>
class ComputeMLIAPPropertyAtomKokkos : public Compute {
 public:
  ComputeMLIAPPropertyAtomKokkos(class LAMMPS *, int, char **);
  ~ComputeMLIAPPropertyAtomKokkos() override;
  void init() override;
  void compute_peratom() override;

 private:
  int nlocal;
  std::string property_name;
  MLIAPDataKokkos<DeviceType> *k_data;
};

}    // namespace LAMMPS_NS

#endif
#endif
