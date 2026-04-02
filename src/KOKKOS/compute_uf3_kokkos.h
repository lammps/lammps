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
ComputeStyle(uf3/kk, ComputeUF3Kokkos<LMPDeviceType>);
ComputeStyle(uf3/kk/device, ComputeUF3Kokkos<LMPDeviceType>);
ComputeStyle(uf3/kk/host, ComputeUF3Kokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_COMPUTE_UF3_KOKKOS_H
#define LMP_COMPUTE_UF3_KOKKOS_H

#include "compute.h"
#include "kokkos.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template <class DeviceType> class ComputeUF3Kokkos : public Compute {
 public:
  typedef DeviceType device_type;

  ComputeUF3Kokkos(class LAMMPS *, int, char **);
  ~ComputeUF3Kokkos() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_array() override;
  double memory_usage() override;

 private:
  class PairUF3 *uf3_reader;
  class AtomKokkos *atomKK;
  double cutmax;
  int lastcol;
  int ncoeff;
  int **col_off_2b;
  class NeighList *list;
  class Compute *c_pe, *c_virial;
  std::string id_virial;
  double **uf3local;
  double **uf3all;

  void free_local();
  void build_column_offsets();
};

}    // namespace LAMMPS_NS
#endif
#endif
