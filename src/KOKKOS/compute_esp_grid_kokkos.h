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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(esp/grid/kk,ComputeESPGridKokkos<LMPDeviceType>);
ComputeStyle(esp/grid/kk/device,ComputeESPGridKokkos<LMPDeviceType>);
ComputeStyle(esp/grid/kk/host,ComputeESPGridKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_COMPUTE_ESP_GRID_KOKKOS_H
#define LMP_COMPUTE_ESP_GRID_KOKKOS_H

#include "compute_esp_grid.h"
#include "kokkos_type.h"
#include "grid3d.h"

namespace LAMMPS_NS {

template<class DeviceType>
class ComputeESPGridKokkos : public ComputeESPGrid {
 public:
  ComputeESPGridKokkos(class LAMMPS *, int, char **);
  ~ComputeESPGridKokkos() override;

  void init() override;
  void compute_pergrid() override;
  double compute_scalar() override;

 private:
  // typedefs for managing Kokkos data
  typedef ArrayTypes<DeviceType> AT;
  
  // DualViews for grid data
  typedef Kokkos::DualView<double***, Kokkos::LayoutRight, DeviceType> tdual_3d_double;
  typedef Kokkos::DualView<double**, Kokkos::LayoutRight, DeviceType> tdual_2d_double;
  typedef Kokkos::DualView<double*, DeviceType> tdual_1d_double;
  typedef Kokkos::DualView<int*, DeviceType> tdual_1d_int;
  
  tdual_3d_double k_esp;          // ESP values on grid
  tdual_3d_double k_weight;       // Weight values on grid
  tdual_3d_double k_reference;    // Reference values on grid
  tdual_1d_double k_bcut_acks2;   // ReaxFF cutoff values
  
  // Device views for direct kernel access
  typename AT::t_x_array d_x;     // atom positions
  typename AT::t_efloat_1d d_q;   // atom charges
  typename AT::t_int_1d d_type;   // atom types

  KOKKOS_INLINE_FUNCTION
  double compute_weight(double r, double rcut) const;
};

}    // namespace LAMMPS_NS

#endif
#endif