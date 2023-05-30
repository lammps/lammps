// clang-format off
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

#ifndef LMP_NEIGHBOR_KOKKOS_H
#define LMP_NEIGHBOR_KOKKOS_H

#include "neighbor.h"           // IWYU pragma: export
#include "neigh_list_kokkos.h"
#include "neigh_bond_kokkos.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
struct TagNeighborCheckDistance{};

template<class DeviceType>
struct TagNeighborXhold{};

class NeighborKokkos : public Neighbor {
 public:
  typedef int value_type;

  NeighborKokkos(class LAMMPS *);
  ~NeighborKokkos() override;
  void init() override;
  void init_topology() override;
  void build_topology() override;

  template<class DeviceType>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagNeighborCheckDistance<DeviceType>, const int&, int&) const;

  template<class DeviceType>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagNeighborXhold<DeviceType>, const int&) const;

  DAT::tdual_xfloat_2d k_cutneighsq;

  DAT::tdual_int_1d k_ex1_type,k_ex2_type;
  DAT::tdual_int_2d k_ex_type;
  DAT::tdual_int_1d k_ex1_bit,k_ex2_bit;
  DAT::tdual_int_1d k_ex_mol_group;
  DAT::tdual_int_1d k_ex_mol_bit;
  DAT::tdual_int_1d k_ex_mol_intra;

  NeighBondKokkos<LMPHostType> neighbond_host;
  NeighBondKokkos<LMPDeviceType> neighbond_device;

  DAT::tdual_int_2d k_bondlist;
  DAT::tdual_int_2d k_anglelist;
  DAT::tdual_int_2d k_dihedrallist;
  DAT::tdual_int_2d k_improperlist;

  int device_flag;

 private:

  DAT::tdual_x_array x;
  DAT::tdual_x_array xhold;

  X_FLOAT deltasq;

  void init_cutneighsq_kokkos(int) override;
  void create_kokkos_list(int) override;
  void init_ex_type_kokkos(int) override;
  void init_ex_bit_kokkos() override;
  void init_ex_mol_bit_kokkos() override;
  void grow_ex_mol_intra_kokkos() override;
  int check_distance() override;
  template<class DeviceType> int check_distance_kokkos();
  void build(int) override;
  template<class DeviceType> void build_kokkos(int);
  void setup_bins_kokkos(int);
  void modify_ex_type_grow_kokkos();
  void modify_ex_group_grow_kokkos();
  void modify_mol_group_grow_kokkos();
  void modify_mol_intra_grow_kokkos();
  void set_binsize_kokkos() override;
};

}

#endif

