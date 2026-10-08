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

// writes stamp to moved if any atom moved more than sqrt(deltasq) since the last build
template<class DeviceType>
struct NeighborCheckDistanceFlag {
  Kokkos::View<int,LMPPinnedHostType> moved;    // pinned host memory, written by the device
  typename ArrayTypes<DeviceType>::t_kkfloat_1d_3_lr x, xhold;
  double deltasq;
  int stamp;
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    const double delx = static_cast<double>(x(i,0) - xhold(i,0));
    const double dely = static_cast<double>(x(i,1) - xhold(i,1));
    const double delz = static_cast<double>(x(i,2) - xhold(i,2));
    if (delx*delx + dely*dely + delz*delz > deltasq) moved() = stamp;
  }
};

class NeighborKokkos : public Neighbor {
 public:
  typedef int value_type;

  NeighborKokkos(class LAMMPS *);
  ~NeighborKokkos() override;
  void init() override;
  void init_topology() override;
  void build_topology() override;

  template<class DeviceType>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagNeighborCheckDistance<DeviceType>, const int&, int&) const;

  template<class DeviceType>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagNeighborXhold<DeviceType>, const int&) const;

  DAT::ttransform_kkfloat_2d k_cutneighsq;
  DAT::ttransform_kkfloat_2d k_cutneighghostsq;

  DAT::tdual_int_1d k_ex1_type,k_ex2_type;
  DAT::ttransform_int_2d k_ex_type;
  DAT::tdual_int_1d k_ex1_bit,k_ex2_bit;
  DAT::tdual_int_1d k_ex_mol_group;
  DAT::tdual_int_1d k_ex_mol_bit;
  DAT::tdual_int_1d k_ex_mol_intra;

  NeighBondKokkos<LMPHostType> neighbond_host;
  NeighBondKokkos<LMPDeviceType> neighbond_device;

  DAT::tdual_int_2d_lr k_bondlist;
  DAT::tdual_int_2d_lr k_anglelist;
  DAT::tdual_int_2d_lr k_dihedrallist;
  DAT::tdual_int_2d_lr k_improperlist;

  int device_flag;

  int check_distance() override;
  void build(int) override;

 private:

  DAT::ttransform_kkfloat_1d_3_lr x;
  DAT::ttransform_kkfloat_1d_3_lr xhold;

  double deltasq;

  // flag of check_distance(): atoms that moved too far write the stamp of
  // the current call, so the flag needs no reset between calls
  Kokkos::View<int,LMPPinnedHostType> h_moved;
  int moved_stamp;

  void init_cutneighsq_kokkos(int) override;
  void init_cutneighghostsq_kokkos(int) override;
  void create_kokkos_list(int) override;
  void init_ex_type_kokkos(int) override;
  void init_ex_bit_kokkos() override;
  void init_ex_mol_bit_kokkos() override;
  void grow_ex_mol_intra_kokkos() override;
  template<class DeviceType> int check_distance_kokkos();
  template<class DeviceType> void build_kokkos(int);
  void modify_ex_type_grow_kokkos();
  void modify_ex_group_grow_kokkos();
  void modify_mol_group_grow_kokkos();
  void modify_mol_intra_grow_kokkos();
  void set_binsize_kokkos() override;
};

}

#endif

