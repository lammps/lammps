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

#ifdef ATOM_CLASS
// clang-format off
AtomStyle(ellipsoid/kk,AtomVecEllipsoidKokkos);
AtomStyle(ellipsoid/kk/device,AtomVecEllipsoidKokkos);
AtomStyle(ellipsoid/kk/host,AtomVecEllipsoidKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_ATOM_VEC_ELLIPSOID_KOKKOS_H
#define LMP_ATOM_VEC_ELLIPSOID_KOKKOS_H

#include "atom_vec_kokkos.h"
#include "atom_vec_ellipsoid.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

/* ---------------------------------------------------------------------- */
// DualViews for Bonus struct - shape,quat,ilocal.

template <class DeviceType>
struct AtomVecEllipsoidKokkosBonusArray;

template <>
struct AtomVecEllipsoidKokkosBonusArray<LMPDeviceType> {
  typedef Kokkos::
    DualView<AtomVecEllipsoid::Bonus*,
    LMPDeviceType::array_layout,LMPDeviceType> tdual_bonus_1d;
  typedef tdual_bonus_1d::t_dev t_bonus_1d;
  typedef tdual_bonus_1d::t_dev_const_randomread t_bonus_1d_randomread;
};
#ifdef LMP_KOKKOS_GPU
template <>
struct AtomVecEllipsoidKokkosBonusArray<LMPHostType> {
  typedef Kokkos::
    DualView<AtomVecEllipsoid::Bonus*,
    LMPDeviceType::array_layout,LMPDeviceType> tdual_bonus_1d;
  typedef tdual_bonus_1d::t_host t_bonus_1d;
  typedef tdual_bonus_1d::t_host_const_randomread t_bonus_1d_randomread;
};
#endif

typedef AtomVecEllipsoidKokkosBonusArray<LMPDeviceType> DEllipsoidBonusAT;
typedef AtomVecEllipsoidKokkosBonusArray<LMPHostType> HEllipsoidBonusAT;

/* ---------------------------------------------------------------------- */

class AtomVecEllipsoidKokkos : public AtomVecKokkos, public AtomVecEllipsoid {
 public:
  AtomVecEllipsoidKokkos(class LAMMPS *);
  ~AtomVecEllipsoidKokkos() override;
  void init() override;

  void grow(int) override;
  void grow_pointers() override;
  void sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter) override;
  void sync(ExecutionSpace space, uint64_t mask) override;
  void modified(ExecutionSpace space, uint64_t mask) override;
  void sync_pinned(ExecutionSpace space, uint64_t mask, int async_flag = 0) override;

  /* Bonus functions */

  void pack_comm_bonus_kokkos(const int &n, const DAT::tdual_int_1d &list,
                              const DAT::tdual_double_2d_lr &buf, int vel_flag = 0) override;

  void unpack_comm_bonus_kokkos(const int &n, const int &nfirst,
                                const DAT::tdual_double_2d_lr &buf, int vel_flag = 0) override;

  void pack_comm_self_bonus_kokkos(const int &n, const DAT::tdual_int_1d &list,
                                   const int nfirst) override;

  void pack_comm_self_fused_bonus_kokkos(const int &n,
                                         const DAT::tdual_int_2d_lr &list,
                                         const DAT::tdual_int_1d &sendnum_scan,
                                         const DAT::tdual_int_1d &firstrecv,
                                         const DAT::tdual_int_1d &g2l) override;

  void pack_border_bonus_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                                DAT::tdual_double_2d_lr &buf,
                                ExecutionSpace space, int vel_flag = 0) override;
  void unpack_border_bonus_kokkos(const int &n, const int &nfirst,
                                  const DAT::tdual_double_2d_lr &buf,
                                  ExecutionSpace space, int vel_flag = 0) override;

  void pack_exchange_bonus_kokkos(const int &nsend, DAT::tdual_double_2d_lr &buf,
                                  DAT::tdual_int_1d k_sendlist,
                                  DAT::tdual_int_1d k_copylist,
                                  DAT::tdual_int_1d k_copylist_bonus,
                                  ExecutionSpace space) override;

  void unpack_exchange_bonus_kokkos(DAT::tdual_double_2d_lr &k_buf,
                                    int nrecv,
                                    ExecutionSpace space,
                                    DAT::tdual_int_1d &k_indices) override;

  int get_status_nlocal_bonus() override;     // Using these for use in
  void set_status_nlocal_bonus(int) override; // CommKokkos::exchange_device()

  // Bonus struct

  void grow_bonus() override;

  DEllipsoidBonusAT::tdual_bonus_1d k_bonus;
  DEllipsoidBonusAT::t_bonus_1d d_bonus;
  HEllipsoidBonusAT::t_bonus_1d h_bonus;

  void set_size_exchange() override;

 private:
  double **torque;

  DAT::t_tagint_1d d_tag;
  HAT::t_tagint_1d h_tag;
  DAT::t_imageint_1d d_image;
  HAT::t_imageint_1d h_image;
  DAT::t_int_1d d_type, d_mask;
  HAT::t_int_1d h_type, h_mask;

  DAT::t_kkfloat_1d_3_lr d_x;
  DAT::t_kkfloat_1d_3 d_v;
  DAT::t_kkacc_1d_3 d_f;

  DAT::t_kkfloat_1d d_rmass;
  HAT::t_kkfloat_1d h_rmass;
  DAT::t_kkfloat_1d_3 d_angmom;
  HAT::t_kkfloat_1d_3 h_angmom;
  DAT::t_kkacc_1d_3 d_torque;
  HAT::t_kkacc_1d_3 h_torque;
  DAT::t_int_1d d_ellipsoid;
  HAT::t_int_1d h_ellipsoid;

  DAT::tdual_int_scalar k_nghost_bonus;
  DAT::tdual_int_scalar k_nlocal_bonus;
};

}    // namespace LAMMPS_NS

#endif
#endif
