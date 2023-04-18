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

#ifdef FIX_CLASS
// clang-format off
FixStyle(shake/kk,FixShakeKokkos<LMPDeviceType>);
FixStyle(shake/kk/device,FixShakeKokkos<LMPDeviceType>);
FixStyle(shake/kk/host,FixShakeKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_SHAKE_KOKKOS_H
#define LMP_FIX_SHAKE_KOKKOS_H

#include "fix_shake.h"
#include "kokkos_type.h"
#include "kokkos_base.h"
#include <Kokkos_UnorderedMap.hpp>

namespace LAMMPS_NS {

struct TagFixShakePreNeighbor{};

template<int NEIGHFLAG, int EVFLAG>
struct TagFixShakePostForce{};

template<int PBC_FLAG>
struct TagFixShakePackForwardComm{};

struct TagFixShakeUnpackForwardComm{};
struct TagFixShakeUnpackExchange{};

template<class DeviceType>
class FixShakeKokkos : public FixShake, public KokkosBase {

 //friend class FixEHEX;

 public:
  typedef DeviceType device_type;
  typedef EV_FLOAT value_type;
  typedef ArrayTypes<DeviceType> AT;

  FixShakeKokkos(class LAMMPS *, int, char **);
  ~FixShakeKokkos() override;
  void init() override;
  void min_setup(int) override;
  void pre_neighbor() override;
  void post_force(int) override;
  void min_post_force(int) override;

  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  void update_arrays(int, int) override;
  void set_molecule(int, tagint, int, double *, double *, double *) override;

  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_forward_comm_kokkos(int, DAT::tdual_int_2d, int, DAT::tdual_xfloat_1d&,
                       int, int *) override;
  void unpack_forward_comm_kokkos(int, int, DAT::tdual_xfloat_1d&) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

  void shake_end_of_step(int vflag) override;
  void correct_coordinates(int vflag) override;

  int dof(int) override;

  void unconstrained_update() override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixShakePreNeighbor, const int&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixShakePostForce<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixShakePostForce<NEIGHFLAG,EVFLAG>, const int&) const;

  template<int PBC_FLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixShakePackForwardComm<PBC_FLAG>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixShakeUnpackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void pack_exchange_item(const int&, int &, const bool &) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixShakeUnpackExchange, const int&) const;

  int pack_exchange_kokkos(const int &nsend,DAT::tdual_xfloat_2d &buf,
                           DAT::tdual_int_1d k_sendlist,
                           DAT::tdual_int_1d k_copylist,
                           ExecutionSpace space) override;

  void unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf,
                              DAT::tdual_int_1d &indices,int nrecv,
                              ExecutionSpace space) override;

 protected:
  typename AT::t_x_array d_x;
  typename AT::t_v_array d_v;
  typename AT::t_f_array d_f;
  typename AT::t_float_1d d_rmass;
  typename AT::t_float_1d d_mass;
  typename AT::t_tagint_1d_randomread d_tag;
  typename AT::t_int_1d d_type;
  typename AT::t_int_1d d_mask;

  DAT::tdual_efloat_1d k_eatom;
  typename AT::t_efloat_1d d_eatom;

  DAT::tdual_virial_array k_vatom;
  typename AT::t_virial_array d_vatom;

  DAT::tdual_float_1d k_bond_distance; // constraint distances
  typename AT::t_float_1d d_bond_distance;
  DAT::tdual_float_1d k_angle_distance;
  typename AT::t_float_1d d_angle_distance;

                                         // atom-based arrays
  DAT::tdual_int_1d k_shake_flag;
  typename AT::t_int_1d d_shake_flag; // 0 if atom not in SHAKE cluster
                                         // 1 = size 3 angle cluster
                                         // 2,3,4 = size of bond-only cluster
  DAT::tdual_tagint_2d k_shake_atom;
  typename AT::t_tagint_2d d_shake_atom; // global IDs of atoms in cluster
                                         // central atom is 1st
                                         // lowest global ID is 1st for size 2
  DAT::tdual_int_2d k_shake_type;
  typename AT::t_int_2d d_shake_type; // bondtype of each bond in cluster
                                         // for angle cluster, 3rd value
                                         //   is angletype
  DAT::tdual_x_array k_xshake;
  typename AT::t_x_array d_xshake; // unconstrained atom coords

  DAT::tdual_int_1d k_list;
  typename AT::t_int_1d d_list; // list of clusters to SHAKE

  DAT::tdual_int_2d k_closest_list;
  typename AT::t_int_2d d_closest_list; // list of closest atom indices in SHAKE clusters

  DAT::tdual_int_scalar k_error_flag;
  DAT::tdual_int_scalar k_nlist;

  typename AT::t_int_scalar d_count;
  HAT::t_int_scalar h_count;

  void stats() override;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void shake(int, EV_FLOAT&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void shake3(int, EV_FLOAT&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void shake4(int, EV_FLOAT&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void shake3angle(int, EV_FLOAT&) const;

  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template<typename DataType, typename Layout>
  using DupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;

  template<typename DataType, typename Layout>
  using NonDupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<F_FLOAT*[3], typename DAT::t_f_array::array_layout> dup_f;
  DupScatterView<E_FLOAT*, typename DAT::t_efloat_1d::array_layout> dup_eatom;
  DupScatterView<F_FLOAT*[6], typename DAT::t_virial_array::array_layout> dup_vatom;

  NonDupScatterView<F_FLOAT*[3], typename DAT::t_f_array::array_layout> ndup_f;
  NonDupScatterView<E_FLOAT*, typename DAT::t_efloat_1d::array_layout> ndup_eatom;
  NonDupScatterView<F_FLOAT*[6], typename DAT::t_virial_array::array_layout> ndup_vatom;

  int neighflag,need_dup;

  typename AT::t_int_1d d_scalars;
  HAT::t_int_1d h_scalars;
  typename AT::t_int_scalar d_error_flag;
  typename AT::t_int_scalar d_nlist;
  HAT::t_int_scalar h_error_flag;
  HAT::t_int_scalar h_nlist;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void v_tally(EV_FLOAT&, int, int *, double, double *) const;

  int iswap,first,nsend;

  typename AT::t_int_2d d_sendlist;
  typename AT::t_xfloat_1d_um d_buf;

  typename AT::t_int_1d d_exchange_sendlist;
  typename AT::t_int_1d d_copylist;
  typename AT::t_int_1d d_indices;

  X_FLOAT dx,dy,dz;

  int *shake_flag_tmp;
  tagint **shake_atom_tmp;
  int **shake_type_tmp;

  DAT::tdual_int_1d k_sametag;
  typename AT::t_int_1d d_sametag;
  int map_style;
  DAT::tdual_int_1d k_map_array;
  dual_hash_type k_map_hash;

  // copied from Domain

  KOKKOS_INLINE_FUNCTION
  int closest_image(const int, int) const;

  int triclinic;
  int xperiodic,yperiodic,zperiodic;
  X_FLOAT xprd_half,yprd_half,zprd_half;
  X_FLOAT xprd,yprd,zprd;
  X_FLOAT xy,xz,yz;
};

template <class DeviceType>
struct FixShakeKokkosPackExchangeFunctor {
  typedef DeviceType device_type;
  typedef int value_type;
  FixShakeKokkos<DeviceType> c;
  FixShakeKokkosPackExchangeFunctor(FixShakeKokkos<DeviceType>* c_ptr):c(*c_ptr) {};
  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i, int &offset, const bool &final) const {
    c.pack_exchange_item(i, offset, final);
  }
};

}

#endif
#endif

