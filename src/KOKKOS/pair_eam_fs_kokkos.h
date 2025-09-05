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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(eam/fs/kk,PairEAMFSKokkos<LMPDeviceType>);
PairStyle(eam/fs/kk/device,PairEAMFSKokkos<LMPDeviceType>);
PairStyle(eam/fs/kk/host,PairEAMFSKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_EAM_FS_KOKKOS_H
#define LMP_PAIR_EAM_FS_KOKKOS_H

#include "kokkos_base.h"
#include "pair_kokkos.h"
#include "pair_eam.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

struct TagPairEAMFSPackForwardComm{};
struct TagPairEAMFSUnpackForwardComm{};
struct TagPairEAMFSInitialize{};

template<int NEIGHFLAG, int NEWTON_PAIR>
struct TagPairEAMFSKernelA{};

template<int EFLAG>
struct TagPairEAMFSKernelB{};

template<int EFLAG>
struct TagPairEAMFSKernelAB{};

template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairEAMFSKernelC{};

// Cannot use virtual inheritance on the GPU

template<class DeviceType>
class PairEAMFSKokkos : public PairEAM, public KokkosBase {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairEAMFSKokkos(class LAMMPS *);
  ~PairEAMFSKokkos() override;
  void compute(int, int) override;
  void init_style() override;
  void coeff(int, char **) override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMFSPackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMFSUnpackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMFSInitialize, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMFSKernelA<NEIGHFLAG,NEWTON_PAIR>, const int&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMFSKernelB<EFLAG>, const int&, EV_FLOAT&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMFSKernelB<EFLAG>, const int&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMFSKernelAB<EFLAG>, const int&, EV_FLOAT&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMFSKernelAB<EFLAG>, const int&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMFSKernelAB<EFLAG>, const typename Kokkos::TeamPolicy<DeviceType>::member_type&, EV_FLOAT&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMFSKernelAB<EFLAG>, const typename Kokkos::TeamPolicy<DeviceType>::member_type&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMFSKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMFSKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMFSKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const typename Kokkos::TeamPolicy<DeviceType>::member_type&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMFSKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const typename Kokkos::TeamPolicy<DeviceType>::member_type&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                  const KK_FLOAT &dely, const KK_FLOAT &delz) const;

  int pack_forward_comm_kokkos(int, DAT::tdual_int_1d, DAT::tdual_double_1d&,
                       int, int *) override;
  void unpack_forward_comm_kokkos(int, int, DAT::tdual_double_1d&) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 protected:
  typename AT::t_kkfloat_1d_3_lr x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d type;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  int need_dup,inum;

  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template<typename DataType, typename Layout>
  using DupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;

  template<typename DataType, typename Layout>
  using NonDupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<KK_FLOAT*, typename DAT::t_kkfloat_1d::array_layout> dup_rho;
  DupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> dup_f;
  DupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> dup_eatom;
  DupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> dup_vatom;
  NonDupScatterView<KK_FLOAT*, typename DAT::t_kkfloat_1d::array_layout> ndup_rho;
  NonDupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> ndup_f;
  NonDupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> ndup_eatom;
  NonDupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> ndup_vatom;

  DAT::tdual_kkfloat_1d k_rho;
  DAT::tdual_kkfloat_1d k_fp;
  typename AT::t_kkfloat_1d d_rho;
  typename AT::t_kkfloat_1d d_fp;
  HAT::t_kkfloat_1d h_rho;
  HAT::t_kkfloat_1d h_fp;

  typename AT::t_int_1d d_type2frho;
  typename AT::t_int_2d_dl d_type2rhor;
  typename AT::t_int_2d_dl d_type2z2r;

  typedef Kokkos::DualView<KK_FLOAT**[7],DeviceType> tdual_kkfloat_2d_n7;
  typedef typename tdual_kkfloat_2d_n7::t_dev_const t_kkfloat_2d_n7;
  typedef typename tdual_kkfloat_2d_n7::t_host t_hostkkfloat_2d_n7;

  t_kkfloat_2d_n7 d_frho_spline;
  t_kkfloat_2d_n7 d_rhor_spline;
  t_kkfloat_2d_n7 d_z2r_spline;
  void interpolate(int, double, double *, t_hostkkfloat_2d_n7, int);
  void file2array() override;
  void file2array_fs();
  void array2spline() override;
  void read_file(char *) override;

  template<class TAG>
  struct policyInstance;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d d_ilist;
  typename AT::t_int_1d d_numneigh;

  int first;
  typename AT::t_int_1d d_sendlist;
  typename AT::t_double_1d_um v_buf;

  int neighflag,newton_pair;
  int nlocal,nall,eflag,vflag;

  friend void pair_virial_fdotr_compute<PairEAMFSKokkos>(PairEAMFSKokkos*);
};

}
#endif
#endif

