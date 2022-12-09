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
   Contributing authors: Vladislav Galigerov (HSE),  Vsevolod Nikolskiy (HSE)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(adp/kk,PairADPKokkos<LMPDeviceType>);
PairStyle(adp/kk/device,PairADPKokkos<LMPDeviceType>);
PairStyle(adp/kk/host,PairADPKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_ADP_KOKKOS_H
#define LMP_PAIR_ADP_KOKKOS_H

#include "kokkos_base.h"
#include "pair_kokkos.h"
#include "pair_adp.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

struct TagPairADPPackForwardComm{};
struct TagPairADPUnpackForwardComm{};
struct TagPairADPInitialize{};

template<int NEIGHFLAG, int NEWTON_PAIR>
struct TagPairADPKernelA{};

template<int EFLAG>
struct TagPairADPKernelB{};

template<int EFLAG>
struct TagPairADPKernelAB{};

template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairADPKernelC{};

template<class DeviceType>
class PairADPKokkos : public PairADP, public KokkosBase
{
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairADPKokkos(class LAMMPS *);
  ~PairADPKokkos() override;
  void compute(int, int) override;
  void init_style() override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairADPPackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairADPUnpackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairADPInitialize, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairADPKernelA<NEIGHFLAG,NEWTON_PAIR>, const int&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairADPKernelB<EFLAG>, const int&, EV_FLOAT&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairADPKernelB<EFLAG>, const int&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairADPKernelAB<EFLAG>, const int&, EV_FLOAT&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairADPKernelAB<EFLAG>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairADPKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairADPKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &epair, const F_FLOAT &fx, const F_FLOAT &fy, const F_FLOAT &fz, const F_FLOAT &delx,
                  const F_FLOAT &dely, const F_FLOAT &delz) const;

  int pack_forward_comm_kokkos(int, DAT::tdual_int_2d, int, DAT::tdual_xfloat_1d&,
                       int, int *) override;
  void unpack_forward_comm_kokkos(int, int, DAT::tdual_xfloat_1d&) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 protected:
  typename AT::t_x_array x;
  typename AT::t_f_array f;
  typename AT::t_int_1d type;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename AT::t_efloat_1d d_eatom;
  typename AT::t_virial_array d_vatom;

  int need_dup;

  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template<typename DataType, typename Layout>
  using DupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;

  template<typename DataType, typename Layout>
  using NonDupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<F_FLOAT*, typename DAT::t_ffloat_1d::array_layout> dup_rho;
  DupScatterView<F_FLOAT*[3], typename DAT::t_f_array::array_layout> dup_mu;
  DupScatterView<F_FLOAT*[6], typename DAT::t_virial_array::array_layout> dup_lambda;
  DupScatterView<F_FLOAT*[3], typename DAT::t_f_array::array_layout> dup_f;
  DupScatterView<E_FLOAT*, typename DAT::t_efloat_1d::array_layout> dup_eatom;
  DupScatterView<F_FLOAT*[6], typename DAT::t_virial_array::array_layout> dup_vatom;
  NonDupScatterView<F_FLOAT*, typename DAT::t_ffloat_1d::array_layout> ndup_rho;
  NonDupScatterView<F_FLOAT*[3], typename DAT::t_f_array::array_layout> ndup_mu;
  NonDupScatterView<F_FLOAT*[6], typename DAT::t_virial_array::array_layout> ndup_lambda;
  NonDupScatterView<F_FLOAT*[3], typename DAT::t_f_array::array_layout> ndup_f;
  NonDupScatterView<E_FLOAT*, typename DAT::t_efloat_1d::array_layout> ndup_eatom;
  NonDupScatterView<F_FLOAT*[6], typename DAT::t_virial_array::array_layout> ndup_vatom;

  DAT::tdual_ffloat_1d k_rho;
  DAT::tdual_ffloat_1d k_fp;
  DAT::tdual_f_array k_mu;
  DAT::tdual_virial_array k_lambda;
  typename AT::t_ffloat_1d d_rho;
  typename AT::t_ffloat_1d d_fp;
  typename AT::t_f_array d_mu;
  typename AT::t_virial_array d_lambda;
  HAT::t_ffloat_1d h_rho;
  HAT::t_ffloat_1d h_fp;
  HAT::t_f_array h_mu;
  HAT::t_virial_array h_lambda;

  typename AT::t_int_1d d_type2frho;
  typename AT::t_int_2d_dl d_type2rhor;
  typename AT::t_int_2d_dl d_type2z2r;
  typename AT::t_int_2d_dl d_type2u2r;
  typename AT::t_int_2d_dl d_type2w2r;

  typedef Kokkos::DualView<F_FLOAT**[7],DeviceType> tdual_ffloat_2d_n7;
  typedef typename tdual_ffloat_2d_n7::t_dev_const t_ffloat_2d_n7;
  typedef typename tdual_ffloat_2d_n7::t_host t_host_ffloat_2d_n7;

  t_ffloat_2d_n7 d_frho_spline;
  t_ffloat_2d_n7 d_rhor_spline;
  t_ffloat_2d_n7 d_z2r_spline;
  t_ffloat_2d_n7 d_u2r_spline, d_w2r_spline;

  void interpolate(int, double, double *, t_host_ffloat_2d_n7, int);
  void file2array() override;
  void array2spline() override;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d d_ilist;
  typename AT::t_int_1d d_numneigh;

  int iswap;
  int first;
  typename AT::t_int_2d d_sendlist;
  typename AT::t_xfloat_1d_um v_buf;

  int neighflag,newton_pair;
  int nlocal,nall,eflag,vflag;

  friend void pair_virial_fdotr_compute<PairADPKokkos>(PairADPKokkos*);
};

}
#endif
#endif

