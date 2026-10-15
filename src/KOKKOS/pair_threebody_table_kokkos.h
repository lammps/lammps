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
PairStyle(threebody/table/kk,PairThreebodyTableKokkos<LMPDeviceType>);
PairStyle(threebody/table/kk/device,PairThreebodyTableKokkos<LMPDeviceType>);
PairStyle(threebody/table/kk/host,PairThreebodyTableKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_THREEBODY_TABLE_KOKKOS_H
#define LMP_PAIR_THREEBODY_TABLE_KOKKOS_H

#include "pair_threebody_table.h"
#include "pair_kokkos.h"

template<int NEIGHFLAG, int EVFLAG>
struct TagPairThreebodyTableCompute{};

struct TagPairThreebodyTableComputeShortNeigh{};

namespace LAMMPS_NS {

template<class DeviceType>
class PairThreebodyTableKokkos : public PairThreebodyTable {
 public:
  enum {EnabledNeighFlags=HALF|HALFTHREAD};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairThreebodyTableKokkos(class LAMMPS *);
  ~PairThreebodyTableKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;

  template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairThreebodyTableCompute<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairThreebodyTableCompute<NEIGHFLAG,EVFLAG>, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairThreebodyTableComputeShortNeigh, const int&) const;

  template<int NEIGHFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally3(EV_FLOAT &ev, const int &i, const int &j, int &k,
            const KK_FLOAT &evdwl, const KK_FLOAT &ecoul,
                       KK_ACC_FLOAT *fj, KK_ACC_FLOAT *fk, KK_FLOAT *drji, KK_FLOAT *drki) const;

 protected:
  typename AT::t_int_3d_randomread d_elem3param;
  typename AT::t_int_1d_randomread d_map;

  // parameters of an I-J-K interaction needed for the table lookup

  struct TableParamKokkos {
    KK_FLOAT cutsq, rmin, rmax, dr, dtheta;
    int ninput, symmetric;
  };

  typedef Kokkos::DualView<TableParamKokkos*,DeviceType> tdual_param_1d;
  typedef typename tdual_param_1d::t_dev t_param_1d;
  typedef typename tdual_param_1d::t_host t_host_param_1d;

  t_param_1d d_params;

  // the six force factors and the energy of all tables, padded to the longest table

  typedef Kokkos::DualView<KK_FLOAT**[7],Kokkos::LayoutRight,DeviceType> tdual_table_2d_n7;
  typename tdual_table_2d_n7::t_dev_const_randomread d_tables;

  void setup_params() override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void threebody_kk(const int&, const KK_FLOAT&, const KK_FLOAT&, KK_FLOAT *, KK_FLOAT *,
                    KK_ACC_FLOAT *, KK_ACC_FLOAT *, KK_ACC_FLOAT *, const int&, KK_FLOAT&) const;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread type;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  int need_dup;

  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template<typename DataType, typename Layout>
  using DupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;

  template<typename DataType, typename Layout>
  using NonDupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> dup_f;
  DupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> dup_eatom;
  DupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> dup_vatom;

  NonDupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> ndup_f;
  NonDupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> ndup_eatom;
  NonDupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> ndup_vatom;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  int neighflag,newton_pair;
  int nlocal,nall,eflag,vflag;

  int inum;
  typename AT::t_int_2d_dl d_neighbors_short;
  typename AT::t_int_1d d_numneigh_short;

  friend void pair_virial_fdotr_compute<PairThreebodyTableKokkos>(PairThreebodyTableKokkos*);
};

}

#endif
#endif
