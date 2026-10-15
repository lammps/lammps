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
PairStyle(sw/angle/table/kk,PairSWAngleTableKokkos<LMPDeviceType>);
PairStyle(sw/angle/table/kk/device,PairSWAngleTableKokkos<LMPDeviceType>);
PairStyle(sw/angle/table/kk/host,PairSWAngleTableKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_SW_ANGLE_TABLE_KOKKOS_H
#define LMP_PAIR_SW_ANGLE_TABLE_KOKKOS_H

#include "pair_sw_angle_table.h"
#include "pair_kokkos.h"

template<int NEIGHFLAG, int EVFLAG>
struct TagPairSWAngleTableCompute{};

struct TagPairSWAngleTableComputeShortNeigh{};

namespace LAMMPS_NS {

template<class DeviceType>
class PairSWAngleTableKokkos : public PairSWAngleTable {
 public:
  enum {EnabledNeighFlags=HALF|HALFTHREAD};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairSWAngleTableKokkos(class LAMMPS *);
  ~PairSWAngleTableKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;

  template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairSWAngleTableCompute<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairSWAngleTableCompute<NEIGHFLAG,EVFLAG>, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairSWAngleTableComputeShortNeigh, const int&) const;

  template<int NEIGHFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                  const KK_FLOAT &dely, const KK_FLOAT &delz) const;

  template<int NEIGHFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally3(EV_FLOAT &ev, const int &i, const int &j, int &k,
            const KK_FLOAT &evdwl, const KK_FLOAT &ecoul,
                       KK_ACC_FLOAT *fj, KK_ACC_FLOAT *fk, KK_FLOAT *drji, KK_FLOAT *drki) const;

 protected:
  typename AT::t_int_3d_randomread d_elem3param;
  typename AT::t_int_1d_randomread d_map;

  typedef Kokkos::DualView<Param*,DeviceType> tdual_param_1d;
  typedef typename tdual_param_1d::t_dev t_param_1d;
  typedef typename tdual_param_1d::t_host t_host_param_1d;

  t_param_1d d_params;

  void setup_params() override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void twobody(const Param&, const KK_FLOAT&, KK_FLOAT&, const int&, KK_FLOAT&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void threebody_table_kk(const Param&, const Param&, const int&, const KK_FLOAT&, const KK_FLOAT&, KK_FLOAT *, KK_FLOAT *,
                          KK_ACC_FLOAT *, KK_ACC_FLOAT *, const int&, KK_FLOAT&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void uf_lookup_kk(const int&, const KK_FLOAT&, KK_FLOAT&, KK_FLOAT&) const;

  // angle tables of all I-J-K parameter sets, padded to the longest table

  struct TableParamKokkos {
    int tabstyle, tablength;
    KK_FLOAT invdelta, deltasq6;
  };

  typedef Kokkos::DualView<TableParamKokkos*,DeviceType> tdual_tableparam_1d;
  typename tdual_tableparam_1d::t_dev d_table_params;

  DAT::tdual_kkfloat_2d k_ang, k_e, k_de, k_f, k_df, k_e2, k_f2;
  typename AT::t_kkfloat_2d_randomread d_ang, d_e, d_de, d_f, d_df, d_e2, d_f2;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_tagint_1d tag;
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


  friend void pair_virial_fdotr_compute<PairSWAngleTableKokkos>(PairSWAngleTableKokkos*);
};

}

#endif
#endif

