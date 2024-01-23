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
PairStyle(sw/kk,PairSWKokkos<LMPDeviceType>);
PairStyle(sw/kk/device,PairSWKokkos<LMPDeviceType>);
PairStyle(sw/kk/host,PairSWKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_SW_KOKKOS_H
#define LMP_PAIR_SW_KOKKOS_H

#include "pair_sw.h"
#include "pair_kokkos.h"

template<int NEIGHFLAG, int EVFLAG>
struct TagPairSWCompute{};

struct TagPairSWComputeShortNeigh{};

namespace LAMMPS_NS {

template<class DeviceType>
class PairSWKokkos : public PairSW {
 public:
  enum {EnabledNeighFlags=HALF|HALFTHREAD};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairSWKokkos(class LAMMPS *);
  ~PairSWKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairSWCompute<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairSWCompute<NEIGHFLAG,EVFLAG>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairSWComputeShortNeigh, const int&) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &epair, const F_FLOAT &fpair, const F_FLOAT &delx,
                  const F_FLOAT &dely, const F_FLOAT &delz) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void ev_tally3(EV_FLOAT &ev, const int &i, const int &j, int &k,
            const F_FLOAT &evdwl, const F_FLOAT &ecoul,
                       F_FLOAT *fj, F_FLOAT *fk, F_FLOAT *drji, F_FLOAT *drki) const;

  KOKKOS_INLINE_FUNCTION
  void ev_tally3_atom(EV_FLOAT &ev, const int &i,
            const F_FLOAT &evdwl, const F_FLOAT &ecoul,
                       F_FLOAT *fj, F_FLOAT *fk, F_FLOAT *drji, F_FLOAT *drki) const;

 protected:
  typedef Kokkos::DualView<int***,DeviceType> tdual_int_3d;
  typedef typename tdual_int_3d::t_dev_const_randomread t_int_3d_randomread;
  typedef typename tdual_int_3d::t_host t_host_int_3d;

  t_int_3d_randomread d_elem3param;
  typename AT::t_int_1d_randomread d_map;

  typedef Kokkos::DualView<Param*,DeviceType> tdual_param_1d;
  typedef typename tdual_param_1d::t_dev t_param_1d;
  typedef typename tdual_param_1d::t_host t_host_param_1d;

  t_param_1d d_params;

  void setup_params() override;

  KOKKOS_INLINE_FUNCTION
  void twobody(const Param&, const F_FLOAT&, F_FLOAT&, const int&, F_FLOAT&) const;

  KOKKOS_INLINE_FUNCTION
  void threebody_kk(const Param&, const Param&, const Param&, const F_FLOAT&, const F_FLOAT&, F_FLOAT *, F_FLOAT *,
                    F_FLOAT *, F_FLOAT *, const int&, F_FLOAT&) const;

  typename AT::t_x_array_randomread x;
  typename AT::t_f_array f;
  typename AT::t_tagint_1d tag;
  typename AT::t_int_1d_randomread type;

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

  DupScatterView<F_FLOAT*[3], typename DAT::t_f_array::array_layout> dup_f;
  DupScatterView<E_FLOAT*, typename DAT::t_efloat_1d::array_layout> dup_eatom;
  DupScatterView<F_FLOAT*[6], typename DAT::t_virial_array::array_layout> dup_vatom;

  NonDupScatterView<F_FLOAT*[3], typename DAT::t_f_array::array_layout> ndup_f;
  NonDupScatterView<E_FLOAT*, typename DAT::t_efloat_1d::array_layout> ndup_eatom;
  NonDupScatterView<F_FLOAT*[6], typename DAT::t_virial_array::array_layout> ndup_vatom;

  typename AT::t_int_1d_randomread d_type2frho;
  typename AT::t_int_2d_randomread d_type2rhor;
  typename AT::t_int_2d_randomread d_type2z2r;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;
  //NeighListKokkos<DeviceType> k_list;

  int neighflag,newton_pair;
  int nlocal,nall,eflag,vflag;

  int inum;
  Kokkos::View<int**,DeviceType> d_neighbors_short;
  Kokkos::View<int*,DeviceType> d_numneigh_short;


  friend void pair_virial_fdotr_compute<PairSWKokkos>(PairSWKokkos*);
};

}

#endif
#endif

