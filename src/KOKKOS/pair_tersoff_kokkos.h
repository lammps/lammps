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
PairStyle(tersoff/kk,PairTersoffKokkos<LMPDeviceType>);
PairStyle(tersoff/kk/device,PairTersoffKokkos<LMPDeviceType>);
PairStyle(tersoff/kk/host,PairTersoffKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_TERSOFF_KOKKOS_H
#define LMP_PAIR_TERSOFF_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_tersoff.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<int NEIGHFLAG, int EVFLAG>
struct TagPairTersoffCompute{};

struct TagPairTersoffComputeShortNeigh{};

template<class DeviceType>
class PairTersoffKokkos : public PairTersoff {
 public:
  enum {EnabledNeighFlags=HALF|HALFTHREAD};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  // Static blocking size for PairTersoffCompute, EVFLAG == 0
  static constexpr int block_size_compute_tersoff_force = 128;
  // EVFLAG == 1, intentionally different due to how Kokkos implements
  // reductions vs simple parallel_for
  static constexpr int block_size_compute_tersoff_energy = 256;

  PairTersoffKokkos(class LAMMPS *);
  ~PairTersoffKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;

  // RangePolicy versions
  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTersoffCompute<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTersoffCompute<NEIGHFLAG,EVFLAG>, const int&) const;

  // MDRangePolicy versions
  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTersoffCompute<NEIGHFLAG,EVFLAG>, const int&, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTersoffCompute<NEIGHFLAG,EVFLAG>, const int&, const int&) const;

  // TeamPolicy versions
  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTersoffCompute<NEIGHFLAG,EVFLAG>, const typename Kokkos::TeamPolicy<DeviceType, TagPairTersoffCompute<NEIGHFLAG,EVFLAG> >::member_type&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTersoffCompute<NEIGHFLAG,EVFLAG>, const typename Kokkos::TeamPolicy<DeviceType, TagPairTersoffCompute<NEIGHFLAG,EVFLAG> >::member_type&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTersoffComputeShortNeigh, const int&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void tersoff_compute(const int&, EV_FLOAT&) const;

  KOKKOS_INLINE_FUNCTION
  double ters_fc_k(const Param& param, const F_FLOAT &r) const;

  KOKKOS_INLINE_FUNCTION
  double ters_dfc(const Param& param, const F_FLOAT &r) const;

  KOKKOS_INLINE_FUNCTION
  void ters_fc_k_and_ters_dfc(const Param& param, const F_FLOAT &r, double &fc, double &dfc) const;

  KOKKOS_INLINE_FUNCTION
  double ters_fa_k(const Param& param, const F_FLOAT &r) const;

  KOKKOS_INLINE_FUNCTION
  double ters_dfa(const Param& param, const F_FLOAT &r) const;

  KOKKOS_INLINE_FUNCTION
  void ters_fa_k_and_ters_dfa(const Param& param, const F_FLOAT &r, double &fa, double &dfa) const;

  KOKKOS_INLINE_FUNCTION
  double ters_bij_k(const Param& param, const F_FLOAT &bo) const;

  KOKKOS_INLINE_FUNCTION
  double ters_dbij(const Param& param, const F_FLOAT &bo) const;

  KOKKOS_INLINE_FUNCTION
  void ters_bij_k_and_ters_dbij(const Param& param, const F_FLOAT &bo, double &bij, double &prefactor) const;

  KOKKOS_INLINE_FUNCTION
  double bondorder(const Param& param,
              const F_FLOAT &rij, const F_FLOAT &dx1, const F_FLOAT &dy1, const F_FLOAT &dz1,
              const F_FLOAT &rik, const F_FLOAT &dx2, const F_FLOAT &dy2, const F_FLOAT &dz2) const;

  KOKKOS_INLINE_FUNCTION
  double ters_gijk(const Param& param, const F_FLOAT &cos) const;

  KOKKOS_INLINE_FUNCTION
  double ters_dgijk(const Param& param, const F_FLOAT &cos) const;

  KOKKOS_INLINE_FUNCTION
  void ters_gijk_and_ters_dgijk(const Param& param, const F_FLOAT &cos, double& gijk, double& dgijk) const;

  KOKKOS_INLINE_FUNCTION
  void ters_dthb(const Param& param, const F_FLOAT &prefactor,
              const F_FLOAT &rij, const F_FLOAT &dx1, const F_FLOAT &dy1, const F_FLOAT &dz1,
              const F_FLOAT &rik, const F_FLOAT &dx2, const F_FLOAT &dy2, const F_FLOAT &dz2,
              F_FLOAT *fi, F_FLOAT *fj, F_FLOAT *fk) const;

  KOKKOS_INLINE_FUNCTION
  void ters_dthbj(const Param& param, const F_FLOAT &prefactor,
              const F_FLOAT &rij, const F_FLOAT &dx1, const F_FLOAT &dy1, const F_FLOAT &dz1,
              const F_FLOAT &rik, const F_FLOAT &dx2, const F_FLOAT &dy2, const F_FLOAT &dz2,
              F_FLOAT *fj, F_FLOAT *fk) const;

  KOKKOS_INLINE_FUNCTION
  void ters_dthbk(const Param& param, const F_FLOAT &prefactor,
              const F_FLOAT &rij, const F_FLOAT &dx1, const F_FLOAT &dy1, const F_FLOAT &dz1,
              const F_FLOAT &rik, const F_FLOAT &dx2, const F_FLOAT &dy2, const F_FLOAT &dz2,
              F_FLOAT *fk) const;

  KOKKOS_INLINE_FUNCTION
  double vec3_dot(const F_FLOAT x[3], const double y[3]) const {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
  }

  KOKKOS_INLINE_FUNCTION
  void vec3_add(const F_FLOAT x[3], const double y[3], double * const z) const {
    z[0] = x[0]+y[0]; z[1] = x[1]+y[1]; z[2] = x[2]+y[2];
  }

  KOKKOS_INLINE_FUNCTION
  void vec3_scale(const F_FLOAT k, const double x[3], double y[3]) const {
    y[0] = k*x[0]; y[1] = k*x[1]; y[2] = k*x[2];
  }

  KOKKOS_INLINE_FUNCTION
  void vec3_scaleadd(const F_FLOAT k, const double x[3], const double y[3], double * const z) const {
    z[0] = k*x[0]+y[0]; z[1] = k*x[1]+y[1]; z[2] = k*x[2]+y[2];
  }

  KOKKOS_INLINE_FUNCTION
  int sbmask(const int& j) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &epair, const F_FLOAT &fpair, const F_FLOAT &delx,
                  const F_FLOAT &dely, const F_FLOAT &delz) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void v_tally3(EV_FLOAT &ev, const int &i, const int &j, const int &k,
                F_FLOAT *fj, F_FLOAT *fk, F_FLOAT *drij, F_FLOAT *drik) const;

  KOKKOS_INLINE_FUNCTION
  void v_tally3_atom(EV_FLOAT &ev, const int &i, const int &j, const int &k,
                F_FLOAT *fj, F_FLOAT *fk, F_FLOAT *drji, F_FLOAT *drjk) const;

  void setup_params() override;

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

  int inum;
  typename AT::t_x_array_randomread x;
  typename AT::t_f_array f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_tagint_1d tag;

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

  typedef Kokkos::DualView<F_FLOAT**[7],Kokkos::LayoutRight,DeviceType> tdual_ffloat_2d_n7;
  typedef typename tdual_ffloat_2d_n7::t_dev_const_randomread t_ffloat_2d_n7_randomread;
  typedef typename tdual_ffloat_2d_n7::t_host t_host_ffloat_2d_n7;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;
  //NeighListKokkos<DeviceType> k_list;

  int neighflag,newton_pair;
  int nlocal,nall,eflag,vflag;

  Kokkos::View<int**,DeviceType> d_neighbors_short;
  Kokkos::View<int*,DeviceType> d_numneigh_short;

  friend void pair_virial_fdotr_compute<PairTersoffKokkos>(PairTersoffKokkos*);
};

}

#endif
#endif

