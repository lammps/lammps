/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(rx/kk,FixRxKokkos<LMPDeviceType>);
FixStyle(rx/kk/device,FixRxKokkos<LMPDeviceType>);
FixStyle(rx/kk/host,FixRxKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_RX_KOKKOS_H
#define LMP_FIX_RX_KOKKOS_H

#include "fix_rx.h"
#include "pair_dpd_fdt_energy_kokkos.h"
#include "kokkos_type.h"
#include "neigh_list.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

struct Tag_FixRxKokkos_zeroTemperatureViews {};
struct Tag_FixRxKokkos_zeroCounterViews {};

template <int WT_FLAG, bool NEWTON_PAIR, int NEIGHFLAG>
struct Tag_FixRxKokkos_firstPairOperator {};

template <int WT_FLAG, int LOCAL_TEMP_FLAG>
struct Tag_FixRxKokkos_2ndPairOperator {};

template <bool ZERO_RATES>
struct Tag_FixRxKokkos_solveSystems {};

struct s_CounterType
{
  int nSteps, nIters, nFuncs, nFails;

  KOKKOS_INLINE_FUNCTION
  s_CounterType() : nSteps(0), nIters(0), nFuncs(0), nFails(0) {};

  KOKKOS_INLINE_FUNCTION
  s_CounterType& operator+=(const s_CounterType &rhs)
  {
    nSteps += rhs.nSteps;
    nIters += rhs.nIters;
    nFuncs += rhs.nFuncs;
    nFails += rhs.nFails;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const volatile s_CounterType &rhs) volatile
  {
    nSteps += rhs.nSteps;
    nIters += rhs.nIters;
    nFuncs += rhs.nFuncs;
    nFails += rhs.nFails;
  }
};
typedef struct s_CounterType CounterType;

template <class DeviceType>
class FixRxKokkos : public FixRX {
 public:
  typedef ArrayTypes<DeviceType> AT;

  FixRxKokkos(class LAMMPS *, int, char **);
  ~FixRxKokkos() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void post_constructor() override;
  void setup_pre_force(int) override;
  void pre_force(int) override;

  // Define a value_type here for the reduction operator on CounterType.
  typedef CounterType value_type;

  KOKKOS_INLINE_FUNCTION
  void operator()(Tag_FixRxKokkos_zeroCounterViews, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(Tag_FixRxKokkos_zeroTemperatureViews, const int&) const;

  template <int WT_FLAG, bool NEWTON_PAIR, int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(Tag_FixRxKokkos_firstPairOperator<WT_FLAG,NEWTON_PAIR,NEIGHFLAG>, const int&) const;

  template <int WT_FLAG, int LOCAL_TEMP_FLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(Tag_FixRxKokkos_2ndPairOperator<WT_FLAG,LOCAL_TEMP_FLAG>, const int&) const;

  template <bool ZERO_RATES>
  KOKKOS_INLINE_FUNCTION
  void operator()(Tag_FixRxKokkos_solveSystems<ZERO_RATES>, const int&, CounterType&) const;

 //protected:
  PairDPDfdtEnergyKokkos<DeviceType>* pairDPDEKK;
  double VDPD;

  double boltz;
  double t_stop;

  template <typename T, int stride = 1>
  struct StridedArrayType
  {
    typedef T value_type;
    enum { Stride = stride };

    value_type *m_data;

    KOKKOS_INLINE_FUNCTION
    StridedArrayType() : m_data(nullptr) {}
    KOKKOS_INLINE_FUNCTION
    StridedArrayType(value_type *ptr) : m_data(ptr) {}

    KOKKOS_INLINE_FUNCTION       value_type& operator()(const int idx)       { return m_data[Stride*idx]; }
    KOKKOS_INLINE_FUNCTION const value_type& operator()(const int idx) const { return m_data[Stride*idx]; }
    KOKKOS_INLINE_FUNCTION       value_type& operator[](const int idx)       { return m_data[Stride*idx]; }
    KOKKOS_INLINE_FUNCTION const value_type& operator[](const int idx) const { return m_data[Stride*idx]; }
  };

  template <int stride = 1>
  struct UserRHSDataKokkos
  {
    StridedArrayType<double,1> kFor;
    StridedArrayType<double,1> rxnRateLaw;
  };

  void solve_reactions(const int vflag, const bool isPreForce);

  int rhs       (double, const double *, double *, void *) const;
  int rhs_dense (double, const double *, double *, void *) const;
  int rhs_sparse(double, const double *, double *, void *) const;

  template <typename VectorType, typename UserDataType>
    KOKKOS_INLINE_FUNCTION
  int k_rhs       (double, const VectorType&, VectorType&, UserDataType& ) const;

  template <typename VectorType, typename UserDataType>
    KOKKOS_INLINE_FUNCTION
  int k_rhs_dense (double, const VectorType&, VectorType&, UserDataType& ) const;

  template <typename VectorType, typename UserDataType>
    KOKKOS_INLINE_FUNCTION
  int k_rhs_sparse(double, const VectorType&, VectorType&, UserDataType& ) const;

  //!< Classic Runge-Kutta 4th-order stepper.
  void rk4(const double t_stop, double *y, double *rwork, void *v_params) const;

  //!< Runge-Kutta-Fehlberg ODE Solver.
  void rkf45(const int neq, const double t_stop, double *y, double *rwork, void *v_params, CounterType& counter) const;

  //!< Runge-Kutta-Fehlberg ODE stepper function.
  void rkf45_step (const int neq, const double h, double y[], double y_out[],
                   double rwk[], void *) const;

  //!< Initial step size estimation for the Runge-Kutta-Fehlberg ODE solver.
  int rkf45_h0 (const int neq, const double t, const double t_stop,
                     const double hmin, const double hmax,
                     double& h0, double y[], double rwk[], void *v_params) const;

  //!< Classic Runge-Kutta 4th-order stepper.
  template <typename VectorType, typename UserDataType>
    KOKKOS_INLINE_FUNCTION
  void k_rk4(const double t_stop, VectorType& y, VectorType& rwork, UserDataType& userData) const;

  //!< Runge-Kutta-Fehlberg ODE Solver.
  template <typename VectorType, typename UserDataType>
    KOKKOS_INLINE_FUNCTION
  void k_rkf45(const int neq, const double t_stop, VectorType& y, VectorType& rwork, UserDataType& userData, CounterType& counter) const;

  //!< Runge-Kutta-Fehlberg ODE stepper function.
  template <typename VectorType, typename UserDataType>
    KOKKOS_INLINE_FUNCTION
  void k_rkf45_step (const int neq, const double h, VectorType& y, VectorType& y_out,
                     VectorType& rwk, UserDataType& userData) const;

  //!< Initial step size estimation for the Runge-Kutta-Fehlberg ODE solver.
  template <typename VectorType, typename UserDataType>
    KOKKOS_INLINE_FUNCTION
  int k_rkf45_h0 (const int neq, const double t, const double t_stop,
                  const double hmin, const double hmax,
                  double& h0, VectorType& y, VectorType& rwk, UserDataType& userData) const;

  //!< ODE Solver diagnostics.
  void odeDiagnostics();

  //!< Special counters per-ode.
  int *diagnosticCounterPerODEnSteps;
  int *diagnosticCounterPerODEnFuncs;
  DAT::tdual_int_1d k_diagnosticCounterPerODEnSteps;
  DAT::tdual_int_1d k_diagnosticCounterPerODEnFuncs;
  //typename ArrayTypes<DeviceType>::t_int_1d d_diagnosticCounterPerODEnSteps;
  //typename ArrayTypes<DeviceType>::t_int_1d d_diagnosticCounterPerODEnFuncs;
  typename AT::t_int_1d d_diagnosticCounterPerODEnSteps;
  typename AT::t_int_1d d_diagnosticCounterPerODEnFuncs;
  HAT::t_int_1d h_diagnosticCounterPerODEnSteps;
  HAT::t_int_1d h_diagnosticCounterPerODEnFuncs;

  template <typename KokkosDeviceType>
  struct KineticsType
  {
    // Arrhenius rate coefficients.
    typename ArrayTypes<KokkosDeviceType>::t_float_1d Arr, nArr, Ea;

    // Dense versions.
    typename ArrayTypes<KokkosDeviceType>::t_float_2d stoich, stoichReactants, stoichProducts;

    // Sparse versions.
    typename ArrayTypes<KokkosDeviceType>::t_int_2d   nuk, inu;
    typename ArrayTypes<KokkosDeviceType>::t_float_2d nu;
    typename ArrayTypes<KokkosDeviceType>::t_int_1d   isIntegral;
  };

  //!< Kokkos versions of the kinetics data.
  KineticsType<LMPHostType> h_kineticsData;
  KineticsType<DeviceType>  d_kineticsData;

  bool update_kinetics_data;

  void create_kinetics_data();

  // Need a dual-view and device-view for dpdThetaLocal and sumWeights since they're used in several callbacks.
  DAT::tdual_efloat_1d k_dpdThetaLocal, k_sumWeights;
  //typename ArrayTypes<DeviceType>::t_efloat_1d d_dpdThetaLocal, d_sumWeights;
  typename AT::t_efloat_1d d_dpdThetaLocal, d_sumWeights;
  HAT::t_efloat_1d h_dpdThetaLocal, h_sumWeights;

  typename AT::t_x_array_randomread d_x;
  typename AT::t_int_1d_randomread  d_type;
  typename AT::t_efloat_1d          d_dpdTheta;

  typename AT::tdual_ffloat_2d k_cutsq;
  typename AT::t_ffloat_2d     d_cutsq;
  //double **h_cutsq;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d       d_ilist;
  typename AT::t_int_1d       d_numneigh;

  typename AT::t_float_2d  d_dvector;
  typename AT::t_int_1d    d_mask;

  typename AT::t_double_1d d_scratchSpace;
  size_t scratchSpaceSize;

  // Error flag for any failures.
  DAT::tdual_int_scalar k_error_flag;

  template <int WT_FLAG, int LOCAL_TEMP_FLAG, bool NEWTON_PAIR, int NEIGHFLAG>
  void computeLocalTemperature();

  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  int pack_forward_comm(int , int *, double *, int, int *) override;
  void unpack_forward_comm(int , int , double *) override;

 //private: // replicate a few from FixRX
  int my_restartFlag;
  int nlocal;
};

}

#endif
#endif

