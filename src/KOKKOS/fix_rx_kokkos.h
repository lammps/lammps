/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(rx/kk,FixRxKokkos<LMPDeviceType>)
FixStyle(rx/kk/device,FixRxKokkos<LMPDeviceType>)
FixStyle(rx/kk/host,FixRxKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_RX_KOKKOS_H
#define LMP_FIX_RX_KOKKOS_H

#include "fix_rx.h"
#include "pair_dpd_fdt_energy_kokkos.h"
#include "kokkos_type.h"
#include "neigh_list.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template <typename DeviceType>
class FixRxKokkos : public FixRX {
 public:
  FixRxKokkos(class LAMMPS *, int, char **);
  virtual ~FixRxKokkos();
  virtual void init();
  void init_list(int, class NeighList *);
  void post_constructor();
  virtual void setup_pre_force(int);
  virtual void pre_force(int);

  //template <typename SolverTag>
  //  KOKKOS_INLINE_FUNCTION
  //void operator()(SolverTag, const int&) const;

  struct CounterType
  {
    int nSteps, nIters, nFuncs, nFails;

    KOKKOS_INLINE_FUNCTION
    CounterType() : nSteps(0), nIters(0), nFuncs(0), nFails(0) {};

    KOKKOS_INLINE_FUNCTION
    CounterType& operator+=(const CounterType &rhs)
    {
      nSteps += rhs.nSteps;
      nIters += rhs.nIters;
      nFuncs += rhs.nFuncs;
      nFails += rhs.nFails;
      return *this;
    }

    KOKKOS_INLINE_FUNCTION
    volatile CounterType& operator+=(const volatile CounterType &rhs) volatile
    {
      nSteps += rhs.nSteps;
      nIters += rhs.nIters;
      nFuncs += rhs.nFuncs;
      nFails += rhs.nFails;
      return *this;
    }
  };

 //protected:
  PairDPDfdtEnergyKokkos<DeviceType>* pairDPDEKK;
  double VDPD;

  template <typename T, int stride = 1>
  struct StridedArrayType
  {
    typedef T value_type;
    enum { Stride = stride };

    value_type *m_data;

    KOKKOS_INLINE_FUNCTION
    StridedArrayType() : m_data(NULL) {}
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
  void odeDiagnostics(void);

  //!< Special counters per-ode.
  int *diagnosticCounterPerODEnSteps;
  int *diagnosticCounterPerODEnFuncs;
  DAT::tdual_int_1d k_diagnosticCounterPerODEnSteps;
  DAT::tdual_int_1d k_diagnosticCounterPerODEnFuncs;
  //typename ArrayTypes<DeviceType>::t_int_1d d_diagnosticCounterPerODEnSteps;
  //typename ArrayTypes<DeviceType>::t_int_1d d_diagnosticCounterPerODEnFuncs;
  typename DAT::t_int_1d d_diagnosticCounterPerODEnSteps;
  typename DAT::t_int_1d d_diagnosticCounterPerODEnFuncs;
  typename HAT::t_int_1d h_diagnosticCounterPerODEnSteps;
  typename HAT::t_int_1d h_diagnosticCounterPerODEnFuncs;

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

  void create_kinetics_data(void);

  // Need a dual-view and device-view for dpdThetaLocal and sumWeights since they're used in several callbacks.
  DAT::tdual_efloat_1d k_dpdThetaLocal, k_sumWeights;
  //typename ArrayTypes<DeviceType>::t_efloat_1d d_dpdThetaLocal, d_sumWeights;
  typename DAT::t_efloat_1d d_dpdThetaLocal, d_sumWeights;
  typename HAT::t_efloat_1d h_dpdThetaLocal, h_sumWeights;

  template <int WT_FLAG, int LOCAL_TEMP_FLAG, bool NEWTON_PAIR, int NEIGHFLAG>
  void computeLocalTemperature();

  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  int pack_forward_comm(int , int *, double *, int, int *);
  void unpack_forward_comm(int , int , double *);

 //private: // replicate a few from FixRX
  int my_restartFlag;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
