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

FixStyle(rx/kk,FixRxKokkos<Device>)
FixStyle(rx/kk/device,FixRxKokkos<Device>)
FixStyle(rx/kk/host,FixRxKokkos<Host>)

#else

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
  volatile s_CounterType& operator+=(const volatile s_CounterType &rhs) volatile
  {
    nSteps += rhs.nSteps;
    nIters += rhs.nIters;
    nFuncs += rhs.nFuncs;
    nFails += rhs.nFails;
    return *this;
  }
};
typedef struct s_CounterType CounterType;

template <ExecutionSpace Space>
class FixRxKokkos : public FixRX {
 public:
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;
  typedef typename GetFloatType<Space>::type SPACE_FLOAT;

  FixRxKokkos(class LAMMPS *, int, char **);
  virtual ~FixRxKokkos();
  virtual void init();
  void init_list(int, class NeighList *);
  void post_constructor();
  virtual void setup_pre_force(int);
  virtual void pre_force(int);

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
  PairDPDfdtEnergyKokkos<Space>* pairDPDEKK;
  KK_FLOAT VDPD;

  KK_FLOAT boltz;
  KK_FLOAT t_stop;

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
    StridedArrayType<SPACE_FLOAT,1> kFor;
    StridedArrayType<SPACE_FLOAT,1> rxnRateLaw;
  };

  void solve_reactions(const int vflag, const bool isPreForce);

  int rhs       (KK_FLOAT, const KK_FLOAT *, KK_FLOAT *, void *) const;
  int rhs_dense (KK_FLOAT, const KK_FLOAT *, KK_FLOAT *, void *) const;
  int rhs_sparse(KK_FLOAT, const KK_FLOAT *, KK_FLOAT *, void *) const;

  template <typename VectorType, typename UserDataType>
    KOKKOS_INLINE_FUNCTION
  int k_rhs       (KK_FLOAT, const VectorType&, VectorType&, UserDataType& ) const;

  template <typename VectorType, typename UserDataType>
    KOKKOS_INLINE_FUNCTION
  int k_rhs_dense (KK_FLOAT, const VectorType&, VectorType&, UserDataType& ) const;

  template <typename VectorType, typename UserDataType>
    KOKKOS_INLINE_FUNCTION
  int k_rhs_sparse(KK_FLOAT, const VectorType&, VectorType&, UserDataType& ) const;

  //!< Classic Runge-Kutta 4th-order stepper.
  void rk4(const KK_FLOAT t_stop, KK_FLOAT *y, KK_FLOAT *rwork, void *v_params) const;

  //!< Runge-Kutta-Fehlberg ODE Solver.
  void rkf45(const int neq, const KK_FLOAT t_stop, KK_FLOAT *y, KK_FLOAT *rwork, void *v_params, CounterType& counter) const;

  //!< Runge-Kutta-Fehlberg ODE stepper function.
  void rkf45_step (const int neq, const KK_FLOAT h, KK_FLOAT y[], KK_FLOAT y_out[],
                   KK_FLOAT rwk[], void *) const;

  //!< Initial step size estimation for the Runge-Kutta-Fehlberg ODE solver.
  int rkf45_h0 (const int neq, const KK_FLOAT t, const KK_FLOAT t_stop,
                     const KK_FLOAT hmin, const KK_FLOAT hmax,
                     KK_FLOAT& h0, KK_FLOAT y[], KK_FLOAT rwk[], void *v_params) const;

  //!< Classic Runge-Kutta 4th-order stepper.
  template <typename VectorType, typename UserDataType>
    KOKKOS_INLINE_FUNCTION
  void k_rk4(const KK_FLOAT t_stop, VectorType& y, VectorType& rwork, UserDataType& userData) const;

  //!< Runge-Kutta-Fehlberg ODE Solver.
  template <typename VectorType, typename UserDataType>
    KOKKOS_INLINE_FUNCTION
  void k_rkf45(const int neq, const KK_FLOAT t_stop, VectorType& y, VectorType& rwork, UserDataType& userData, CounterType& counter) const;

  //!< Runge-Kutta-Fehlberg ODE stepper function.
  template <typename VectorType, typename UserDataType>
    KOKKOS_INLINE_FUNCTION
  void k_rkf45_step (const int neq, const KK_FLOAT h, VectorType& y, VectorType& y_out,
                     VectorType& rwk, UserDataType& userData) const;

  //!< Initial step size estimation for the Runge-Kutta-Fehlberg ODE solver.
  template <typename VectorType, typename UserDataType>
    KOKKOS_INLINE_FUNCTION
  int k_rkf45_h0 (const int neq, const KK_FLOAT t, const KK_FLOAT t_stop,
                  const KK_FLOAT hmin, const KK_FLOAT hmax,
                  KK_FLOAT& h0, VectorType& y, VectorType& rwk, UserDataType& userData) const;

  //!< ODE Solver diagnostics.
  void odeDiagnostics(void);

  //!< Special counters per-ode.
  int *diagnosticCounterPerODEnSteps;
  int *diagnosticCounterPerODEnFuncs;
  DAT::tdual_int_1d k_diagnosticCounterPerODEnSteps;
  DAT::tdual_int_1d k_diagnosticCounterPerODEnFuncs;
  //typename AT::t_int_1d d_diagnosticCounterPerODEnSteps;
  //typename AT::t_int_1d d_diagnosticCounterPerODEnFuncs;
  typename AT::t_int_1d d_diagnosticCounterPerODEnSteps;
  typename AT::t_int_1d d_diagnosticCounterPerODEnFuncs;
  HAT::t_int_1d h_diagnosticCounterPerODEnSteps;
  HAT::t_int_1d h_diagnosticCounterPerODEnFuncs;

  struct KineticsDual
  {
    // Arrhenius rate coefficients.
    DAT::tdual_float_1d k_Arr, k_nArr, k_Ea;

    // Dense versions.
    DAT::tdual_float_2d k_stoich, k_stoichReactants, k_stoichProducts;

    // Sparse versions.
    DAT::tdual_int_2d   k_nuk, k_inu;
    DAT::tdual_float_2d k_nu;
    DAT::tdual_int_1d   k_isIntegral;
  };

  struct KineticsDevice
  {
    // Arrhenius rate coefficients.
    typename AT::t_float_1d Arr, nArr, Ea;

    // Dense versions.
    typename AT::t_float_2d stoich, stoichReactants, stoichProducts;

    // Sparse versions.
    typename AT::t_int_2d   nuk, inu;
    typename AT::t_float_2d nu;
    typename AT::t_int_1d   isIntegral;

    KineticsDevice() {}

    KineticsDevice(const KineticsDual &rhs) {
      Arr = DualViewHelper<Space>::view(rhs.k_Arr);
      nArr = DualViewHelper<Space>::view(rhs.k_nArr);
      Ea = DualViewHelper<Space>::view(rhs.k_Ea);
      stoich = DualViewHelper<Space>::view(rhs.k_stoich);
      stoichReactants = DualViewHelper<Space>::view(rhs.k_stoichReactants);
      stoichProducts = DualViewHelper<Space>::view(rhs.k_stoichProducts);
      nuk = DualViewHelper<Space>::view(rhs.k_nuk);
      inu = DualViewHelper<Space>::view(rhs.k_inu);
      nu = DualViewHelper<Space>::view(rhs.k_nu);
      isIntegral = DualViewHelper<Space>::view(rhs.k_isIntegral);
    }
  };

  struct KineticsHost
  { 
    
    // Arrhenius rate coefficients.
    HAT::t_float_1d Arr, nArr, Ea;
    
    // Dense versions.
    HAT::t_float_2d stoich, stoichReactants, stoichProducts;
    
    // Sparse versions.
    HAT::t_int_2d   nuk, inu;
    HAT::t_float_2d nu;
    HAT::t_int_1d   isIntegral;

    KineticsHost() {}
    
    KineticsHost(const KineticsDual &rhs) {
      Arr = rhs.k_Arr.h_view;
      nArr = rhs.k_nArr.h_view;
      Ea = rhs.k_Ea.h_view;
      stoich = rhs.k_stoich.h_view;
      stoichReactants = rhs.k_stoichReactants.h_view;
      stoichProducts = rhs.k_stoichProducts.h_view;
      nuk = rhs.k_nuk.h_view;
      inu = rhs.k_inu.h_view;
      nu = rhs.k_nu.h_view;
      isIntegral = rhs.k_isIntegral.h_view;
    }
  };

  //!< Kokkos versions of the kinetics data.
  KineticsDual k_kineticsData;
  KineticsHost h_kineticsData;
  KineticsDevice d_kineticsData;

  bool update_kinetics_data;

  void create_kinetics_data(void);

  // Need a dual-view and device-view for dpdThetaLocal and sumWeights since they're used in several callbacks.
  DAT::tdual_float_1d k_dpdThetaLocal, k_sumWeights;
  //typename AT::t_float_1d d_dpdThetaLocal, d_sumWeights;
  typename AT::t_float_1d d_dpdThetaLocal, d_sumWeights;
  HAT::t_float_1d h_dpdThetaLocal, h_sumWeights;

  typename AT::t_float_1d_3_randomread d_x       ;
  typename AT::t_int_1d_randomread  d_type    ;
  typename AT::t_float_1d          d_dpdTheta;

  DAT::tdual_float_2d k_cutsq;
  typename AT::t_float_2d     d_cutsq;
  //KK_FLOAT **h_cutsq;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d       d_ilist    ;
  typename AT::t_int_1d       d_numneigh ;

  typename AT::t_float_2d  d_dvector;
  typename AT::t_int_1d    d_mask   ;

  typename AT::t_float_1d d_scratchSpace;
  size_t scratchSpaceSize;

  // Error flag for any failures.
  DAT::tdual_int_scalar k_error_flag;

  template <int WT_FLAG, int LOCAL_TEMP_FLAG, bool NEWTON_PAIR, int NEIGHFLAG>
  void computeLocalTemperature();

  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  int pack_forward_comm(int , int *, double *, int, int *);
  void unpack_forward_comm(int , int , double *);

 //private: // replicate a few from FixRX
  int my_restartFlag;
  int nlocal;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
