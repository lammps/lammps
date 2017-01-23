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

namespace LAMMPS_NS {

template <bool _setToZero>
struct TagFixRxKokkosSolver
{
  enum { setToZero = (_setToZero == true) ? 1 : 0 };
};

template <typename DeviceType>
class FixRxKokkos : public FixRX {
 public:
  FixRxKokkos(class LAMMPS *, int, char **);
  virtual ~FixRxKokkos();
  virtual void init();
  virtual void pre_force(int);

  //template <typename SolverTag>
  //  KOKKOS_INLINE_FUNCTION
  //void operator()(SolverTag, const int&) const;

  struct CounterType
  {
    int nSteps, nIters, nFuncs, nFails;

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

 protected:
  PairDPDfdtEnergyKokkos<DeviceType>* pairDPDEKK;
  double VDPD;

  void solve_reactions(void);

  int rhs(double, const double *, double *, void *) const;
  int rhs_dense (double, const double *, double *, void *) const;
  int rhs_sparse(double, const double *, double *, void *) const;

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

};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
