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

#ifdef FIX_CLASS
// clang-format off
FixStyle(langevin/kk,FixLangevinKokkos<LMPDeviceType>);
FixStyle(langevin/kk/device,FixLangevinKokkos<LMPDeviceType>);
FixStyle(langevin/kk/host,FixLangevinKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_LANGEVIN_KOKKOS_H
#define LMP_FIX_LANGEVIN_KOKKOS_H

#include "fix_langevin.h"
#include "kokkos_type.h"
#include "Kokkos_Random.hpp"
#include "comm_kokkos.h"

namespace LAMMPS_NS {

  struct s_FSUM {
    double fx, fy, fz;
    KOKKOS_INLINE_FUNCTION
    s_FSUM() {
      fx = fy = fz = 0.0;
    }
    KOKKOS_INLINE_FUNCTION
    s_FSUM& operator+=(const s_FSUM &rhs) {
      fx += rhs.fx;
      fy += rhs.fy;
      fz += rhs.fz;
      return *this;
    }
  };
  typedef s_FSUM FSUM;

  template<class DeviceType>
  class FixLangevinKokkos;

  template <class DeviceType>
  struct FixLangevinKokkosInitialIntegrateFunctor;

  template<class DeviceType,int Tp_TSTYLEATOM, int Tp_GJF, int Tp_TALLY,
    int Tp_BIAS, int Tp_RMASS, int Tp_ZERO>
  struct FixLangevinKokkosPostForceFunctor;

  template<class DeviceType> struct FixLangevinKokkosZeroForceFunctor;

  template<class DeviceType> struct FixLangevinKokkosTallyEnergyFunctor;

  template<class DeviceType>
  class FixLangevinKokkos : public FixLangevin {
   public:
    FixLangevinKokkos(class LAMMPS *, int, char **);
    ~FixLangevinKokkos() override;

    void cleanup_copy();
    void init() override;
    void initial_integrate(int) override;
    void post_force(int) override;
    void reset_dt() override;
    void grow_arrays(int) override;
    void copy_arrays(int i, int j, int delflag) override;
    double compute_scalar() override;
    void end_of_step() override;

    KOKKOS_INLINE_FUNCTION
      void initial_integrate_item(int) const;

    KOKKOS_INLINE_FUNCTION
      void initial_integrate_rmass_item(int) const;

    template<int Tp_TSTYLEATOM, int Tp_GJF, int Tp_TALLY,
      int Tp_BIAS, int Tp_RMASS, int Tp_ZERO>
      KOKKOS_INLINE_FUNCTION
      FSUM post_force_item(int) const;

    KOKKOS_INLINE_FUNCTION
      void zero_force_item(int) const;

    KOKKOS_INLINE_FUNCTION
      double compute_energy_item(int) const;

    KOKKOS_INLINE_FUNCTION
      void end_of_step_item(int) const;

    KOKKOS_INLINE_FUNCTION
      void end_of_step_rmass_item(int) const;

  private:
    class CommKokkos *commKK;

    typename ArrayTypes<DeviceType>::t_float_1d rmass;
    typename ArrayTypes<DeviceType>::t_float_1d mass;
    typename ArrayTypes<DeviceType>::tdual_double_2d k_franprev;
    typename ArrayTypes<DeviceType>::t_double_2d d_franprev;
    HAT::t_double_2d h_franprev;

    typename ArrayTypes<DeviceType>::tdual_double_2d k_lv;
    typename ArrayTypes<DeviceType>::t_double_2d d_lv;
    HAT::t_double_2d h_lv;

    typename ArrayTypes<DeviceType>::tdual_double_2d k_flangevin;
    typename ArrayTypes<DeviceType>::t_double_2d d_flangevin;
    HAT::t_double_2d h_flangevin;

    typename ArrayTypes<DeviceType>::tdual_double_1d k_tforce;
    typename ArrayTypes<DeviceType>::t_double_1d d_tforce;
    HAT::t_double_1d h_tforce;

    typename ArrayTypes<DeviceType>::t_v_array v;
    typename ArrayTypes<DeviceType>::t_f_array f;
    typename ArrayTypes<DeviceType>::t_int_1d type;
    typename ArrayTypes<DeviceType>::t_int_1d mask;

    typename ArrayTypes<DeviceType>::tdual_double_1d k_gfactor1, k_gfactor2, k_ratio;
    typename ArrayTypes<DeviceType>::t_double_1d d_gfactor1, d_gfactor2, d_ratio;
    HAT::t_double_1d h_gfactor1, h_gfactor2, h_ratio;

    typedef Kokkos::DualView<double[3], DeviceType>
      tdual_double_1d_3n;
    tdual_double_1d_3n k_fsumall;
    typename tdual_double_1d_3n::t_dev d_fsumall;
    typename tdual_double_1d_3n::t_host h_fsumall;

    double boltz,dt,mvv2e,ftm2v,fran_prop_const;

    void compute_target();

    Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
    typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;

  };

  template <class DeviceType>
  struct FixLangevinKokkosInitialIntegrateFunctor  {
    typedef DeviceType  device_type ;
    FixLangevinKokkos<DeviceType> c;

  FixLangevinKokkosInitialIntegrateFunctor(FixLangevinKokkos<DeviceType>* c_ptr):
    c(*c_ptr) {c.cleanup_copy();};

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i) const {
      c.initial_integrate_item(i);
    }
  };


  template <class DeviceType,int Tp_TSTYLEATOM, int Tp_GJF, int Tp_TALLY,
    int Tp_BIAS, int Tp_RMASS, int Tp_ZERO>
    struct FixLangevinKokkosPostForceFunctor {
      typedef DeviceType  device_type;
      typedef FSUM value_type;
      FixLangevinKokkos<DeviceType> c;

    FixLangevinKokkosPostForceFunctor(FixLangevinKokkos<DeviceType>* c_ptr):
      c(*c_ptr) {}
      ~FixLangevinKokkosPostForceFunctor() {c.cleanup_copy();}

      KOKKOS_INLINE_FUNCTION
      void operator()(const int i) const {
        c.template post_force_item<Tp_TSTYLEATOM,Tp_GJF, Tp_TALLY,
          Tp_BIAS,Tp_RMASS,Tp_ZERO>(i);
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const int i, value_type &fsum) const {

        fsum += c.template post_force_item<Tp_TSTYLEATOM,Tp_GJF, Tp_TALLY,
          Tp_BIAS,Tp_RMASS,Tp_ZERO>(i);
      }

      KOKKOS_INLINE_FUNCTION
      static void init(value_type &update) {
        update.fx = 0.0;
        update.fy = 0.0;
        update.fz = 0.0;
      }
      KOKKOS_INLINE_FUNCTION
      static void join(value_type &update,
                       const value_type &source) {
        update.fx += source.fx;
        update.fy += source.fy;
        update.fz += source.fz;
      }
    };

  template <class DeviceType>
    struct FixLangevinKokkosZeroForceFunctor {
      typedef DeviceType  device_type ;
      FixLangevinKokkos<DeviceType> c;

    FixLangevinKokkosZeroForceFunctor(FixLangevinKokkos<DeviceType>* c_ptr):
      c(*c_ptr) {c.cleanup_copy();}

      KOKKOS_INLINE_FUNCTION
      void operator()(const int i) const {
        c.zero_force_item(i);
      }
    };

  template<class DeviceType>
    struct FixLangevinKokkosTallyEnergyFunctor {
      typedef DeviceType  device_type ;
      FixLangevinKokkos<DeviceType> c;
      typedef double value_type;
    FixLangevinKokkosTallyEnergyFunctor(FixLangevinKokkos<DeviceType>* c_ptr):
      c(*c_ptr) {c.cleanup_copy();}

      KOKKOS_INLINE_FUNCTION
      void operator()(const int i, value_type &energy) const {
        energy += c.compute_energy_item(i);
      }
      KOKKOS_INLINE_FUNCTION
      static void init(value_type &update) {
        update = 0.0;
      }
      KOKKOS_INLINE_FUNCTION
      static void join(value_type &update,
                       const value_type &source) {
        update += source;
      }
    };

  template <class DeviceType, int RMass>
  struct FixLangevinKokkosEndOfStepFunctor {
    typedef DeviceType  device_type ;
    FixLangevinKokkos<DeviceType> c;

    FixLangevinKokkosEndOfStepFunctor(FixLangevinKokkos<DeviceType>* c_ptr):
      c(*c_ptr) {c.cleanup_copy();}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i) const {
      if (RMass) c.end_of_step_rmass_item(i);
      else c.end_of_step_item(i);
    }
  };
}

#endif
#endif

