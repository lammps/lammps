/* -*- c++ -*- ----------------------------------------------------------
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

FixStyle(langevin/kk,FixLangevinKokkos<Device>)
FixStyle(langevin/kk/device,FixLangevinKokkos<Device>)
FixStyle(langevin/kk/host,FixLangevinKokkos<Host>)

#else

#ifndef LMP_FIX_LANGEVIN_KOKKOS_H
#define LMP_FIX_LANGEVIN_KOKKOS_H

#include "fix_langevin.h"
#include "kokkos_type.h"
#include "Kokkos_Random.hpp"
#include "comm_kokkos.h"

namespace LAMMPS_NS {

  struct s_FSUM {
    KK_FLOAT fx, fy, fz;
    KOKKOS_INLINE_FUNCTION
    s_FSUM() {
      fx = fy = fz = 0.0;
    }
    KOKKOS_INLINE_FUNCTION
    s_FSUM& operator+=(const s_FSUM &rhs){
      fx += rhs.fx;
      fy += rhs.fy;
      fz += rhs.fz;
      return *this;
    }

    KOKKOS_INLINE_FUNCTION
    volatile s_FSUM& operator+=(const volatile s_FSUM &rhs) volatile {
      fx += rhs.fx;
      fy += rhs.fy;
      fz += rhs.fz;
      return *this;
    }
  };
  typedef s_FSUM FSUM;

  template<ExecutionSpace Space>
    class FixLangevinKokkos;

  template <ExecutionSpace Space>
  class FixLangevinKokkosInitialIntegrateFunctor;

  template<ExecutionSpace Space,int Tp_TSTYLEATOM, int Tp_GJF, int Tp_TALLY,
    int Tp_BIAS, int Tp_RMASS, int Tp_ZERO>
    class FixLangevinKokkosPostForceFunctor;

  template<ExecutionSpace Space> class FixLangevinKokkosZeroForceFunctor;

  template<ExecutionSpace Space> class FixLangevinKokkosTallyEnergyFunctor;

  template<ExecutionSpace Space>
    class FixLangevinKokkos : public FixLangevin {
  public:
    typedef typename GetDeviceType<Space>::value DeviceType;
    typedef ArrayTypes<Space> AT;

    FixLangevinKokkos(class LAMMPS *, int, char **);
    ~FixLangevinKokkos();

    void cleanup_copy();
    void init();
    void initial_integrate(int);
    void post_force(int);
    void reset_dt();
    void grow_arrays(int);
    void copy_arrays(int i, int j, int delflag);
    double compute_scalar();
    void end_of_step();

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
      KK_FLOAT compute_energy_item(int) const;

    KOKKOS_INLINE_FUNCTION
      void end_of_step_item(int) const;

    KOKKOS_INLINE_FUNCTION
      void end_of_step_rmass_item(int) const;

  private:
    class CommKokkos *commKK;

    typename AT::t_float_1d rmass;
    typename AT::t_float_1d mass;
    DAT::tdual_float_2d k_franprev;
    typename AT::t_float_2d d_franprev;
    HAT::t_float_2d h_franprev;

    DAT::tdual_float_2d k_lv;
    typename AT::t_float_2d d_lv;
    HAT::t_float_2d h_lv;

    DAT::tdual_float_2d k_flangevin;
    typename AT::t_float_2d d_flangevin;
    HAT::t_float_2d h_flangevin;

    DAT::tdual_float_1d k_tforce;
    typename AT::t_float_1d d_tforce;
    HAT::t_float_1d h_tforce;

    typename AT::t_float_1d_3 v;
    typename AT::t_float_1d_3 f;
    typename AT::t_int_1d type;
    typename AT::t_int_1d mask;

    DAT::tdual_float_1d k_gfactor1, k_gfactor2, k_ratio;
    typename AT::t_float_1d d_gfactor1, d_gfactor2, d_ratio;
    HAT::t_float_1d h_gfactor1, h_gfactor2, h_ratio;

    typedef Kokkos::DualView<KK_FLOAT[3], DeviceType>
      tdual_float_1d_3n;
    tdual_float_1d_3n k_fsumall;
    typename tdual_float_1d_3n::t_dev d_fsumall;
    typename tdual_float_1d_3n::t_host h_fsumall;

    KK_FLOAT boltz,dt,mvv2e,ftm2v,fran_prop_const;

    void compute_target();

    Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
    typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;

  };

  template <ExecutionSpace Space>
  struct FixLangevinKokkosInitialIntegrateFunctor  {
    typedef typename GetDeviceType<Space>::value DeviceType;
    FixLangevinKokkos<Space> c;

  FixLangevinKokkosInitialIntegrateFunctor(FixLangevinKokkos<Space>* c_ptr):
    c(*c_ptr) {c.cleanup_copy();};

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i) const {
      c.initial_integrate_item(i);
    }
  };


  template <ExecutionSpace Space,int Tp_TSTYLEATOM, int Tp_GJF, int Tp_TALLY,
    int Tp_BIAS, int Tp_RMASS, int Tp_ZERO>
    struct FixLangevinKokkosPostForceFunctor {

      typedef typename GetDeviceType<Space>::value DeviceType;
      typedef FSUM value_type;
      FixLangevinKokkos<Space> c;

    FixLangevinKokkosPostForceFunctor(FixLangevinKokkos<Space>* c_ptr):
      c(*c_ptr) {}
      ~FixLangevinKokkosPostForceFunctor(){c.cleanup_copy();}

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
      static void init(volatile value_type &update) {
        update.fx = 0.0;
        update.fy = 0.0;
        update.fz = 0.0;
      }
      KOKKOS_INLINE_FUNCTION
      static void join(volatile value_type &update,
                       const volatile value_type &source) {
        update.fx += source.fx;
        update.fy += source.fy;
        update.fz += source.fz;
      }

    };

  template <ExecutionSpace Space>
    struct FixLangevinKokkosZeroForceFunctor {
      typedef typename GetDeviceType<Space>::value DeviceType;
      FixLangevinKokkos<Space> c;

    FixLangevinKokkosZeroForceFunctor(FixLangevinKokkos<Space>* c_ptr):
      c(*c_ptr) {c.cleanup_copy();}

      KOKKOS_INLINE_FUNCTION
      void operator()(const int i) const {
        c.zero_force_item(i);
      }
    };

  template<ExecutionSpace Space>
    struct FixLangevinKokkosTallyEnergyFunctor {
      typedef typename GetDeviceType<Space>::value DeviceType;
      FixLangevinKokkos<Space> c;
      typedef KK_FLOAT value_type;
    FixLangevinKokkosTallyEnergyFunctor(FixLangevinKokkos<Space>* c_ptr):
      c(*c_ptr) {c.cleanup_copy();}

      KOKKOS_INLINE_FUNCTION
      void operator()(const int i, value_type &energy) const {
        energy += c.compute_energy_item(i);
      }
      KOKKOS_INLINE_FUNCTION
      static void init(volatile value_type &update) {
        update = 0.0;
      }
      KOKKOS_INLINE_FUNCTION
      static void join(volatile value_type &update,
                       const volatile value_type &source) {
        update += source;
      }
    };

  template <ExecutionSpace Space, int RMass>
  struct FixLangevinKokkosEndOfStepFunctor {
    typedef typename GetDeviceType<Space>::value DeviceType;
    FixLangevinKokkos<Space> c;

    FixLangevinKokkosEndOfStepFunctor(FixLangevinKokkos<Space>* c_ptr):
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

/* ERROR/WARNING messages:

E: Fix langevin omega is not yet implemented with kokkos

This option is not yet available.

E: Fix langevin angmom is not yet implemented with kokkos

This option is not yet available.

E: Cannot zero Langevin force of 0 atoms

The group has zero atoms, so you cannot request its force
be zeroed.

E: Fix langevin variable returned negative temperature

Self-explanatory.

E: Fix langevin gjf with tbias is not yet implemented with kokkos

This option is not yet available.

W: Fix langevin gjf using random gaussians is not implemented with kokkos

This will most likely cause errors in kinetic fluctuations.

*/
