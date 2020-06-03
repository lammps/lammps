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

#ifdef PAIR_CLASS

PairStyle(dpd/fdt/energy/kk,PairDPDfdtEnergyKokkos<Device>)
PairStyle(dpd/fdt/energy/kk/device,PairDPDfdtEnergyKokkos<Device>)
PairStyle(dpd/fdt/energy/kk/host,PairDPDfdtEnergyKokkos<Host>)

#else

#ifndef LMP_PAIR_DPD_FDT_ENERGY_KOKKOS_H
#define LMP_PAIR_DPD_FDT_ENERGY_KOKKOS_H

#if !defined(DPD_USE_RAN_MARS) && !defined(DPD_USE_Random_XorShift64) && !defined(Random_XorShift1024)
#define DPD_USE_Random_XorShift64
#endif

#include "pair_dpd_fdt_energy.h"
#include "pair_kokkos.h"
#include "kokkos_type.h"
#ifdef DPD_USE_RAN_MARS
#include "rand_pool_wrap_kokkos.h"
#else
#include "Kokkos_Random.hpp"
#endif

namespace LAMMPS_NS {

struct TagPairDPDfdtEnergyZero{};

template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool STACKPARAMS>
struct TagPairDPDfdtEnergyComputeSplit{};

template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool STACKPARAMS>
struct TagPairDPDfdtEnergyComputeNoSplit{};

template<ExecutionSpace Space>
class PairDPDfdtEnergyKokkos : public PairDPDfdtEnergy {
 public:
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;
  typedef EV_FLOAT value_type;

  PairDPDfdtEnergyKokkos(class LAMMPS *);
  virtual ~PairDPDfdtEnergyKokkos();
  virtual void compute(int, int);
  void init_style();
  double init_one(int, int);

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairDPDfdtEnergyZero, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool STACKPARAMS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairDPDfdtEnergyComputeSplit<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool STACKPARAMS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairDPDfdtEnergyComputeSplit<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool STACKPARAMS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairDPDfdtEnergyComputeNoSplit<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool STACKPARAMS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairDPDfdtEnergyComputeNoSplit<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                  const KK_FLOAT &dely, const KK_FLOAT &delz) const;

  KOKKOS_INLINE_FUNCTION
  int sbmask(const int& j) const;

  struct params_dpd {
    KOKKOS_INLINE_FUNCTION
    params_dpd(){cut=0;a0=0;sigma=0;kappa=0;alpha=0;};
    KOKKOS_INLINE_FUNCTION
    params_dpd(int i){cut=0;a0=0;sigma=0;kappa=0;alpha=0;};
    KK_FLOAT cut,a0,sigma,kappa,alpha;
  };

  DAT::tdual_float_1d k_duCond,k_duMech;

#ifdef DPD_USE_RAN_MARS
  RandPoolWrap rand_pool;
  typedef RandWrap rand_type;
#elif defined(DPD_USE_Random_XorShift64)
  Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
  typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;
#elif defined(DPD_USE_Random_XorShift1024)
  Kokkos::Random_XorShift1024_Pool<DeviceType> rand_pool;
  typedef typename Kokkos::Random_XorShift1024_Pool<DeviceType>::generator_type rand_type;
#endif

  DAT::tdual_float_2d k_cutsq;
  typename AT::t_float_2d d_cutsq;

 protected:
  int eflag,vflag;
  int nlocal,neighflag;
  KK_FLOAT dtinvsqrt;
  KK_FLOAT boltz,ftm2v;
  KK_FLOAT special_lj[4];

  virtual void allocate();

  Kokkos::DualView<params_dpd**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_dpd**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  // hardwired to space for MAX_TYPES_STACKPARAMS (12) atom types
  params_dpd m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_float_1d_3_randomread x;
  typename AT::t_float_1d_3 c_x;
  typename AT::t_float_1d_3_randomread v;
  typename AT::t_float_1d_3 f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_float_1d_randomread mass;
  typename AT::t_float_1d rmass;
  typename AT::t_float_1d dpdTheta;
  typename AT::t_float_1d d_duCond,d_duMech;
  HAT::t_float_1d h_duCond,h_duMech;

  DAT::tdual_float_1d k_eatom;
  DAT::tdual_float_1d_6 k_vatom;
  typename AT::t_float_1d d_eatom;
  typename AT::t_float_1d_6 d_vatom;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  friend void pair_virial_fdotr_compute<Space,PairDPDfdtEnergyKokkos>(PairDPDfdtEnergyKokkos*);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair dpd/fdt/energy requires ghost atoms store velocity

Use the communicate vel yes command to enable this.

E: Pair dpd/fdt/energy requires newton pair on

Self-explanatory.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

*/
