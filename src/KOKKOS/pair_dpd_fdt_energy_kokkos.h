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

PairStyle(dpd/fdt/energy/kk,PairDPDfdtEnergyKokkos<LMPDeviceType>)
PairStyle(dpd/fdt/energy/kk/device,PairDPDfdtEnergyKokkos<LMPDeviceType>)
PairStyle(dpd/fdt/energy/kk/host,PairDPDfdtEnergyKokkos<LMPHostType>)

#else

#ifndef LMP_PAIR_DPD_FDT_ENERGY_KOKKOS_H
#define LMP_PAIR_DPD_FDT_ENERGY_KOKKOS_H

#include "pair_dpd_fdt_energy.h"
#include "pair_kokkos.h"
#include "kokkos_type.h"
#include "Kokkos_Random.hpp"
#include "rand_pool_wrap_kokkos.h"

namespace LAMMPS_NS {

struct TagPairDPDfdtEnergyZero{};

template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairDPDfdtEnergyComputeSplit{};

template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairDPDfdtEnergyComputeNoSplit{};

template<class DeviceType>
class PairDPDfdtEnergyKokkos : public PairDPDfdtEnergy {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairDPDfdtEnergyKokkos(class LAMMPS *);
  virtual ~PairDPDfdtEnergyKokkos();
  virtual void compute(int, int);
  void init_style();
  double init_one(int, int);

  void operator()(TagPairDPDfdtEnergyZero, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairDPDfdtEnergyComputeSplit<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairDPDfdtEnergyComputeSplit<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairDPDfdtEnergyComputeNoSplit<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairDPDfdtEnergyComputeNoSplit<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &epair, const F_FLOAT &fpair, const F_FLOAT &delx,
                  const F_FLOAT &dely, const F_FLOAT &delz) const;

  KOKKOS_INLINE_FUNCTION
  int sbmask(const int& j) const;

  struct params_dpd {
    params_dpd(){cut=0;a0=0;sigma=0;kappa=0;};
    params_dpd(int i){cut=0;a0=0;sigma=0;kappa=0;};
    F_FLOAT cut,a0,sigma,kappa;
  };

 protected:
  int eflag,vflag;
  int nlocal,neighflag;
  int STACKPARAMS;
  double dtinvsqrt;
  double boltz;

  virtual void allocate();

  Kokkos::DualView<params_dpd**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_dpd**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  // hardwired to space for 15 atom types
  params_dpd m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  F_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  F_FLOAT m_cut[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename ArrayTypes<DeviceType>::t_x_array_randomread x;
  typename ArrayTypes<DeviceType>::t_x_array c_x;
  typename ArrayTypes<DeviceType>::t_v_array_randomread v;
  typename ArrayTypes<DeviceType>::t_f_array f;
  typename ArrayTypes<DeviceType>::t_int_1d_randomread type;
  typename ArrayTypes<DeviceType>::t_float_1d_randomread mass;
  double *rmass;
  typename AT::t_efloat_1d dpdTheta;
  DAT::tdual_efloat_1d k_duCond,k_duMech;
  typename AT::t_efloat_1d d_duCond,d_duMech;
  HAT::t_efloat_1d h_duCond,h_duMech;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  DAT::t_efloat_1d d_eatom;
  DAT::t_virial_array d_vatom;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  typename ArrayTypes<DeviceType>::tdual_ffloat_2d k_cutsq;
  typename ArrayTypes<DeviceType>::t_ffloat_2d d_cutsq;

  /**/Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
  typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;/**/

  /**RandPoolWrap rand_pool;
  typedef RandWrap rand_type;/**/

  friend void pair_virial_fdotr_compute<PairDPDfdtEnergyKokkos>(PairDPDfdtEnergyKokkos*);
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
