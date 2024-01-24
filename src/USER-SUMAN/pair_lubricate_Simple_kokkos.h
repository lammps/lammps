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
PairStyle(lubricate/Simple/kk,PairLubricateSimpleKokkos<LMPDeviceType>);
PairStyle(lubricate/Simple/kk/device,PairLubricateSimpleKokkos<LMPDeviceType>);
PairStyle(lubricate/Simple/kk/host,PairLubricateSimpleKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LUBRICATE_SIMPLE_KOKKOS_H
#define LMP_PAIR_LUBRICATE_SIMPLE_KOKKOS_H

#include "pair_lubricate_Simple.h"
#include "pair_kokkos.h"
//#include "kokkos_type.h"
//#include "kokkos_base.h"
//#include "Kokkos_Random.hpp"
#include "comm_kokkos.h"

namespace LAMMPS_NS {

template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int FLAGFLD>
struct TagPairLubricateSimpleCompute {};

template<int NEIGHFLAG, int NEWTON_PAIR, int SHEARING>
struct TagPairLubricateSimpleComputeFLD {};
  
template<class DeviceType>
class PairLubricateSimpleKokkos : public PairLubricateSimple {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;
  
  PairLubricateSimpleKokkos(class LAMMPS *);
  ~PairLubricateSimpleKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void settings(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int FLAGFLD>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairLubricateSimpleCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,FLAGFLD>, const int, EV_FLOAT &ev) const;
  template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int FLAGFLD>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairLubricateSimpleCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,FLAGFLD>, const int) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int SHEARING>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairLubricateSimpleComputeFLD<NEIGHFLAG, NEWTON_PAIR, SHEARING>, const int, EV_FLOAT &ev) const;
 
  template<int NEIGHFLAG, int NEWTON_PAIR, int SHEARING>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairLubricateSimpleComputeFLD<NEIGHFLAG, NEWTON_PAIR, SHEARING>, const int) const;
  
  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally_xyz(EV_FLOAT &ev, int i, int j,
                    F_FLOAT fx, F_FLOAT fy, F_FLOAT fz,
                    X_FLOAT delx, X_FLOAT dely, X_FLOAT delz) const;
  
  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void v_tally_tensor(EV_FLOAT &ev, int i, int j,
		      F_FLOAT vxx, F_FLOAT vyy, F_FLOAT vzz,
		      F_FLOAT vxy, F_FLOAT vxz, F_FLOAT vyz) const;

 protected:
  typename AT::t_x_array_randomread x;
  typename AT::t_x_array c_x;
  typename AT::t_v_array_randomread v;
  typename AT::t_v_array_randomread omega;
  typename AT::t_f_array f;
  typename AT::t_f_array torque;
  typename AT::t_int_1d_randomread type;
  typename AT::t_float_1d_randomread radius;
  
  DAT::tdual_virial_array k_vatom;
  typename AT::t_virial_array d_vatom;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;
  
  int newton_pair;
  double special_lj[4];

  typename AT::tdual_ffloat_2d k_cutsq;
  typename AT::t_ffloat_2d d_cutsq;
  typename AT::tdual_ffloat_2d k_cut_inner;
  typename AT::t_ffloat_2d d_cut_inner;
  
  int neighflag;
  int nlocal,nall,eflag,vflag;
  LMP_FLOAT vxmu2f;
  
  LMP_FLOAT R0, RT0, RS0;
  LMP_FLOAT xprd, yprd, zprd;
  Few<double, 6> h_rate;
  Few<double, 3> h_ratelo;
  
  class DomainKokkos *domainKK;
  
  void allocate();

  friend void pair_virial_fdotr_compute<PairLubricateSimpleKokkos>(PairLubricateSimpleKokkos*);
};

}

#endif
#endif

