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
FixStyle(shardlow/kk,FixShardlowKokkos<LMPDeviceType>);
FixStyle(shardlow/kk/device,FixShardlowKokkos<LMPDeviceType>);
FixStyle(shardlow/kk/host,FixShardlowKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_SHARDLOW_KOKKOS_H
#define LMP_FIX_SHARDLOW_KOKKOS_H

#include "fix_shardlow.h"
#include "kokkos_type.h"
#include "neigh_list_kokkos.h"
#ifdef ENABLE_KOKKOS_DPD_CONSTANT_TEMPERATURE
#include "pair_dpd_fdt_kokkos.h"
#endif
#include "pair_dpd_fdt_energy_kokkos.h"

namespace LAMMPS_NS {

template<bool STACKPARAMS>
struct TagFixShardlowSSAUpdateDPDE{};

template<bool STACKPARAMS>
struct TagFixShardlowSSAUpdateDPDEGhost{};

template<class DeviceType>
class FixShardlowKokkos : public FixShardlow {
 public:
  typedef ArrayTypes<DeviceType> AT;
  NeighListKokkos<DeviceType> *k_list; // The SSA specific neighbor list

  FixShardlowKokkos(class LAMMPS *, int, char **);
  ~FixShardlowKokkos();
  int setmask();
  virtual void init();
  virtual void init_list(int, class NeighList *);
  virtual void initial_integrate(int);
  void setup_pre_neighbor();
  void pre_neighbor();

  double memory_usage();

  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  int pack_forward_comm(int , int *, double *, int, int *);
  void unpack_forward_comm(int , int , double *);

  struct params_ssa {
    KOKKOS_INLINE_FUNCTION
    params_ssa() {cutinv=FLT_MAX;halfsigma=0;kappa=0;alpha=0;};
    KOKKOS_INLINE_FUNCTION
    params_ssa(int /*i*/) {cutinv=FLT_MAX;halfsigma=0;kappa=0;alpha=0;};
    F_FLOAT cutinv,halfsigma,kappa,alpha;
  };

  template<bool STACKPARAMS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixShardlowSSAUpdateDPDE<STACKPARAMS>, const int&) const;

  template<bool STACKPARAMS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixShardlowSSAUpdateDPDEGhost<STACKPARAMS>, const int&) const;

#ifdef DEBUG_SSA_PAIR_CT
  typename AT::t_int_2d d_counters;
  typename HAT::t_int_2d h_counters;
  typename AT::t_int_1d d_hist;
  typename HAT::t_int_1d h_hist;
#endif

 protected:
  int workPhase;
  double theta_ij_inv,boltz_inv,ftm2v,dt;

#ifdef ENABLE_KOKKOS_DPD_CONSTANT_TEMPERATURE
//  class PairDPDfdt *pairDPD; FIXME as per k_pairDPDE below
#endif
  PairDPDfdtEnergyKokkos<DeviceType> *k_pairDPDE;

  Kokkos::DualView<params_ssa**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_ssa**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  // hardwired to space for MAX_TYPES_STACKPARAMS (12) atom types
  params_ssa m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  F_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename ArrayTypes<DeviceType>::t_ffloat_2d d_cutsq;

  typename DAT::tdual_v_array k_v_t0;
  // typename AT::t_v_array d_v_t0; v_t0 only used in comm routines (on host)
  typename HAT::t_v_array h_v_t0;

  typename AT::t_x_array x;
  typename AT::t_v_array v;
  typename HAT::t_v_array h_v;
  typename AT::t_efloat_1d uCond, uMech;
  typename HAT::t_efloat_1d h_uCond, h_uMech;
  typename AT::t_int_1d type;
  bool massPerI;
  typename AT::t_float_1d_randomread masses;
  typename AT::t_efloat_1d dpdTheta;

  // Storage for the es_RNG state variables
  typedef Kokkos::View<random_external_state::es_RNG_t*,DeviceType> es_RNGs_type;
  es_RNGs_type d_rand_state;

  double dtsqrt; // = sqrt(update->dt);
  int ghostmax;
  int nlocal, nghost;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist, d_numneigh;

  int ssa_phaseCt;
  typename AT::t_int_1d ssa_phaseLen;
  typename AT::t_int_2d ssa_itemLoc, ssa_itemLen;

  int ssa_gphaseCt;
  typename AT::t_int_1d ssa_gphaseLen;
  typename AT::t_int_2d ssa_gitemLoc, ssa_gitemLen;


#ifdef ENABLE_KOKKOS_DPD_CONSTANT_TEMPERATURE
  template<bool STACKPARAMS>
  KOKKOS_INLINE_FUNCTION
  void ssa_update_dpd(int, int, int) const;  // Constant Temperature
#endif
  template<bool STACKPARAMS>
  KOKKOS_INLINE_FUNCTION
  void ssa_update_dpde(int, int, int) const; // Constant Energy

};

}

#endif
#endif

