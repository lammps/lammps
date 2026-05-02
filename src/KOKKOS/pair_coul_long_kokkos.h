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

PairStyle(coul/long/kk,        PairCoulLongKokkos<LMPDeviceType,false,false>);
PairStyle(coul/long/kk/device, PairCoulLongKokkos<LMPDeviceType,false,false>);
PairStyle(coul/long/kk/host,   PairCoulLongKokkos<LMPHostType,  false,false>);

PairStyle(tip4p/long/kk,        PairCoulLongKokkos<LMPDeviceType,true,false>);
PairStyle(tip4p/long/kk/device, PairCoulLongKokkos<LMPDeviceType,true,false>);
PairStyle(tip4p/long/kk/host,   PairCoulLongKokkos<LMPHostType,  true,false>);

/*
PairStyle(coul/long/soft/kk,        PairCoulLongKokkos<LMPDeviceType,false,true>);
PairStyle(coul/long/soft/kk/device, PairCoulLongKokkos<LMPDeviceType,false,true>);
PairStyle(coul/long/soft/kk/host,   PairCoulLongKokkos<LMPHostType,  false,true>);


PairStyle(tip4p/long/soft/kk,        PairCoulLongKokkos<LMPDeviceType,true,true>);
PairStyle(tip4p/long/soft/kk/device, PairCoulLongKokkos<LMPDeviceType,true,true>);
PairStyle(tip4p/long/soft/kk/host,   PairCoulLongKokkos<LMPHostType,  true,true>);
*/

// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_COUL_LONG_KOKKOS_H
#define LMP_PAIR_COUL_LONG_KOKKOS_H

#include "pair_kokkos.h"
#include "neigh_list_kokkos.h"
#include "fix.h"

#include "pair_coul_long.h"
#include "pair_tip4p_long.h"
#ifdef LMP_FEP
  #include "pair_coul_long_soft.h"
  #include "pair_tip4p_long_soft.h"
#endif

namespace LAMMPS_NS {

template<class DeviceType, bool TIP4P, bool SOFT>
class PairCoulLongKokkos :
#ifdef LMP_FEP
  public std::conditional_t<TIP4P,
    std::conditional_t<SOFT,PairTIP4PLongSoft,PairTIP4PLong>,
    std::conditional_t<SOFT,PairCoulLongSoft,PairCoulLong>
  >
#else
  public std::conditional_t<TIP4P,PairTIP4PLong,PairCoulLong>
#endif
{
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=1};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairCoulLongKokkos(class LAMMPS *);
  ~PairCoulLongKokkos() override;

  void compute(int, int) override;

  void init_tables(double cut_coul, double *cut_respa) override;
  void init_style() override;
  double init_one(int, int) override;

  struct params_coul{
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_coul() {cut_coulsq=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_coul(int /*i*/) {cut_coulsq=0;};
    KK_FLOAT cut_coulsq;
  };

 protected:

  using Pointers::atom;
  using Pointers::atomKK;
  using Pointers::error;
  using Pointers::force;
  using Pointers::lmp;
  using Pointers::memory;
  using Pointers::memoryKK;
  using Pointers::neighbor;
  using Pointers::update;

  using Pair::copymode;
  using Pair::execution_space;

  using Pair::eatom;
  using Pair::eflag_atom;
  using Pair::vatom;
  using Pair::cutsq;

  using Pair::ncoulshiftbits;
using Pair::vflag_atom;
using Pair::vflag_fdotr;

using Pair::vflag_global;
using Pair::virial;
using Pair::eng_vdwl;
using Pair::eng_coul;

using Pair::list;

  using Pair::rtable;
  using Pair::drtable;
  using Pair::ftable;
  using Pair::dftable;
  using Pair::ctable;
  using Pair::dctable;
  using Pair::etable;
  using Pair::detable;
  using Pair::ptable;
  using Pair::dptable;
  using Pair::vtable;
  using Pair::dvtable;

  using Pair::ev_init;
  using Pair::maxeatom;
  using Pair::maxvatom;
  using Pair::ncoultablebits;
  using Pair::no_virial_fdotr_compute;
using Pair::union_int_float_t;

  using Pair::respa_enable;
  using Pair::kokkosable;
  using Pair::datamask_modify;
  using Pair::datamask_read;
  using Pair::allocated;
  using PairCoulLong::cut_coulsq;
  using PairCoulLong::g_ewald;

  using Pair::ncoulmask;
  using Pair::tabinnersq;




  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fpair(const KK_FLOAT& /*rsq*/, const int& /*i*/, const int& /*j*/,
                        const int& /*itype*/, const int& /*jtype*/) const { return 0.0; }

  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fcoul(const KK_FLOAT& rsq, const int& i, const int&j, const int& itype,
                        const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const;

  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_evdwl(const KK_FLOAT& /*rsq*/, const int& /*i*/, const int& /*j*/,
                        const int& /*itype*/, const int& /*jtype*/) const { return 0; }

  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_ecoul(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const;

  Kokkos::DualView<params_coul**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_coul**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  // hardwired to space for 12 atom types
  params_coul m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cut_ljsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cut_coulsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_3_lr c_x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_1d_randomread q;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  int newton_pair;

  DAT::ttransform_kkfloat_2d k_cutsq;
  typename AT::t_kkfloat_2d d_cutsq;
  typename AT::t_kkfloat_2d d_cut_ljsq;
  typename AT::t_kkfloat_2d d_cut_coulsq;

  typename AT::t_kkfloat_1d_randomread
    d_rtable, d_drtable, d_ftable, d_dftable,
    d_ctable, d_dctable, d_etable, d_detable;

  int neighflag;
  int nlocal,nall,eflag,vflag;

  KK_FLOAT special_lj[4], special_coul[4];
  KK_FLOAT qqrd2e;

  void allocate() override;

  friend struct PairComputeFunctor<PairCoulLongKokkos,FULL,true,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,FULL,true,1,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,HALF,true,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,HALFTHREAD,true,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,FULL,false,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,FULL,false,1,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,HALF,false,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,HALFTHREAD,false,0,CoulLongTable<1>>;
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongKokkos,FULL,0,CoulLongTable<1>>(PairCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongKokkos,FULL,1,CoulLongTable<1>>(PairCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongKokkos,HALF,0,CoulLongTable<1>>(PairCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongKokkos,HALFTHREAD,0,CoulLongTable<1>>(PairCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairCoulLongKokkos,CoulLongTable<1>>(PairCoulLongKokkos*,
                                                            NeighListKokkos<DeviceType>*);
  friend struct PairComputeFunctor<PairCoulLongKokkos,FULL,true,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,FULL,true,1,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,HALF,true,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,HALFTHREAD,true,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,FULL,false,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,FULL,false,1,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,HALF,false,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,HALFTHREAD,false,0,CoulLongTable<0>>;
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongKokkos,FULL,0,CoulLongTable<0>>(PairCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongKokkos,FULL,1,CoulLongTable<0>>(PairCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongKokkos,HALF,0,CoulLongTable<0>>(PairCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongKokkos,HALFTHREAD,0,CoulLongTable<0>>(PairCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairCoulLongKokkos,CoulLongTable<0>>(PairCoulLongKokkos*,
                                                            NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairCoulLongKokkos>(PairCoulLongKokkos*);

};

}

#endif // !LMP_PAIR_COUL_LONG_KOKKOS_H
#endif

