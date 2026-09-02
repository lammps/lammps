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
PairStyle(coul/long/cs/kk,PairCoulLongCSKokkos<LMPDeviceType>);
PairStyle(coul/long/cs/kk/device,PairCoulLongCSKokkos<LMPDeviceType>);
PairStyle(coul/long/cs/kk/host,PairCoulLongCSKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_COUL_LONG_CS_KOKKOS_H
#define LMP_PAIR_COUL_LONG_CS_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_coul_long_cs.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairCoulLongCSKokkos : public PairCoulLongCS {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=1};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairCoulLongCSKokkos(class LAMMPS *);
  ~PairCoulLongCSKokkos() override;

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

  friend struct PairComputeFunctor<PairCoulLongCSKokkos,FULL,true,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairCoulLongCSKokkos,FULL,true,1,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairCoulLongCSKokkos,HALF,true,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairCoulLongCSKokkos,HALFTHREAD,true,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairCoulLongCSKokkos,FULL,false,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairCoulLongCSKokkos,FULL,false,1,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairCoulLongCSKokkos,HALF,false,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairCoulLongCSKokkos,HALFTHREAD,false,0,CoulLongTable<1>>;
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongCSKokkos,FULL,0,CoulLongTable<1>>(PairCoulLongCSKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongCSKokkos,FULL,1,CoulLongTable<1>>(PairCoulLongCSKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongCSKokkos,HALF,0,CoulLongTable<1>>(PairCoulLongCSKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongCSKokkos,HALFTHREAD,0,CoulLongTable<1>>(PairCoulLongCSKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairCoulLongCSKokkos,CoulLongTable<1>>(PairCoulLongCSKokkos*,
                                                            NeighListKokkos<DeviceType>*);
  friend struct PairComputeFunctor<PairCoulLongCSKokkos,FULL,true,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairCoulLongCSKokkos,FULL,true,1,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairCoulLongCSKokkos,HALF,true,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairCoulLongCSKokkos,HALFTHREAD,true,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairCoulLongCSKokkos,FULL,false,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairCoulLongCSKokkos,FULL,false,1,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairCoulLongCSKokkos,HALF,false,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairCoulLongCSKokkos,HALFTHREAD,false,0,CoulLongTable<0>>;
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongCSKokkos,FULL,0,CoulLongTable<0>>(PairCoulLongCSKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongCSKokkos,FULL,1,CoulLongTable<0>>(PairCoulLongCSKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongCSKokkos,HALF,0,CoulLongTable<0>>(PairCoulLongCSKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongCSKokkos,HALFTHREAD,0,CoulLongTable<0>>(PairCoulLongCSKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairCoulLongCSKokkos,CoulLongTable<0>>(PairCoulLongCSKokkos*,
                                                            NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairCoulLongCSKokkos>(PairCoulLongCSKokkos*);

};

}

#endif
#endif

