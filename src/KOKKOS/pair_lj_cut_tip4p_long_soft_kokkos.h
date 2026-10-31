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

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/tip4p/long/soft/kk,PairLJCutTIP4PLongSoftKokkos<LMPDeviceType>);
PairStyle(lj/cut/tip4p/long/soft/kk/device,PairLJCutTIP4PLongSoftKokkos<LMPDeviceType>);
PairStyle(lj/cut/tip4p/long/soft/kk/host,PairLJCutTIP4PLongSoftKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LJ_CUT_TIP4P_LONG_SOFT_KOKKOS_H
#define LMP_PAIR_LJ_CUT_TIP4P_LONG_SOFT_KOKKOS_H

#include "pair_lj_cut_tip4p_long_soft.h"
#include "pair_tip4p_kokkos.h"

namespace LAMMPS_NS {

template<int EVFLAG>
struct TagPairLJCutTIP4PLongSoftCompute{};

template<class DeviceType>
class PairLJCutTIP4PLongSoftKokkos : public PairTIP4PKokkos<DeviceType,PairLJCutTIP4PLongSoft> {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;
  typedef PairTIP4PKokkos<DeviceType,PairLJCutTIP4PLongSoft> Base;

  PairLJCutTIP4PLongSoftKokkos(class LAMMPS *lmp) : Base(lmp) { this->respa_enable = 0; }

  void compute(int, int) override;
  double init_one(int, int) override;

  using Base::operator();

  template<int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairLJCutTIP4PLongSoftCompute<EVFLAG>, const int&, EV_FLOAT&) const;

  template<int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairLJCutTIP4PLongSoftCompute<EVFLAG>, const int&) const;

 protected:
  // per-type-pair soft-core coefficients (device); re-synced when init_one() ran
  typename AT::t_kkfloat_2d d_lj1, d_lj2, d_lj3, d_lj4, d_epsilon, d_offset, d_cut_ljsq;
  int soft_coeffs_synced = 0;

  using Base::x; using Base::f; using Base::q; using Base::type;
  using Base::d_newsite; using Base::d_hneigh; using Base::d_h_missing;
  using Base::d_neighbors; using Base::d_numneigh; using Base::d_ilist;
  using Base::m_typeO;
  using Base::m_cut_coulsq; using Base::m_cut_coulsqplus;
  using Base::qqrd2e; using Base::special_coul; using Base::special_lj; using Base::g_ewald_kk;
  using Base::sbmask; using Base::ev_tally_tip4p; using Base::apply_site_force;
  using Base::ev_tally;

  // soft-core long-range (Ewald) Coulomb force prefactor and energy
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT coul_long_soft(const KK_FLOAT rsq, const KK_FLOAT qtmp, const KK_FLOAT qj,
                          const KK_FLOAT factor_coul, const KK_FLOAT lam1, const KK_FLOAT lam2,
                          const bool want_ecoul, KK_FLOAT &ecoul) const
  {
    const KK_FLOAT r = Kokkos::sqrt(rsq);
    const KK_FLOAT grij = g_ewald_kk * r;
    const KK_FLOAT expm2 = Kokkos::exp(-grij*grij);
    const KK_FLOAT t = (KK_FLOAT)1.0 / ((KK_FLOAT)1.0 + (KK_FLOAT)EwaldConst::EWALD_P*grij);
    const KK_FLOAT erfc = t*((KK_FLOAT)EwaldConst::A1+t*((KK_FLOAT)EwaldConst::A2+
                          t*((KK_FLOAT)EwaldConst::A3+t*((KK_FLOAT)EwaldConst::A4+
                          t*(KK_FLOAT)EwaldConst::A5))))*expm2;
    const KK_FLOAT denc = Kokkos::sqrt(lam2 + rsq);
    const KK_FLOAT prefactor = qqrd2e * lam1 * qtmp * qj / (denc*denc*denc);

    // the soft core replaces only the 1/r factor, while the Ewald damping keeps
    // the true distance, so differentiating erfc(g*r)/denc leaves the exponential
    // term scaled by denc^2/r^2 against a plain 1/r kernel.  already a force/r
    // quantity, as in the base style

    KK_FLOAT forcecoul = prefactor * (erfc + (KK_FLOAT)EwaldConst::EWALD_F*grij*expm2*denc*denc/rsq);
    if (factor_coul < (KK_FLOAT)1.0) forcecoul -= ((KK_FLOAT)1.0-factor_coul)*prefactor;

    if (want_ecoul) {
      const KK_FLOAT eprefactor = qqrd2e * lam1 * qtmp * qj / denc;
      ecoul = eprefactor * erfc;
      if (factor_coul < (KK_FLOAT)1.0) ecoul -= ((KK_FLOAT)1.0-factor_coul)*eprefactor;
    }
    return forcecoul;
  }
};

}

#endif
#endif
