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

#ifdef BOND_CLASS
// clang-format off
BondStyle(fene/kk,BondFENEKokkos<LMPDeviceType>);
BondStyle(fene/kk/device,BondFENEKokkos<LMPDeviceType>);
BondStyle(fene/kk/host,BondFENEKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_BOND_FENE_KOKKOS_H
#define LMP_BOND_FENE_KOKKOS_H

#include "bond_fene.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int NEWTON_BOND, int EVFLAG>
struct TagBondFENECompute{};

template<class DeviceType>
class BondFENEKokkos : public BondFENE {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  BondFENEKokkos(class LAMMPS *);
  ~BondFENEKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void read_restart(FILE *) override;

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagBondFENECompute<NEWTON_BOND,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagBondFENECompute<NEWTON_BOND,EVFLAG>, const int&) const;

  //template<int NEWTON_BOND>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &ebond, const KK_FLOAT &fbond, const KK_FLOAT &delx,
                  const KK_FLOAT &dely, const KK_FLOAT &delz) const;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;

 protected:

  class NeighborKokkos *neighborKK;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_2d_lr bondlist;

  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  typename AT::t_int_scalar d_flag;
  HAT::t_int_scalar h_flag;

  int nlocal,newton_bond;
  int eflag,vflag;

  DAT::tdual_kkfloat_1d k_k;
  DAT::tdual_kkfloat_1d k_r0;
  DAT::tdual_kkfloat_1d k_epsilon;
  DAT::tdual_kkfloat_1d k_sigma;

  typename AT::t_kkfloat_1d d_k;
  typename AT::t_kkfloat_1d d_r0;
  typename AT::t_kkfloat_1d d_epsilon;
  typename AT::t_kkfloat_1d d_sigma;

  void allocate() override;
};

}

#endif
#endif

