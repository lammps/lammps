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

#ifdef DIHEDRAL_CLASS
// clang-format off
DihedralStyle(cosine/shift/exp/kk,DihedralCosineShiftExpKokkos<LMPDeviceType>);
DihedralStyle(cosine/shift/exp/kk/device,DihedralCosineShiftExpKokkos<LMPDeviceType>);
DihedralStyle(cosine/shift/exp/kk/host,DihedralCosineShiftExpKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_DIHEDRAL_COSINE_SHIFT_EXP_KOKKOS_H
#define LMP_DIHEDRAL_COSINE_SHIFT_EXP_KOKKOS_H

#include "dihedral_cosine_shift_exp.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int NEWTON_BOND, int EVFLAG>
struct TagDihedralCosineShiftExpCompute{};

template<class DeviceType>
class DihedralCosineShiftExpKokkos : public DihedralCosineShiftExp {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  DihedralCosineShiftExpKokkos(class LAMMPS *);
  ~DihedralCosineShiftExpKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void read_restart(FILE *) override;

  template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagDihedralCosineShiftExpCompute<NEWTON_BOND,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagDihedralCosineShiftExpCompute<NEWTON_BOND,EVFLAG>, const int&) const;

  //template<int NEWTON_BOND>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int i1, const int i2, const int i3, const int i4,
                          KK_FLOAT &edihedral, KK_FLOAT *f1, KK_FLOAT *f3, KK_FLOAT *f4,
                          const KK_FLOAT &vb1x, const KK_FLOAT &vb1y, const KK_FLOAT &vb1z,
                          const KK_FLOAT &vb2x, const KK_FLOAT &vb2y, const KK_FLOAT &vb2z,
                          const KK_FLOAT &vb3x, const KK_FLOAT &vb3y, const KK_FLOAT &vb3z) const;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;

 protected:

  class NeighborKokkos *neighborKK;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_2d_lr dihedrallist;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  int nlocal,newton_bond;
  int eflag,vflag;

  DAT::tdual_int_scalar k_warning_flag;
  typename AT::t_int_scalar d_warning_flag;
  HAT::t_int_scalar h_warning_flag;

  DAT::tdual_kkfloat_1d k_umin;
  DAT::tdual_kkfloat_1d k_a;
  DAT::tdual_kkfloat_1d k_opt1;
  DAT::tdual_kkfloat_1d k_cost;
  DAT::tdual_kkfloat_1d k_sint;
  DAT::tdual_int_1d k_doExpansion;

  typename AT::t_kkfloat_1d d_umin;
  typename AT::t_kkfloat_1d d_a;
  typename AT::t_kkfloat_1d d_opt1;
  typename AT::t_kkfloat_1d d_cost;
  typename AT::t_kkfloat_1d d_sint;
  typename AT::t_int_1d d_doExpansion;

  void allocate() override;
};

}

#endif
#endif
