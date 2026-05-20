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
DihedralStyle(spherical/kk,DihedralSphericalKokkos<LMPDeviceType>);
DihedralStyle(spherical/kk/device,DihedralSphericalKokkos<LMPDeviceType>);
DihedralStyle(spherical/kk/host,DihedralSphericalKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_DIHEDRAL_SPHERICAL_KOKKOS_H
#define LMP_DIHEDRAL_SPHERICAL_KOKKOS_H

#include "dihedral_spherical.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int NEWTON_BOND, int EVFLAG>
struct TagDihedralSphericalCompute{};

template<class DeviceType>
class DihedralSphericalKokkos : public DihedralSpherical {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  DihedralSphericalKokkos(class LAMMPS *);
  ~DihedralSphericalKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void read_restart(FILE *) override;

  template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagDihedralSphericalCompute<NEWTON_BOND,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagDihedralSphericalCompute<NEWTON_BOND,EVFLAG>, const int&) const;

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

  DAT::tdual_int_1d k_nterms;
  DAT::tdual_kkfloat_2d k_Ccoeff;
  DAT::tdual_kkfloat_2d k_phi_mult;
  DAT::tdual_kkfloat_2d k_phi_shift;
  DAT::tdual_kkfloat_2d k_phi_offset;
  DAT::tdual_kkfloat_2d k_theta1_mult;
  DAT::tdual_kkfloat_2d k_theta1_shift;
  DAT::tdual_kkfloat_2d k_theta1_offset;
  DAT::tdual_kkfloat_2d k_theta2_mult;
  DAT::tdual_kkfloat_2d k_theta2_shift;
  DAT::tdual_kkfloat_2d k_theta2_offset;

  typename AT::t_int_1d d_nterms;
  typename AT::t_kkfloat_2d d_Ccoeff;
  typename AT::t_kkfloat_2d d_phi_mult;
  typename AT::t_kkfloat_2d d_phi_shift;
  typename AT::t_kkfloat_2d d_phi_offset;
  typename AT::t_kkfloat_2d d_theta1_mult;
  typename AT::t_kkfloat_2d d_theta1_shift;
  typename AT::t_kkfloat_2d d_theta1_offset;
  typename AT::t_kkfloat_2d d_theta2_mult;
  typename AT::t_kkfloat_2d d_theta2_shift;
  typename AT::t_kkfloat_2d d_theta2_offset;

  int allocated_kokkos;

  void allocate_kokkos();
  void allocate() override;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT CalcGeneralizedForcesKK(int type, KK_FLOAT phi, KK_FLOAT theta1, KK_FLOAT theta2,
                                   KK_FLOAT &m_du_dth1, KK_FLOAT &m_du_dth2,
                                   KK_FLOAT &m_du_dphi) const;
};

}

#endif
#endif
