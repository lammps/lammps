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

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(pppm/tip4p/kk,PPPMTIP4PKokkos<LMPDeviceType>);
KSpaceStyle(pppm/tip4p/kk/device,PPPMTIP4PKokkos<LMPDeviceType>);
KSpaceStyle(pppm/tip4p/kk/host,PPPMTIP4PKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_PPPM_TIP4P_KOKKOS_H
#define LMP_PPPM_TIP4P_KOKKOS_H

#include "kokkos_type.h"
#include "pppm_kokkos.h"

namespace LAMMPS_NS {

struct TagPPPMTIP4P_make_rho_atomic {};
struct TagPPPMTIP4P_make_rho {};
struct TagPPPMTIP4P_slabcorr1 {};
struct TagPPPMTIP4P_slabcorr2 {};
struct TagPPPMTIP4P_slabcorr3 {};
struct TagPPPMTIP4P_slabcorr4 {};

template<class DeviceType>
class PPPMTIP4PKokkos : public PPPMKokkos<DeviceType> {
 public:
  typedef ArrayTypes<DeviceType> AT;
  typedef PPPMKokkos<DeviceType> Base;

  PPPMTIP4PKokkos(class LAMMPS *);
  ~PPPMTIP4PKokkos() override;
  void init() override;
  double memory_usage() override;

  void pp_pre_particle_map() override;

  void particle_map() override;
  void make_rho() override;
  void fieldforce_ik() override;
  void fieldforce_peratom() override;
  void slabcorr() override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_particle_map, const int &i) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_make_rho_atomic, const int &i) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_make_rho,
                  typename Kokkos::TeamPolicy<DeviceType, TagPPPMTIP4P_make_rho>::member_type dev) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_fieldforce_ik, const int &i) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_fieldforce_peratom, const int &i) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_slabcorr1, const int &i, double &dipole) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_slabcorr2, const int &i, double &dipole_r2) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_slabcorr3, const int &i) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_slabcorr4, const int &i) const;

 private:
  void find_M_host(int, int &, int &, double *);
  void tip4p_preprocess_host();

  DAT::tdual_kkfloat_1d_3 k_xM;
  DAT::tdual_int_1d k_ih1, k_ih2;

  typename AT::t_kkfloat_1d_3_lr_randomread d_xM;
  typename AT::t_int_1d_randomread d_ih1, d_ih2;
  typename AT::t_int_1d_randomread d_type;

  KK_FLOAT alpha_kk;
};

}    // namespace LAMMPS_NS

#endif
#endif
