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

/* ---------------------------------------------------------------------- */
/* Tags for TIP4P-specific kernels.                                        */
/* The two-phase split avoids data races without requiring atomics:        */
/*   phase1 - each atom writes only to its own index                       */
/*   phase2 - O atoms write to their H atoms' indices                     */
/*            (safe because each H belongs to exactly one O)               */
/* ---------------------------------------------------------------------- */
struct TagPPPMTIP4P_make_rho_atomic {};
struct TagPPPMTIP4P_make_rho {};

struct TagPPPMTIP4P_fieldforce_ik_phase1 {};   // own force: O gets (1-alpha), H gets full
struct TagPPPMTIP4P_fieldforce_ik_phase2 {};   // O scatters alpha/2 to each H

struct TagPPPMTIP4P_fieldforce_peratom_phase2 {}; // O scatters eatom/vatom alpha/2 to each H

struct TagPPPMTIP4P_slabcorr1 {};              // dipole sum (TIP4P: M-site for O)
struct TagPPPMTIP4P_slabcorr2 {};              // dipole_r2 sum
struct TagPPPMTIP4P_slabcorr3 {};              // per-atom slab energy (uses atom z)
struct TagPPPMTIP4P_slabcorr4_phase1 {};       // own z-force: O gets (1-alpha), H gets full
struct TagPPPMTIP4P_slabcorr4_phase2 {};       // O scatters alpha/2 z-force to each H

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

  /* particle_map functor: uses M-site positions from d_xM */
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_particle_map, const int &i) const;

  /* make_rho functors: use d_xM for dx/dy/dz */
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_make_rho_atomic, const int &i) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_make_rho,
                  typename Kokkos::TeamPolicy<DeviceType, TagPPPMTIP4P_make_rho>::member_type) const;

  /* fieldforce_ik: two-phase, no atomics needed */
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_fieldforce_ik_phase1, const int &i) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_fieldforce_ik_phase2, const int &i) const;

  /* fieldforce_peratom: phase1 inherited from base (handles all atoms' own contribution);
     phase2 is TIP4P-specific to scatter O's contribution to H atoms */
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPM_fieldforce_peratom, const int &i) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_fieldforce_peratom_phase2, const int &i) const;

  /* slabcorr functors */
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_slabcorr1, const int &i, double &dipole) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_slabcorr2, const int &i, double &dipole_r2) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_slabcorr3, const int &i) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_slabcorr4_phase1, const int &i) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_slabcorr4_phase2, const int &i) const;

 private:
  void find_M_host(int i, int &iH1, int &iH2, double *xM);
  void tip4p_preprocess_host();

  /* host-side dual views for M-site positions and H-atom indices */
  DAT::tdual_kkfloat_1d_3 k_xM;
  DAT::tdual_int_1d k_ih1, k_ih2;

  /* device-side read-only views set each compute() call */
  typename AT::t_kkfloat_1d_3_lr_randomread d_xM;
  typename AT::t_int_1d_randomread d_ih1, d_ih2;
  typename AT::t_int_1d_randomread d_type;

  KK_FLOAT alpha_kk;  // geometric factor (device-precision copy of alpha)
};

}    // namespace LAMMPS_NS

#endif // !LMP_PPPM_TIP4P_KOKKOS_H
#endif
