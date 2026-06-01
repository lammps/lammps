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

// clang-format off
#ifndef LMP_PPPM_TIP4P_KOKKOS_H
#define LMP_PPPM_TIP4P_KOKKOS_H

#include "pppm_kokkos.h"
#include <Kokkos_UnorderedMap.hpp>

namespace LAMMPS_NS {

// TIP4P-specific kernels: identical grid math to the base PPPMKokkos, but the
// O atom's charge lives on the off-atom M (virtual) site, so the position used
// for the grid is d_xM and forces/energies are redistributed onto O and its
// two H atoms (Feenstra).

struct TagPPPMTIP4P_findM{};
struct TagPPPMTIP4P_particle_map{};
struct TagPPPMTIP4P_make_rho_atomic{};
struct TagPPPMTIP4P_make_rho{};
struct TagPPPMTIP4P_fieldforce_ik{};
struct TagPPPMTIP4P_fieldforce_peratom{};
struct TagPPPMTIP4P_slabcorr1{};
struct TagPPPMTIP4P_slabcorr4{};

template<class DeviceType>
class PPPMTIP4PKokkos : public PPPMKokkos<DeviceType> {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef FFTArrayTypes<DeviceType> FFT_AT;
  typedef PPPMKokkos<DeviceType> Base;

  PPPMTIP4PKokkos(class LAMMPS *);
  void init() override;

  // un-hide the base-class kernels (declaring our own operator()s would
  // otherwise hide all inherited operator() overloads)
  using Base::operator();

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_findM, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_particle_map, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_make_rho_atomic, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_make_rho, typename Kokkos::TeamPolicy<DeviceType, TagPPPMTIP4P_make_rho>::member_type) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_fieldforce_ik, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_fieldforce_peratom, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_slabcorr1, const int&, double&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPPPMTIP4P_slabcorr4, const int&) const;

 protected:
  void particle_map() override;
  void make_rho() override;
  void fieldforce_ik() override;
  void fieldforce_peratom() override;
  void slabcorr() override;

  // bind TIP4P device data + (re)compute the M sites; shared by particle_map()
  void compute_newsites();

  // find the periodic image of j closest to i (orthogonal box only)
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int closest_image(const int i, int j) const
  {
    if (j < 0) return j;
    const KK_FLOAT xi0 = x(i,0), xi1 = x(i,1), xi2 = x(i,2);
    int closest = j;
    KK_FLOAT delx = xi0 - x(j,0), dely = xi1 - x(j,1), delz = xi2 - x(j,2);
    KK_FLOAT rsqmin = delx*delx + dely*dely + delz*delz;
    while (d_sametag[j] >= 0) {
      j = d_sametag[j];
      delx = xi0 - x(j,0); dely = xi1 - x(j,1); delz = xi2 - x(j,2);
      const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < rsqmin) { rsqmin = rsq; closest = j; }
    }
    return closest;
  }

  // TIP4P device state
  typename AT::t_int_1d_randomread d_type;
  typename AT::t_tagint_1d_randomread d_tag;
  typename AT::t_int_1d d_sametag;
  typename AT::t_kkfloat_1d_3 d_xM;        // grid position per atom (M site for O)
  typename AT::t_int_1d d_iH1, d_iH2;      // the two H of each O atom

  int map_style;
  DAT::tdual_int_1d k_map_array;
  dual_hash_type k_map_hash;

  DAT::tdual_int_scalar k_tip4p_flag;      // set on device if an H is missing

  KK_FLOAT alpha_kk;
  int typeO_kk, typeH_kk;
  int nmax_tip4p;

  // base-class members referenced by the device kernels
  using Base::x; using Base::f; using Base::q;
  using Base::d_part2grid; using Base::d_rho1d;
  using Base::d_density_brick;
  using Base::d_vdx_brick; using Base::d_vdy_brick; using Base::d_vdz_brick;
  using Base::d_u_brick;
  using Base::d_v0_brick; using Base::d_v1_brick; using Base::d_v2_brick;
  using Base::d_v3_brick; using Base::d_v4_brick; using Base::d_v5_brick;
  using Base::d_eatom; using Base::d_vatom;
  using Base::k_flag;
  using Base::compute_rho1d;
  using Base::boxlo_kk; using Base::shift_kk; using Base::shiftone_kk;
  using Base::delxinv_kk; using Base::delyinv_kk; using Base::delzinv_kk;
  using Base::delvolinv_kk; using Base::qscale_kk;
  using Base::nlower; using Base::nupper; using Base::order;
  using Base::nxlo_out; using Base::nxhi_out; using Base::nylo_out;
  using Base::nyhi_out; using Base::nzlo_out; using Base::nzhi_out;
  using Base::ngrid; using Base::ix; using Base::iy; using Base::nlocal;
  using Base::slabflag; using Base::eflag_atom; using Base::vflag_atom;
  using Base::qsum; using Base::dipole_all; using Base::zprd_slab; using Base::ffact;
  using Base::numx_out; using Base::numy_out; using Base::numz_out;
};

}    // namespace LAMMPS_NS

#endif
#endif
