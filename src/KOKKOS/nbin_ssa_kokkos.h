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

#ifdef NBIN_CLASS
// clang-format off
NBinStyle(ssa/kk/host,
          NBinSSAKokkos<LMPHostType>,
          NB_STANDARD | NB_SSA | NB_KOKKOS_HOST);

NBinStyle(ssa/kk/device,
          NBinSSAKokkos<LMPDeviceType>,
          NB_STANDARD | NB_SSA | NB_KOKKOS_DEVICE);
// clang-format on
#else

// clang-format off
#ifndef LMP_NBIN_SSA_KOKKOS_H
#define LMP_NBIN_SSA_KOKKOS_H

#include "nbin_standard.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class NBinSSAKokkos : public NBinStandard {
 public:
  typedef ArrayTypes<DeviceType> AT;

  NBinSSAKokkos(class LAMMPS *);

  void bin_atoms_setup(int) override;
  void bin_atoms() override;

   // temporary array to hold the binID for each atom
  DAT::tdual_int_1d k_binID;
  typename AT::t_int_1d binID;
  typename AT::t_int_1d_const c_binID;

  int atoms_per_bin;
  DAT::tdual_int_1d k_bincount;
  DAT::tdual_int_2d k_bins;
  typename AT::t_int_1d bincount;
  typename AT::t_int_2d bins;
  typename AT::t_int_2d_const c_bins;

  int ghosts_per_gbin;
  DAT::tdual_int_1d k_gbincount;
  DAT::tdual_int_2d k_gbins;
  typename AT::t_int_1d gbincount;
  typename AT::t_int_2d gbins;
  typename AT::t_int_2d_const c_gbins;

  typename AT::t_int_scalar d_resize;
  typename ArrayTypes<LMPHostType>::t_int_scalar h_resize;
  typename AT::t_x_array_randomread x;

  // Bounds of the local atoms in the bins array
  typename AT::t_int_scalar d_lbinxlo;  // lowest local bin x-dim coordinate
  typename AT::t_int_scalar d_lbinylo;  // lowest local bin y-dim coordinate
  typename AT::t_int_scalar d_lbinzlo;  // lowest local bin z-dim coordinate
  typename AT::t_int_scalar d_lbinxhi;  // highest local bin x-dim coordinate
  typename AT::t_int_scalar d_lbinyhi;  // highest local bin y-dim coordinate
  typename AT::t_int_scalar d_lbinzhi;  // highest local bin z-dim coordinate
  typename ArrayTypes<LMPHostType>::t_int_scalar h_lbinxlo;
  typename ArrayTypes<LMPHostType>::t_int_scalar h_lbinylo;
  typename ArrayTypes<LMPHostType>::t_int_scalar h_lbinzlo;
  typename ArrayTypes<LMPHostType>::t_int_scalar h_lbinxhi;
  typename ArrayTypes<LMPHostType>::t_int_scalar h_lbinyhi;
  typename ArrayTypes<LMPHostType>::t_int_scalar h_lbinzhi;


  KOKKOS_INLINE_FUNCTION
  void binAtomsItem(const int &i) const;

  KOKKOS_INLINE_FUNCTION
  void binIDAtomsItem(const int &i, int &update) const;

  KOKKOS_INLINE_FUNCTION
  void binIDGhostsItem(const int &i, int &update) const;

  static KOKKOS_INLINE_FUNCTION
  void sortBin(
      typename AT::t_int_1d gbincount,
      typename AT::t_int_2d gbins,
      const int &ibin);

/* ----------------------------------------------------------------------
   convert atom coords into the ssa active interaction region number
------------------------------------------------------------------------- */
  KOKKOS_INLINE_FUNCTION
  int coord2ssaAIR(const X_FLOAT & x,const X_FLOAT & y,const X_FLOAT & z) const
  {
    int ix, iy, iz;
    ix = iy = iz = 0;
    if (z < sublo_[2]) iz = -1;
    if (z >= subhi_[2]) iz = 1;
    if (y < sublo_[1]) iy = -1;
    if (y >= subhi_[1]) iy = 1;
    if (x < sublo_[0]) ix = -1;
    if (x >= subhi_[0]) ix = 1;
    if (iz < 0) {
      return -1;
    } else if (iz == 0) {
      if (iy<0) return -1; // bottom left/middle/right
      if ((iy==0) && (ix<0) ) return -1; // left atoms
      if ((iy==0) && (ix==0)) return 0; // Locally owned atoms
      if ((iy==0) && (ix>0) ) return 2; // Right atoms
      if ((iy>0)  && (ix==0)) return 1; // Top-middle atoms
      if ((iy>0)  && (ix!=0)) return 3; // Top-right and top-left atoms
    } else { // iz > 0
      if ((ix==0) && (iy==0)) return 4; // Back atoms
      if ((ix==0) && (iy!=0)) return 5; // Top-back and bottom-back atoms
      if ((ix!=0) && (iy==0)) return 6; // Left-back and right-back atoms
      if ((ix!=0) && (iy!=0)) return 7; // Back corner atoms
    }
    return -2;
  }

  KOKKOS_INLINE_FUNCTION
  int coord2bin(const X_FLOAT & x,const X_FLOAT & y,const X_FLOAT & z, int* i) const
  {
    int ix,iy,iz;

    if (x >= bboxhi_[0])
      ix = static_cast<int> ((x-bboxhi_[0])*bininvx) + nbinx;
    else if (x >= bboxlo_[0]) {
      ix = static_cast<int> ((x-bboxlo_[0])*bininvx);
      ix = MIN(ix,nbinx-1);
    } else
      ix = static_cast<int> ((x-bboxlo_[0])*bininvx) - 1;

    if (y >= bboxhi_[1])
      iy = static_cast<int> ((y-bboxhi_[1])*bininvy) + nbiny;
    else if (y >= bboxlo_[1]) {
      iy = static_cast<int> ((y-bboxlo_[1])*bininvy);
      iy = MIN(iy,nbiny-1);
    } else
      iy = static_cast<int> ((y-bboxlo_[1])*bininvy) - 1;

    if (z >= bboxhi_[2])
      iz = static_cast<int> ((z-bboxhi_[2])*bininvz) + nbinz;
    else if (z >= bboxlo_[2]) {
      iz = static_cast<int> ((z-bboxlo_[2])*bininvz);
      iz = MIN(iz,nbinz-1);
    } else
      iz = static_cast<int> ((z-bboxlo_[2])*bininvz) - 1;

    i[0] = ix - mbinxlo;
    i[1] = iy - mbinylo;
    i[2] = iz - mbinzlo;

    return (iz-mbinzlo)*mbiny*mbinx + (iy-mbinylo)*mbinx + (ix-mbinxlo);
  }

 private:
  double bboxlo_[3],bboxhi_[3];
  double sublo_[3], subhi_[3];
};

template<class DeviceType>
struct NPairSSAKokkosBinAtomsFunctor {
  typedef DeviceType device_type;

  const NBinSSAKokkos<DeviceType> c;

  NPairSSAKokkosBinAtomsFunctor(const NBinSSAKokkos<DeviceType> &_c):
    c(_c) {};
  ~NPairSSAKokkosBinAtomsFunctor() {}
  KOKKOS_INLINE_FUNCTION
  void operator() (const int & i) const {
    c.binAtomsItem(i);
  }
};

template<class DeviceType>
struct NPairSSAKokkosBinIDAtomsFunctor {
  typedef DeviceType device_type;
  typedef int value_type;

  const NBinSSAKokkos<DeviceType> c;

  NPairSSAKokkosBinIDAtomsFunctor(const NBinSSAKokkos<DeviceType> &_c):
    c(_c) {};
  ~NPairSSAKokkosBinIDAtomsFunctor() {}
  KOKKOS_INLINE_FUNCTION
  void operator() (const int & i, value_type& update) const {
    c.binIDAtomsItem(i, update);
  }

  KOKKOS_INLINE_FUNCTION
  void join (value_type& dst,
             const value_type& src) const {
    if (dst < src) dst = src;
  }

  KOKKOS_INLINE_FUNCTION
  void init (value_type& dst) const {
    dst = INT_MIN;
  }
};

template<class DeviceType>
struct NPairSSAKokkosBinIDGhostsFunctor {
  typedef DeviceType device_type;
  typedef int value_type;

  const NBinSSAKokkos<DeviceType> c;

  NPairSSAKokkosBinIDGhostsFunctor(const NBinSSAKokkos<DeviceType> &_c):
    c(_c) {};
  ~NPairSSAKokkosBinIDGhostsFunctor() {}
  KOKKOS_INLINE_FUNCTION
  void operator() (const int & i, value_type& update) const {
    c.binIDGhostsItem(i, update);
  }

  KOKKOS_INLINE_FUNCTION
  void join (value_type& dst,
             const value_type& src) const {
    if (dst < src) dst = src;
  }

  KOKKOS_INLINE_FUNCTION
  void init (value_type& dst) const {
    dst = INT_MIN;
  }
};

}

#endif
#endif

