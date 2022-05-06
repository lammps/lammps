/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef NBIN_CLASS
// clang-format off
NBinStyle(kk/host,
          NBinKokkos<LMPHostType>,
          NB_STANDARD | NB_KOKKOS_HOST);

NBinStyle(kk/device,
          NBinKokkos<LMPDeviceType>,
          NB_STANDARD | NB_KOKKOS_DEVICE);
// clang-format on
#else

// clang-format off
#ifndef LMP_NBIN_KOKKOS_H
#define LMP_NBIN_KOKKOS_H

#include "nbin_standard.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class NBinKokkos : public NBinStandard {
 public:
  typedef ArrayTypes<DeviceType> AT;

  NBinKokkos(class LAMMPS *);

  void bin_atoms_setup(int) override;
  void bin_atoms() override;

  int atoms_per_bin;
  DAT::tdual_int_1d k_bincount;
  DAT::tdual_int_2d k_bins;
  DAT::tdual_int_1d k_atom2bin;

  typename AT::t_int_1d bincount;
  const typename AT::t_int_1d_const c_bincount;
  typename AT::t_int_2d bins;
  typename AT::t_int_2d_const c_bins;
  typename AT::t_int_1d atom2bin;
  typename AT::t_int_scalar d_resize;
  typename ArrayTypes<LMPHostType>::t_int_scalar h_resize;
  typename AT::t_x_array_randomread x;

  KOKKOS_INLINE_FUNCTION
  void binatomsItem(const int &i) const;

  KOKKOS_INLINE_FUNCTION
  int coord2bin(const X_FLOAT & x,const X_FLOAT & y,const X_FLOAT & z) const
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

    return (iz-mbinzlo)*mbiny*mbinx + (iy-mbinylo)*mbinx + (ix-mbinxlo);
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
};

template<class DeviceType>
struct NPairKokkosBinAtomsFunctor {
  typedef DeviceType device_type;

  const NBinKokkos<DeviceType> c;

  NPairKokkosBinAtomsFunctor(const NBinKokkos<DeviceType> &_c):
    c(_c) {};
  ~NPairKokkosBinAtomsFunctor() {}
  KOKKOS_INLINE_FUNCTION
  void operator() (const int & i) const {
    c.binatomsItem(i);
  }
};

}

#endif
#endif

