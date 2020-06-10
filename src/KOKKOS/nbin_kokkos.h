/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef NBIN_CLASS

NBinStyle(kk/host,
          NBinKokkos<Host>,
          NB_KOKKOS_HOST)

NBinStyle(kk/device,
          NBinKokkos<Device>,
          NB_KOKKOS_DEVICE)

#else

#ifndef LMP_NBIN_KOKKOS_H
#define LMP_NBIN_KOKKOS_H

#include "nbin_standard.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<ExecutionSpace Space>
class NBinKokkos : public NBinStandard {
 public:
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef ArrayTypes<Space> AT;

  NBinKokkos(class LAMMPS *);
  ~NBinKokkos() {}
  void bin_atoms_setup(int);
  void bin_atoms();

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
  HAT::t_int_scalar h_resize;
  typename AT::t_float_1d_3_lr_randomread x;

  KOKKOS_INLINE_FUNCTION
  void binatomsItem(const int &i) const;

  KOKKOS_INLINE_FUNCTION
  int coord2bin(const KK_FLOAT & x,const KK_FLOAT & y,const KK_FLOAT & z) const
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
  int coord2bin(const KK_FLOAT & x,const KK_FLOAT & y,const KK_FLOAT & z, int* i) const
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
  KK_FLOAT bboxlo_[3],bboxhi_[3];
};

template<ExecutionSpace Space>
struct NPairKokkosBinAtomsFunctor {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  const NBinKokkos<Space> c;

  NPairKokkosBinAtomsFunctor(const NBinKokkos<Space> &_c):
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

/* ERROR/WARNING messages:

*/
