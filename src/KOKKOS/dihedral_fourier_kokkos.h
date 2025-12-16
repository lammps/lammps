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
DihedralStyle(fourier/kk,DihedralFourierKokkos<LMPDeviceType>);
DihedralStyle(fourier/kk/device,DihedralFourierKokkos<LMPDeviceType>);
DihedralStyle(fourier/kk/host,DihedralFourierKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_DIHEDRAL_FOURIER_KOKKOS_H
#define LMP_DIHEDRAL_FOURIER_KOKKOS_H

#include "dihedral_fourier.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int NEWTON_BOND, int EVFLAG>
struct TagDihedralFourierCompute{};

template<class DeviceType>
class DihedralFourierKokkos : public DihedralFourier {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  DihedralFourierKokkos(class LAMMPS *);
  ~DihedralFourierKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void read_restart(FILE *) override;

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagDihedralFourierCompute<NEWTON_BOND,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagDihedralFourierCompute<NEWTON_BOND,EVFLAG>, const int&) const;

  //template<int NEWTON_BOND>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int i1, const int i2, const int i3, const int i4,
                          KK_FLOAT &edihedral, KK_FLOAT *f1, KK_FLOAT *f3, KK_FLOAT *f4,
                          const KK_FLOAT &vb1x, const KK_FLOAT &vb1y, const KK_FLOAT &vb1z,
                          const KK_FLOAT &vb2x, const KK_FLOAT &vb2y, const KK_FLOAT &vb2z,
                          const KK_FLOAT &vb3x, const KK_FLOAT &vb3y, const KK_FLOAT &vb3z) const;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;

 protected:
  int allocated_kokkos;

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

  DAT::tdual_kkfloat_2d k_k;
  DAT::tdual_kkfloat_2d k_cos_shift;
  DAT::tdual_kkfloat_2d k_sin_shift;
  DAT::tdual_int_2d k_multiplicity;
  DAT::tdual_int_1d k_nterms;

  typename AT::t_kkfloat_2d d_k;
  typename AT::t_kkfloat_2d d_cos_shift;
  typename AT::t_kkfloat_2d d_sin_shift;
  typename AT::t_int_2d d_multiplicity;
  typename AT::t_int_1d d_nterms;

  void allocate_kokkos();
};

}    // namespace LAMMPS_NS

#endif
#endif
