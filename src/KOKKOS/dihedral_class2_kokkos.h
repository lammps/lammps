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

#ifdef DIHEDRAL_CLASS
// clang-format off
DihedralStyle(class2/kk,DihedralClass2Kokkos<LMPDeviceType>);
DihedralStyle(class2/kk/device,DihedralClass2Kokkos<LMPDeviceType>);
DihedralStyle(class2/kk/host,DihedralClass2Kokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_DIHEDRAL_CLASS2_KOKKOS_H
#define LMP_DIHEDRAL_CLASS2_KOKKOS_H

#include "dihedral_class2.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int NEWTON_BOND, int EVFLAG>
struct TagDihedralClass2Compute{};

template<class DeviceType>
class DihedralClass2Kokkos : public DihedralClass2 {
 public:
  typedef DeviceType device_type;
  typedef EV_FLOAT value_type;
  typedef ArrayTypes<DeviceType> AT;

  DihedralClass2Kokkos(class LAMMPS *);
  virtual ~DihedralClass2Kokkos();
  void compute(int, int);
  void coeff(int, char **);
  void read_restart(FILE *);

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagDihedralClass2Compute<NEWTON_BOND,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagDihedralClass2Compute<NEWTON_BOND,EVFLAG>, const int&) const;

  //template<int NEWTON_BOND>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int i1, const int i2, const int i3, const int i4,
                          F_FLOAT &edihedral, F_FLOAT *f1, F_FLOAT *f3, F_FLOAT *f4,
                          const F_FLOAT &vb1x, const F_FLOAT &vb1y, const F_FLOAT &vb1z,
                          const F_FLOAT &vb2x, const F_FLOAT &vb2y, const F_FLOAT &vb2z,
                          const F_FLOAT &vb3x, const F_FLOAT &vb3y, const F_FLOAT &vb3z) const;

 protected:

  class NeighborKokkos *neighborKK;

  typename AT::t_x_array_randomread x;
  typename AT::t_f_array f;
  typename AT::t_int_2d dihedrallist;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename AT::t_efloat_1d d_eatom;
  typename AT::t_virial_array d_vatom;

  int nlocal,newton_bond;
  int eflag,vflag;

  DAT::tdual_int_scalar k_warning_flag;
  typename AT::t_int_scalar d_warning_flag;
  HAT::t_int_scalar h_warning_flag;

  typename AT::tdual_ffloat_1d k_k2, k_k3, k_k1;
  typename AT::tdual_ffloat_1d k_phi1, k_phi2, k_phi3;
  typename AT::tdual_ffloat_1d k_mbt_f1, k_mbt_f2, k_mbt_f3, k_mbt_r0;
  typename AT::tdual_ffloat_1d k_ebt_f1_1, k_ebt_f2_1, k_ebt_f3_1, k_ebt_r0_1;
  typename AT::tdual_ffloat_1d k_ebt_f1_2, k_ebt_f2_2, k_ebt_f3_2, k_ebt_r0_2;
  typename AT::tdual_ffloat_1d k_at_f1_1, k_at_f2_1, k_at_f3_1, k_at_theta0_1;
  typename AT::tdual_ffloat_1d k_at_f1_2, k_at_f2_2, k_at_f3_2, k_at_theta0_2;
  typename AT::tdual_ffloat_1d k_aat_k, k_aat_theta0_1, k_aat_theta0_2;
  typename AT::tdual_ffloat_1d k_bb13t_k, k_bb13t_r10, k_bb13t_r30;
  typename AT::tdual_ffloat_1d k_setflag_d, k_setflag_mbt, k_setflag_ebt;
  typename AT::tdual_ffloat_1d k_setflag_at, k_setflag_aat, k_setflag_bb13t;

  typename AT::t_ffloat_1d d_k2, d_k3, d_k1;
  typename AT::t_ffloat_1d d_phi1, d_phi2, d_phi3;
  typename AT::t_ffloat_1d d_mbt_f1, d_mbt_f2, d_mbt_f3, d_mbt_r0;
  typename AT::t_ffloat_1d d_ebt_f1_1, d_ebt_f2_1, d_ebt_f3_1, d_ebt_r0_1;
  typename AT::t_ffloat_1d d_ebt_f1_2, d_ebt_f2_2, d_ebt_f3_2, d_ebt_r0_2;
  typename AT::t_ffloat_1d d_at_f1_1, d_at_f2_1, d_at_f3_1, d_at_theta0_1;
  typename AT::t_ffloat_1d d_at_f1_2, d_at_f2_2, d_at_f3_2, d_at_theta0_2;
  typename AT::t_ffloat_1d d_aat_k, d_aat_theta0_1, d_aat_theta0_2;
  typename AT::t_ffloat_1d d_bb13t_k, d_bb13t_r10, d_bb13t_r30;
  typename AT::t_ffloat_1d d_setflag_d, d_setflag_mbt, d_setflag_ebt;
  typename AT::t_ffloat_1d d_setflag_at, d_setflag_aat, d_setflag_bb13t;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

W: Dihedral problem

Conformation of the 4 listed dihedral atoms is extreme; you may want
to check your simulation geometry.

*/
