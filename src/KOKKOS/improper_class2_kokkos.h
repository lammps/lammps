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

#ifdef IMPROPER_CLASS
// clang-format off
ImproperStyle(class2/kk,ImproperClass2Kokkos<LMPDeviceType>);
ImproperStyle(class2/kk/device,ImproperClass2Kokkos<LMPDeviceType>);
ImproperStyle(class2/kk/host,ImproperClass2Kokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_IMPROPER_CLASS2_KOKKOS_H
#define LMP_IMPROPER_CLASS2_KOKKOS_H

#include "improper_class2.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int NEWTON_BOND, int EVFLAG>
struct TagImproperClass2Compute{};

template<int NEWTON_BOND, int EVFLAG>
struct TagImproperClass2AngleAngle{};

template<class DeviceType>
class ImproperClass2Kokkos : public ImproperClass2 {
 public:
  typedef DeviceType device_type;
  typedef EV_FLOAT value_type;
  typedef ArrayTypes<DeviceType> AT;

  ImproperClass2Kokkos(class LAMMPS *);
  ~ImproperClass2Kokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void read_restart(FILE *) override;

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagImproperClass2Compute<NEWTON_BOND,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagImproperClass2Compute<NEWTON_BOND,EVFLAG>, const int&) const;

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagImproperClass2AngleAngle<NEWTON_BOND,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagImproperClass2AngleAngle<NEWTON_BOND,EVFLAG>, const int&) const;

  //template<int NEWTON_BOND>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int i1, const int i2, const int i3, const int i4,
                          F_FLOAT &eimproper, F_FLOAT *f1, F_FLOAT *f3, F_FLOAT *f4,
                          const F_FLOAT &vb1x, const F_FLOAT &vb1y, const F_FLOAT &vb1z,
                          const F_FLOAT &vb2x, const F_FLOAT &vb2y, const F_FLOAT &vb2z,
                          const F_FLOAT &vb3x, const F_FLOAT &vb3y, const F_FLOAT &vb3z) const;

 protected:

  class NeighborKokkos *neighborKK;

  typename AT::t_x_array_randomread x;
  typename Kokkos::View<double*[3],typename AT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic> > f;
  typename AT::t_int_2d improperlist;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename AT::t_efloat_1d d_eatom;
  typename AT::t_virial_array d_vatom;

  int nlocal,newton_bond;
  int eflag,vflag;

  DAT::tdual_int_scalar k_warning_flag;
  typename AT::t_int_scalar d_warning_flag;
  HAT::t_int_scalar h_warning_flag;

  typename AT::tdual_ffloat_1d k_k0,k_chi0;
  typename AT::tdual_ffloat_1d k_aa_k1,k_aa_k2,k_aa_k3,k_aa_theta0_1,k_aa_theta0_2,k_aa_theta0_3;
  typename AT::tdual_ffloat_1d k_setflag_i,k_setflag_aa,k_setflag;

  typename AT::t_ffloat_1d d_k0,d_chi0;
  typename AT::t_ffloat_1d d_aa_k1,d_aa_k2,d_aa_k3,d_aa_theta0_1,d_aa_theta0_2,d_aa_theta0_3;
  typename AT::t_ffloat_1d d_setflag_i,d_setflag_aa,d_setflag;

  void allocate();
};

}

#endif
#endif

