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

#ifdef IMPROPER_CLASS
// clang-format off
ImproperStyle(harmonic/kk,ImproperHarmonicKokkos<LMPDeviceType>);
ImproperStyle(harmonic/kk/device,ImproperHarmonicKokkos<LMPDeviceType>);
ImproperStyle(harmonic/kk/host,ImproperHarmonicKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_IMPROPER_HARMONIC_KOKKOS_H
#define LMP_IMPROPER_HARMONIC_KOKKOS_H

#include "improper_harmonic.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int NEWTON_BOND, int EVFLAG>
struct TagImproperHarmonicCompute{};

template<class DeviceType>
class ImproperHarmonicKokkos : public ImproperHarmonic {
 public:
  typedef DeviceType device_type;
  typedef EV_FLOAT value_type;
  typedef ArrayTypes<DeviceType> AT;

  ImproperHarmonicKokkos(class LAMMPS *);
  virtual ~ImproperHarmonicKokkos();
  void compute(int, int);
  void coeff(int, char **);
  void read_restart(FILE *);

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagImproperHarmonicCompute<NEWTON_BOND,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagImproperHarmonicCompute<NEWTON_BOND,EVFLAG>, const int&) const;

  //template<int NEWTON_BOND>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int i1, const int i2, const int i3, const int i4,
                          F_FLOAT &eimproper, F_FLOAT *f1, F_FLOAT *f3, F_FLOAT *f4,
                          const F_FLOAT &vb1x, const F_FLOAT &vb1y, const F_FLOAT &vb1z,
                          const F_FLOAT &vb2x, const F_FLOAT &vb2y, const F_FLOAT &vb2z,
                          const F_FLOAT &vb3x, const F_FLOAT &vb3y, const F_FLOAT &vb3z) const;

 protected:

  class NeighborKokkos *neighborKK;

  typedef typename KKDevice<DeviceType>::value KKDeviceType;
  typename AT::t_x_array_randomread x;
  typename Kokkos::View<double*[3],typename AT::t_f_array::array_layout,KKDeviceType,Kokkos::MemoryTraits<Kokkos::Atomic> > f;
  typename AT::t_int_2d improperlist;

  Kokkos::DualView<E_FLOAT*,Kokkos::LayoutRight,KKDeviceType> k_eatom;
  Kokkos::DualView<F_FLOAT*[6],Kokkos::LayoutRight,KKDeviceType> k_vatom;
  Kokkos::View<E_FLOAT*,Kokkos::LayoutRight,KKDeviceType,Kokkos::MemoryTraits<Kokkos::Atomic> > d_eatom;
  Kokkos::View<F_FLOAT*[6],Kokkos::LayoutRight,KKDeviceType,Kokkos::MemoryTraits<Kokkos::Atomic> > d_vatom;

  int nlocal,newton_bond;
  int eflag,vflag;

  Kokkos::DualView<int,DeviceType> k_warning_flag;
  typename Kokkos::DualView<int,DeviceType>::t_dev d_warning_flag;
  typename Kokkos::DualView<int,DeviceType>::t_host h_warning_flag;

  Kokkos::DualView<F_FLOAT*,DeviceType> k_k;
  Kokkos::DualView<F_FLOAT*,DeviceType> k_chi;

  typename Kokkos::DualView<F_FLOAT*,DeviceType>::t_dev d_k;
  typename Kokkos::DualView<F_FLOAT*,DeviceType>::t_dev d_chi;

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
