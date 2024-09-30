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

#ifdef FIX_CLASS
// clang-format off
FixStyle(efield/kk,FixEfieldKokkos<LMPDeviceType>);
FixStyle(efield/kk/device,FixEfieldKokkos<LMPDeviceType>);
FixStyle(efield/kk/host,FixEfieldKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_EFIELD_KOKKOS_H
#define LMP_FIX_EFIELD_KOKKOS_H

#include "fix_efield.h"
#include "kokkos_type.h"
#include "kokkos_few.h"

namespace LAMMPS_NS {

template<int QFLAG, int MUFLAG>
struct TagFixEfieldConstant{};

template<int QFLAG, int MUFLAG>
struct TagFixEfieldNonConstant{};

template<class DeviceType>
class FixEfieldKokkos : public FixEfield {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixEfieldKokkos(class LAMMPS *, int, char **);
  ~FixEfieldKokkos() override;
  void init() override;
  void post_force(int) override;

  typedef double value_type[];
  const int value_count = 10;

  template<int QFLAG, int MUFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixEfieldConstant<QFLAG,MUFLAG>, const int&, value_type) const;

  template<int QFLAG, int MUFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixEfieldNonConstant<QFLAG,MUFLAG>, const int&, value_type) const;

 private:

  DAT::tdual_ffloat_2d k_efield;
  typename AT::t_ffloat_2d_randomread d_efield;
  typename AT::t_int_1d d_match;

  typename AT::t_x_array_randomread d_x;
  typename AT::t_float_1d_randomread d_q;
  typename AT::t_mu_array_randomread d_mu;
  typename AT::t_f_array d_f;
  typename AT::t_f_array d_torque;
  typename AT::t_imageint_1d_randomread d_image;
  typename AT::t_int_1d_randomread d_mask;

  Few<double,3> prd;
  Few<double,6> h;
  int triclinic;

  DAT::tdual_virial_array k_vatom;
  typename AT::t_virial_array d_vatom;

  KOKKOS_INLINE_FUNCTION
  void v_tally(value_type, int, double*) const;
};

}

#endif
#endif

