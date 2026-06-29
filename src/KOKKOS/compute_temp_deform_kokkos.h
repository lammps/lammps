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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(temp/deform/kk,ComputeTempDeformKokkos<LMPDeviceType>);
ComputeStyle(temp/deform/kk/device,ComputeTempDeformKokkos<LMPDeviceType>);
ComputeStyle(temp/deform/kk/host,ComputeTempDeformKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_TEMP_DEFORM_KOKKOS_H
#define LMP_COMPUTE_TEMP_DEFORM_KOKKOS_H

#include "compute_temp_deform.h"
#include "kokkos_few.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int RMASS>
struct TagComputeTempDeformScalar{};

template<int RMASS>
struct TagComputeTempDeformVector{};

struct TagComputeTempDeformRemoveBias{};

struct TagComputeTempDeformRestoreBias{};

struct TagComputeTempDeformApplyBias{};

template<class DeviceType>
class ComputeTempDeformKokkos: public ComputeTempDeform {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  ComputeTempDeformKokkos(class LAMMPS *, int, char **);
  void post_constructor() override;
  double compute_scalar() override;
  void compute_vector() override;
  void remove_bias_all() override;
  void remove_bias_all_kk() override;
  void restore_bias_all() override;
  void restore_bias_all_kk() override;

  void remove_deform_bias_all() override;
  void remove_deform_bias_all_kk() override;
  void restore_deform_bias_all() override;
  void restore_deform_bias_all_kk() override;
  void apply_deform_bias_all(double dtv = 0.0) override;
  void apply_deform_bias_all_kk(double dtv = 0.0) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempDeformRemoveBias, const int &i) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempDeformRestoreBias, const int &i) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTempDeformApplyBias, const int &i) const;

 protected:
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_3 v;
  typename AT::t_kkfloat_1d_3 vbiasall;
  typename AT::t_kkfloat_1d_randomread rmass;
  typename AT::t_kkfloat_1d_randomread mass;
  typename AT::t_int_1d_randomread type;
  typename AT::t_int_1d_randomread mask;

  class DomainKokkos *domainKK;

  Few<double, 6> h_rate, h_ratelo, d_grad_u;
  Few<double, 5> d_xref;  // stores xmid[3], ylo, zlo

};

}    // namespace LAMMPS_NS

#endif
#endif
