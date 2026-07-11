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
FixStyle(nvt/sllod/kk,FixNVTSllodKokkos<LMPDeviceType>);
FixStyle(nvt/sllod/kk/device,FixNVTSllodKokkos<LMPDeviceType>);
FixStyle(nvt/sllod/kk/host,FixNVTSllodKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_FIX_NVT_SLLOD_KOKKOS_H
#define LMP_FIX_NVT_SLLOD_KOKKOS_H

#include "fix_nh_kokkos.h"
#include "kokkos_few.h"
#include "kokkos_type.h"

// clang-format off
namespace LAMMPS_NS {

struct TagFixNVTSllod_temp1{};
struct TagFixNVTSllod_temp2{};

template<bool PSLLOD>
struct TagFixNVTSllod_nvex{};

template<class DeviceType>
class FixNVTSllodKokkos : public FixNHKokkos<DeviceType> {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixNVTSllodKokkos(class LAMMPS *, int, char **);

  void init() override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixNVTSllod_temp1, const int& i) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixNVTSllod_temp2, const int& i) const;

  template<bool PSLLOD>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixNVTSllod_nvex<PSLLOD>, const int& i) const;

 private:
  int nondeformbias;
  int psllod_flag;
  int peculiar_flag;  // 0 for lab frame, 1 for peculiar
  int kick_flag;      // 0 for no initial velocity kick, 1 for kick
  enum {REVERSIBLE, LEGACY} integrator;

  void nh_v_temp() override;
  void nve_x() override;
  int size_restart_global() override;
  int pack_restart_data(double *list) override;
  void restart(char *buf) override;
  int modify_param(int narg, char **arg) override;

 protected:
  typename AT::t_kkfloat_1d_3_lr x;
  typename AT::t_kkfloat_1d_3 v;
  typename AT::t_kkfloat_1d_3 vdelu;
  typename AT::t_kkacc_1d_3_const f;
  typename AT::t_kkfloat_1d rmass;
  typename AT::t_kkfloat_1d mass;
  typename AT::t_int_1d type;
  typename AT::t_int_1d mask;

  Few<double, 6> d_h_two;
  Few<double, 3> d_xfac;
  Few<double, 3> d_vfac;
  Few<double, 3> d_xmid;
  Few<double, 5> d_xlo;

  class DomainKokkos *domainKK;
  class AtomKokkos *atomKK;
};

}    // namespace LAMMPS_NS

#endif
#endif

