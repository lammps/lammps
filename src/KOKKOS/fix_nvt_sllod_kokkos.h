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

template<class DeviceType>
class FixNVTSllodKokkos : public FixNHKokkos<DeviceType> {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixNVTSllodKokkos(class LAMMPS *, int, char **);
  ~FixNVTSllodKokkos() {}
  void init();

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixNVTSllod_temp1, const int& i) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixNVTSllod_temp2, const int& i) const;

 private:
  int nondeformbias;

  void nh_v_temp();

 protected:
  typename AT::t_x_array x;
  typename AT::t_v_array v;
  typename AT::t_v_array vdelu;
  typename AT::t_f_array_const f;
  typename AT::t_float_1d rmass;
  typename AT::t_float_1d mass;
  typename AT::t_int_1d type;
  typename AT::t_int_1d mask;

  Few<double, 6> d_h_two;

  class DomainKokkos *domainKK;
  class AtomKokkos *atomKK;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Temperature control must be used with fix nvt/sllod

Self-explanatory.

E: Pressure control can not be used with fix nvt/sllod

Self-explanatory.

E: Temperature for fix nvt/sllod does not have a bias

The specified compute must compute temperature with a bias.

E: Using fix nvt/sllod with inconsistent fix deform remap option

Fix nvt/sllod requires that deforming atoms have a velocity profile
provided by "remap v" as a fix deform option.

E: Using fix nvt/sllod with no fix deform defined

Self-explanatory.

*/
