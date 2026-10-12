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
FixStyle(nvk/kk,FixNVKKokkos<LMPDeviceType>);
FixStyle(nvk/kk/device,FixNVKKokkos<LMPDeviceType>);
FixStyle(nvk/kk/host,FixNVKKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_NVK_KOKKOS_H
#define LMP_FIX_NVK_KOKKOS_H

#include "fix_nvk.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

// a and b of Minary 2003 eqs 4.12/4.13, summed in one pass

struct s_FixNVK_ab {
  double a, b;
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  s_FixNVK_ab() {
    a = b = 0.0;
  }
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  s_FixNVK_ab& operator+=(const s_FixNVK_ab &rhs) {
    a += rhs.a;
    b += rhs.b;
    return *this;
  }
};

// RMASS selects per-atom vs per-type mass; XUPDATE is 1 in the first half of
// the step (which also advances x) and 0 in the second

template<int RMASS> struct TagFixNVKReduce{};
template<int RMASS, int XUPDATE> struct TagFixNVKUpdate{};

template<class DeviceType>
class FixNVKKokkos : public FixNVK {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef s_FixNVK_ab value_type;

  FixNVKKokkos(class LAMMPS *, int, char **);
  void init() override;
  void initial_integrate(int) override;
  void final_integrate() override;

  template<int RMASS>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixNVKReduce<RMASS>, const int &, s_FixNVK_ab &) const;

  template<int RMASS, int XUPDATE>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixNVKUpdate<RMASS,XUPDATE>, const int &) const;

 private:
  void integrate(int);

  typename AT::t_kkfloat_1d_3_lr d_x;
  typename AT::t_kkfloat_1d_3 d_v;
  typename AT::t_kkacc_1d_3_randomread d_f;
  typename AT::t_kkfloat_1d_randomread d_rmass;
  typename AT::t_kkfloat_1d_randomread d_mass;
  typename AT::t_int_1d_randomread d_type;
  typename AT::t_int_1d_randomread d_mask;

  KK_FLOAT l_s, l_sdot, l_dtv, l_ftm2v;
};

}

#endif
#endif
