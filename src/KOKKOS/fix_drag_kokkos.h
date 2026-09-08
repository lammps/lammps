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
FixStyle(drag/kk,FixDragKokkos<LMPDeviceType>);
FixStyle(drag/kk/device,FixDragKokkos<LMPDeviceType>);
FixStyle(drag/kk/host,FixDragKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_DRAG_KOKKOS_H
#define LMP_FIX_DRAG_KOKKOS_H

#include "fix_drag.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixDrag{};

template<class DeviceType>
class FixDragKokkos : public FixDrag {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef double value_type[];
  const int value_count = 3;

  FixDragKokkos(class LAMMPS *, int, char **);
  ~FixDragKokkos() override;
  void init() override;
  void post_force(int) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixDrag, const int &, value_type) const;

 private:
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread mask;

  // domain data stored as members for device access
  int m_triclinic;
  int m_xperiodic, m_yperiodic, m_zperiodic;
  KK_FLOAT m_xprd, m_yprd, m_zprd;
  KK_FLOAT m_xprd_half, m_yprd_half, m_zprd_half;
  KK_FLOAT m_xy, m_xz, m_yz;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void minimum_image(KK_FLOAT &dx, KK_FLOAT &dy, KK_FLOAT &dz) const;

};

}

#endif
#endif
