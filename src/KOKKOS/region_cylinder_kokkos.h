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

#ifdef REGION_CLASS
// clang-format off
RegionStyle(cylinder/kk,RegCylinderKokkos<LMPDeviceType>);
RegionStyle(cylinder/kk/device,RegCylinderKokkos<LMPDeviceType>);
RegionStyle(cylinder/kk/host,RegCylinderKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_REGION_CYLINDER_KOKKOS_H
#define LMP_REGION_CYLINDER_KOKKOS_H

#include "region_cylinder.h"

#include "kokkos_base.h"
#include "kokkos_type.h"
#include "region_remap_kokkos.h"

namespace LAMMPS_NS {

struct TagRegCylinderMatchAll{};

// device-side region matching only: fix wall/region/kk additionally needs
// the surface contact machinery (cf. region block/kk), which is not ported

template<class DeviceType>
class RegCylinderKokkos : public RegCylinder, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  RegCylinderKokkos(class LAMMPS *, int, char **);

  void match_all_kokkos(int, DAT::tdual_int_1d) override;

  void prematch() override
  {
    Region::prematch();
    boxremap.capture(domain);
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagRegCylinderMatchAll, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int match_kokkos(double x, double y, double z) const
  {
    boxremap.remap(x,y,z);
    if (dynamic) inverse_transform(x,y,z);
    if (openflag) return 1;
    return !(k_inside(x,y,z) ^ interior);
  }

 private:
  int groupbit;
  RegionRemapKokkos boxremap;
  typename AT::t_int_1d d_match;
  typename AT::t_kkfloat_1d_3_lr_randomread d_x;
  typename AT::t_int_1d_randomread d_mask;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int k_inside(double x, double y, double z) const
  {
    double del1, del2, dist;

    if (axis == 'x') {
      del1 = y - c1;
      del2 = z - c2;
      dist = sqrt(del1 * del1 + del2 * del2);
      if (dist <= radius && x >= lo && x <= hi) return 1;
    } else if (axis == 'y') {
      del1 = x - c1;
      del2 = z - c2;
      dist = sqrt(del1 * del1 + del2 * del2);
      if (dist <= radius && y >= lo && y <= hi) return 1;
    } else {
      del1 = x - c1;
      del2 = y - c2;
      dist = sqrt(del1 * del1 + del2 * del2);
      if (dist <= radius && z >= lo && z <= hi) return 1;
    }
    return 0;
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void inverse_transform(double &x, double &y, double &z) const
  {
    if (moveflag) {
      x -= dx;
      y -= dy;
      z -= dz;
    }
    if (rotateflag) rotate(x,y,z,-theta);
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void rotate(double &x, double &y, double &z, double angle) const
  {
    double a[3],b[3],c[3],d[3],disp[3];

    double sine = sin(angle);
    double cosine = cos(angle);
    d[0] = x - point[0];
    d[1] = y - point[1];
    d[2] = z - point[2];
    double x0dotr = d[0]*runit[0] + d[1]*runit[1] + d[2]*runit[2];
    c[0] = x0dotr * runit[0];
    c[1] = x0dotr * runit[1];
    c[2] = x0dotr * runit[2];
    a[0] = d[0] - c[0];
    a[1] = d[1] - c[1];
    a[2] = d[2] - c[2];
    b[0] = runit[1]*a[2] - runit[2]*a[1];
    b[1] = runit[2]*a[0] - runit[0]*a[2];
    b[2] = runit[0]*a[1] - runit[1]*a[0];
    disp[0] = a[0]*cosine  + b[0]*sine;
    disp[1] = a[1]*cosine  + b[1]*sine;
    disp[2] = a[2]*cosine  + b[2]*sine;
    x = point[0] + c[0] + disp[0];
    y = point[1] + c[1] + disp[1];
    z = point[2] + c[2] + disp[2];
  }

};

}

#endif
#endif
