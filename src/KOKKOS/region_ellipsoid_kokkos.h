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
RegionStyle(ellipsoid/kk,RegEllipsoidKokkos<LMPDeviceType>);
RegionStyle(ellipsoid/kk/device,RegEllipsoidKokkos<LMPDeviceType>);
RegionStyle(ellipsoid/kk/host,RegEllipsoidKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_REGION_ELLIPSOID_KOKKOS_H
#define LMP_REGION_ELLIPSOID_KOKKOS_H

#include "region_ellipsoid.h"

#include "kokkos_base.h"
#include "kokkos_type.h"
#include "region_remap_kokkos.h"

namespace LAMMPS_NS {

struct TagRegEllipsoidMatchAll{};

// device-side region matching only: fix wall/region/kk additionally needs
// the surface contact machinery (cf. region block/kk), which is not ported

template<class DeviceType>
class RegEllipsoidKokkos : public RegEllipsoid, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  RegEllipsoidKokkos(class LAMMPS *, int, char **);
  void init() override;

  void match_all_kokkos(int, DAT::tdual_int_1d) override;

  void prematch() override
  {
    Region::prematch();
    boxremap.capture(domain);
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagRegEllipsoidMatchAll, const int&) const;

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
  int dimension;             // cached domain->dimension: a device kernel
                             // cannot dereference the Domain pointer
  RegionRemapKokkos boxremap;
  typename AT::t_int_1d d_match;
  typename AT::t_kkfloat_1d_3_lr_randomread d_x;
  typename AT::t_int_1d_randomread d_mask;


// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int k_inside(double x, double y, double z) const
  {
    if (dimension == 3) {
      double delx = b * c * (x - xc);
      double dely = a * c * (y - yc);
      double delz = a * b * (z - zc);
      double r = delx * delx + dely * dely + delz * delz;
      double rc = a * a * b * b * c * c;
      if (r <= rc) return 1;
    } else {
      double delx = b * (x - xc);
      double dely = a * (y - yc);
      double r = delx * delx + dely * dely;
      double rc = a * a * b * b;
      if (r <= rc) return 1;
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
    double va[3],vb[3],vc[3],vd[3],disp[3];

    double sine = sin(angle);
    double cosine = cos(angle);
    vd[0] = x - point[0];
    vd[1] = y - point[1];
    vd[2] = z - point[2];
    double x0dotr = vd[0]*runit[0] + vd[1]*runit[1] + vd[2]*runit[2];
    vc[0] = x0dotr * runit[0];
    vc[1] = x0dotr * runit[1];
    vc[2] = x0dotr * runit[2];
    va[0] = vd[0] - vc[0];
    va[1] = vd[1] - vc[1];
    va[2] = vd[2] - vc[2];
    vb[0] = runit[1]*va[2] - runit[2]*va[1];
    vb[1] = runit[2]*va[0] - runit[0]*va[2];
    vb[2] = runit[0]*va[1] - runit[1]*va[0];
    disp[0] = va[0]*cosine  + vb[0]*sine;
    disp[1] = va[1]*cosine  + vb[1]*sine;
    disp[2] = va[2]*cosine  + vb[2]*sine;
    x = point[0] + vc[0] + disp[0];
    y = point[1] + vc[1] + disp[1];
    z = point[2] + vc[2] + disp[2];
  }

};

}

#endif
#endif
