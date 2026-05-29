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
RegionStyle(sphere/kk,RegSphereKokkos<LMPDeviceType>);
RegionStyle(sphere/kk/device,RegSphereKokkos<LMPDeviceType>);
RegionStyle(sphere/kk/host,RegSphereKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_REGION_SPHERE_KOKKOS_H
#define LMP_REGION_SPHERE_KOKKOS_H

#include "region_sphere.h"

#include "kokkos_base.h"
#include "kokkos_type.h"
#include "region_block_kokkos.h"

namespace LAMMPS_NS {

struct TagRegSphereMatchAll{};

template<class DeviceType>
class RegSphereKokkos : public RegSphere, public KokkosBase  {
  friend class FixPour;

 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  RegSphereKokkos(class LAMMPS *, int, char **);
  ~RegSphereKokkos() override;

  void match_all_kokkos(int, DAT::tdual_int_1d) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagRegSphereMatchAll, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int match_kokkos(double x, double y, double z) const
  {
    if (dynamic) inverse_transform(x,y,z);
    if (openflag) return 1;
    return !(k_inside(x,y,z) ^ interior);
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int surface_kokkos(double x, double y, double z, double cutoff, contacts_t& contacts)
  {
    int ncontact;
    double xs, ys, zs;
    double xnear[3], xorig[3];

    xorig[0] = x; xorig[1] = y; xorig[2] = z;
    if (dynamic)
      inverse_transform(x, y, z);

    xnear[0] = x; xnear[1] = y; xnear[2] = z;

    if (!openflag) {
      if (interior) ncontact = surface_interior_kokkos(xnear, cutoff, contacts);
    else
      ncontact = surface_exterior_kokkos(xnear, cutoff, contacts);
    } else {
      // one of surface_int/ext() will return 0
      // so no need to worry about offset of contact indices
      ncontact = surface_exterior_kokkos(xnear, cutoff, contacts) + surface_interior_kokkos(xnear, cutoff, contacts);
    }

    if (rotateflag && ncontact) {
      for (int i = 0; i < ncontact; i++) {
        xs = xnear[0] - contacts[i].delx;
        ys = xnear[1] - contacts[i].dely;
        zs = xnear[2] - contacts[i].delz;
        forward_transform(xs, ys, zs);
        contacts[i].delx = xorig[0] - xs;
        contacts[i].dely = xorig[1] - ys;
        contacts[i].delz = xorig[2] - zs;
      }
    }

    return ncontact;
  }

 private:
  int groupbit;
  typename AT::t_int_1d d_match;
  typename AT::t_kkfloat_1d_3_lr_randomread d_x;
  typename AT::t_int_1d_randomread d_mask;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int surface_interior_kokkos(double *x, double cutoff, contacts_t& contacts)
  {
    double delx = x[0] - xc;
    double dely = x[1] - yc;
    double delz = x[2] - zc;
    double r = sqrt(delx * delx + dely * dely + delz * delz);
    if (r > radius || r == 0.0) return 0;

    double delta = radius - r;
    if (delta < cutoff) {
      contacts[0].r = delta;
      contacts[0].delx = delx * (1.0 - radius / r);
      contacts[0].dely = dely * (1.0 - radius / r);
      contacts[0].delz = delz * (1.0 - radius / r);
      contacts[0].radius = -radius;
      contacts[0].iwall = 0;
      contacts[0].varflag = 1;
      return 1;
    }
    return 0;
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int surface_exterior_kokkos(double *x, double cutoff, contacts_t& contacts)
  {
    double delx = x[0] - xc;
    double dely = x[1] - yc;
    double delz = x[2] - zc;
    double r = sqrt(delx * delx + dely * dely + delz * delz);
    if (r < radius) return 0;

    double delta = r - radius;
    if (delta < cutoff) {
      contacts[0].r = delta;
      contacts[0].delx = delx * (1.0 - radius / r);
      contacts[0].dely = dely * (1.0 - radius / r);
      contacts[0].delz = delz * (1.0 - radius / r);
      contacts[0].radius = radius;
      contacts[0].iwall = 0;
      contacts[0].varflag = 1;
      return 1;
    }
    return 0;
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void add_contact(int n, double *x, double xp, double yp, double zp, contacts_t& contacts)
  {
    double delx = x[0] - xp;
    double dely = x[1] - yp;
    double delz = x[2] - zp;
    contacts[n].r = sqrt(delx * delx + dely * dely + delz * delz);
    contacts[n].radius = 0;
    contacts[n].delx = delx;
    contacts[n].dely = dely;
    contacts[n].delz = delz;
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int k_inside(double x, double y, double z) const
  {
    const double delx = x - xc;
    const double dely = y - yc;
    const double delz = z - zc;
    const double r = sqrt(delx * delx + dely * dely + delz * delz);

    if (r <= radius) return 1;
    return 0;
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void forward_transform(double &x, double &y, double &z) const
  {
    if (rotateflag) rotate(x, y, z, theta);
    if (moveflag) {
      x += dx;
      y += dy;
      z += dz;
    }
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

