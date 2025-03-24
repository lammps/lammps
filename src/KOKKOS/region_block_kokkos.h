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
RegionStyle(block/kk,RegBlockKokkos<LMPDeviceType>);
RegionStyle(block/kk/device,RegBlockKokkos<LMPDeviceType>);
RegionStyle(block/kk/host,RegBlockKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_REGION_BLOCK_KOKKOS_H
#define LMP_REGION_BLOCK_KOKKOS_H

#include "region_block.h"

#include "kokkos_base.h"
#include "kokkos_type.h"
#include "math_special_kokkos.h"

namespace LAMMPS_NS {

using namespace MathSpecialKokkos;

struct TagRegBlockMatchAll{};

template<class DeviceType>
class RegBlockKokkos : public RegBlock, public KokkosBase  {
  friend class FixPour;

 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  RegBlockKokkos(class LAMMPS *, int, char **);
  ~RegBlockKokkos() override;

  void match_all_kokkos(int, DAT::tdual_int_1d) override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagRegBlockMatchAll, const int&) const;

  KOKKOS_INLINE_FUNCTION
  int match_kokkos(double x, double y, double z) const
  {
    if (dynamic) inverse_transform(x,y,z);
    if (openflag) return 1;
    return !(k_inside(x,y,z) ^ interior);
  }

  KOKKOS_INLINE_FUNCTION
  int surface_kokkos(double x, double y, double z, double cutoff)
  {
    int ncontact;
    double xs, ys, zs;
    double xnear[3], xorig[3];

    if (dynamic) {
      xorig[0] = x; xorig[1] = y; xorig[2] = z;
      inverse_transform(x, y, z);
    }

    xnear[0] = x; xnear[1] = y; xnear[2] = z;

    if (!openflag) {
      if (interior)
        ncontact = surface_interior_kokkos(xnear, cutoff);
      else
        ncontact = surface_exterior_kokkos(xnear, cutoff);
    } else {
      // one of surface_int/ext() will return 0
      // so no need to worry about offset of contact indices
      ncontact = surface_exterior_kokkos(xnear, cutoff) + surface_interior_kokkos(xnear, cutoff);
    }

    if (rotateflag && ncontact) {
      for (int i = 0; i < ncontact; i++) {
        xs = xnear[0] - d_contact[i].delx;
        ys = xnear[1] - d_contact[i].dely;
        zs = xnear[2] - d_contact[i].delz;
        forward_transform(xs, ys, zs);
        d_contact[i].delx = xorig[0] - xs;
        d_contact[i].dely = xorig[1] - ys;
        d_contact[i].delz = xorig[2] - zs;
      }
    }

    return ncontact;
  }

  Kokkos::View<Contact*, DeviceType> d_contact;

 private:
  int groupbit;
  typename AT::t_int_1d d_match;
  typename AT::t_x_array_randomread d_x;
  typename AT::t_int_1d_randomread d_mask;

  KOKKOS_INLINE_FUNCTION
  int surface_interior_kokkos(double *x, double cutoff)
  {
    double delta;

    // x is exterior to block

    if (x[0] < xlo || x[0] > xhi || x[1] < ylo || x[1] > yhi || x[2] < zlo || x[2] > zhi) return 0;

    // x is interior to block or on its surface

    int n = 0;

    delta = x[0] - xlo;
    if (delta < cutoff && !open_faces[0]) {
      d_contact[n].r = delta;
      d_contact[n].delx = delta;
      d_contact[n].dely = d_contact[n].delz = 0.0;
      d_contact[n].radius = 0;
      d_contact[n].iwall = 0;
      n++;
    }
    delta = xhi - x[0];
    if (delta < cutoff && !open_faces[1]) {
      d_contact[n].r = delta;
      d_contact[n].delx = -delta;
      d_contact[n].dely = d_contact[n].delz = 0.0;
      d_contact[n].radius = 0;
      d_contact[n].iwall = 1;
      n++;
    }

    delta = x[1] - ylo;
    if (delta < cutoff && !open_faces[2]) {
      d_contact[n].r = delta;
      d_contact[n].dely = delta;
      d_contact[n].delx = d_contact[n].delz = 0.0;
      d_contact[n].radius = 0;
      d_contact[n].iwall = 2;
      n++;
    }
    delta = yhi - x[1];
    if (delta < cutoff && !open_faces[3]) {
      d_contact[n].r = delta;
      d_contact[n].dely = -delta;
      d_contact[n].delx = d_contact[n].delz = 0.0;
      d_contact[n].radius = 0;
      d_contact[n].iwall = 3;
      n++;
    }

    delta = x[2] - zlo;
    if (delta < cutoff && !open_faces[4]) {
      d_contact[n].r = delta;
      d_contact[n].delz = delta;
      d_contact[n].delx = d_contact[n].dely = 0.0;
      d_contact[n].radius = 0;
      d_contact[n].iwall = 4;
      n++;
    }
    delta = zhi - x[2];
    if (delta < cutoff && !open_faces[5]) {
      d_contact[n].r = delta;
      d_contact[n].delz = -delta;
      d_contact[n].delx = d_contact[n].dely = 0.0;
      d_contact[n].radius = 0;
      d_contact[n].iwall = 5;
      n++;
    }

    return n;
  }

  KOKKOS_INLINE_FUNCTION
  int surface_exterior_kokkos(double *x, double cutoff)
  {
    double xp, yp, zp;
    double xc, yc, zc, dist, mindist;

    // x is far enough from block that there is no contact
    // x is interior to block

    if (x[0] <= xlo - cutoff || x[0] >= xhi + cutoff || x[1] <= ylo - cutoff ||
        x[1] >= yhi + cutoff || x[2] <= zlo - cutoff || x[2] >= zhi + cutoff)
      return 0;
    if (x[0] > xlo && x[0] < xhi && x[1] > ylo && x[1] < yhi && x[2] > zlo && x[2] < zhi) return 0;

    // x is exterior to block or on its surface
    // xp,yp,zp = point on surface of block that x is closest to
    //            could be edge or corner pt of block
    // do not add contact point if r >= cutoff

    if (!openflag) {
      if (x[0] < xlo)
        xp = xlo;
      else if (x[0] > xhi)
        xp = xhi;
      else
        xp = x[0];
      if (x[1] < ylo)
        yp = ylo;
      else if (x[1] > yhi)
        yp = yhi;
      else
        yp = x[1];
      if (x[2] < zlo)
        zp = zlo;
      else if (x[2] > zhi)
        zp = zhi;
      else
        zp = x[2];
    } else {
      mindist = MAXDOUBLEINT;
      for (int i = 0; i < 6; i++) {
        if (open_faces[i]) continue;
        dist = find_closest_point(i, x, xc, yc, zc);
        if (dist < mindist) {
          xp = xc;
          yp = yc;
          zp = zc;
          mindist = dist;
        }
      }
    }

    add_contact(0, x, xp, yp, zp);
    d_contact[0].iwall = 0;
    if (d_contact[0].r < cutoff) return 1;
    return 0;
  }

  KOKKOS_INLINE_FUNCTION
  void add_contact(int n, double *x, double xp, double yp, double zp)
  {
    double delx = x[0] - xp;
    double dely = x[1] - yp;
    double delz = x[2] - zp;
    d_contact[n].r = sqrt(delx * delx + dely * dely + delz * delz);
    d_contact[n].radius = 0;
    d_contact[n].delx = delx;
    d_contact[n].dely = dely;
    d_contact[n].delz = delz;
  }

  KOKKOS_INLINE_FUNCTION
  int k_inside(double x, double y, double z) const
  {
    if (x >= xlo && x <= xhi && y >= ylo && y <= yhi && z >= zlo && z <= zhi)
      return 1;
    return 0;
  }

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

  KOKKOS_INLINE_FUNCTION
  void point_on_line_segment(double *a, double *b, double *c, double *d)
  {
    double ba[3], ca[3];

    sub3(b, a, ba);
    sub3(c, a, ca);
    double t = dot3(ca, ba) / dot3(ba, ba);
    if (t <= 0.0) {
      d[0] = a[0];
      d[1] = a[1];
      d[2] = a[2];
    } else if (t >= 1.0) {
      d[0] = b[0];
      d[1] = b[1];
      d[2] = b[2];
    } else {
      d[0] = a[0] + t * ba[0];
      d[1] = a[1] + t * ba[1];
      d[2] = a[2] + t * ba[2];
    }
  }

  KOKKOS_INLINE_FUNCTION
  double inside_face(double *xproj, int iface)
  {
    if (iface < 2) {
      if (xproj[1] > 0 && (xproj[1] < yhi - ylo) && xproj[2] > 0 && (xproj[2] < zhi - zlo)) return 1;
    } else if (iface < 4) {
      if (xproj[0] > 0 && (xproj[0] < (xhi - xlo)) && xproj[2] > 0 && (xproj[2] < (zhi - zlo)))
        return 1;
    } else {
      if (xproj[0] > 0 && xproj[0] < (xhi - xlo) && xproj[1] > 0 && xproj[1] < (yhi - ylo)) return 1;
    }

    return 0;
  }

  KOKKOS_INLINE_FUNCTION
  double find_closest_point(int i, double *x, double &xc, double &yc, double &zc)
  {
    double dot, d2, d2min;
    double xr[3], xproj[3], p[3];

    xr[0] = x[0] - corners[i][0][0];
    xr[1] = x[1] - corners[i][0][1];
    xr[2] = x[2] - corners[i][0][2];
    dot = face[i][0] * xr[0] + face[i][1] * xr[1] + face[i][2] * xr[2];
    xproj[0] = xr[0] - dot * face[i][0];
    xproj[1] = xr[1] - dot * face[i][1];
    xproj[2] = xr[2] - dot * face[i][2];

    d2min = MAXDOUBLEINT;

    // check if point projects inside of face

    if (inside_face(xproj, i)) {
      d2 = d2min = dot * dot;
      xc = xproj[0] + corners[i][0][0];
      yc = xproj[1] + corners[i][0][1];
      zc = xproj[2] + corners[i][0][2];

      // check each edge

    } else {
      point_on_line_segment(corners[i][0], corners[i][1], x, p);
      d2 = (p[0] - x[0]) * (p[0] - x[0]) + (p[1] - x[1]) * (p[1] - x[1]) +
          (p[2] - x[2]) * (p[2] - x[2]);
      if (d2 < d2min) {
        d2min = d2;
        xc = p[0];
        yc = p[1];
        zc = p[2];
      }

      point_on_line_segment(corners[i][1], corners[i][2], x, p);
      d2 = (p[0] - x[0]) * (p[0] - x[0]) + (p[1] - x[1]) * (p[1] - x[1]) +
          (p[2] - x[2]) * (p[2] - x[2]);
      if (d2 < d2min) {
        d2min = d2;
        xc = p[0];
        yc = p[1];
        zc = p[2];
      }

      point_on_line_segment(corners[i][2], corners[i][3], x, p);
      d2 = (p[0] - x[0]) * (p[0] - x[0]) + (p[1] - x[1]) * (p[1] - x[1]) +
          (p[2] - x[2]) * (p[2] - x[2]);
      if (d2 < d2min) {
        d2min = d2;
        xc = p[0];
        yc = p[1];
        zc = p[2];
      }

      point_on_line_segment(corners[i][3], corners[i][0], x, p);
      d2 = (p[0] - x[0]) * (p[0] - x[0]) + (p[1] - x[1]) * (p[1] - x[1]) +
          (p[2] - x[2]) * (p[2] - x[2]);
      if (d2 < d2min) {
        d2min = d2;
        xc = p[0];
        yc = p[1];
        zc = p[2];
      }
    }

    return d2min;
  }

};

}

#endif
#endif

