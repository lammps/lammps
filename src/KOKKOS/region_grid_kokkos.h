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
RegionStyle(grid/kk,RegGridKokkos<LMPDeviceType>);
RegionStyle(grid/kk/device,RegGridKokkos<LMPDeviceType>);
RegionStyle(grid/kk/host,RegGridKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_REGION_GRID_KOKKOS_H
#define LMP_REGION_GRID_KOKKOS_H

#include "region_grid.h"

#include "kokkos_base.h"
#include "kokkos_type.h"
#include "math_special_kokkos.h"

namespace LAMMPS_NS {

using namespace MathSpecialKokkos;

struct TagRegGridMatchAll{};

template<class DeviceType>
class RegGridKokkos : public RegGrid, public KokkosBase {

 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  RegGridKokkos(class LAMMPS *, int, char **);
  ~RegGridKokkos() override;

  void match_all_kokkos(int, DAT::tdual_int_1d) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagRegGridMatchAll, const int&) const;

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
  static constexpr int KOFFSET = 16384;

  int groupbit;
  typename AT::t_int_1d d_match;
  typename AT::t_kkfloat_1d_3_lr_randomread d_x;
  typename AT::t_int_1d_randomread d_mask;

  Kokkos::View<double***, DeviceType> d_gridvals;

  double k_boxlo0, k_boxlo1, k_boxlo2;
  double k_dxinv, k_dyinv, k_dzinv;
  double k_dx, k_dy, k_dz;

  void sync_grid_to_device();

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int k_inside(double x, double y, double z) const
  {
    int ix = static_cast<int>((x - k_boxlo0) * k_dxinv + KOFFSET) - KOFFSET;
    int iy = static_cast<int>((y - k_boxlo1) * k_dyinv + KOFFSET) - KOFFSET;
    int iz = static_cast<int>((z - k_boxlo2) * k_dzinv + KOFFSET) - KOFFSET;

    if (ix < nxlo_out || ix > nxhi_out ||
        iy < nylo_out || iy > nyhi_out ||
        iz < nzlo_out || iz > nzhi_out) return 0;

    double value = d_gridvals(iz - nzlo_out, iy - nylo_out, ix - nxlo_out);
    return k_evaluate(value);
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  bool k_cell_inside(int ix, int iy, int iz) const
  {
    if (ix < nxlo_out || ix > nxhi_out ||
        iy < nylo_out || iy > nyhi_out ||
        iz < nzlo_out || iz > nzhi_out) return false;

    double value = d_gridvals(iz - nzlo_out, iy - nylo_out, ix - nxlo_out);
    return k_evaluate(value) == 1;
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int k_evaluate(double value) const
  {
    switch (compare_op) {
      case OP_GT: return (value > threshold) ? 1 : 0;
      case OP_GE: return (value >= threshold) ? 1 : 0;
      case OP_LT: return (value < threshold) ? 1 : 0;
      case OP_LE: return (value <= threshold) ? 1 : 0;
      case OP_EQ: return (value == threshold) ? 1 : 0;
      case OP_NE: return (value != threshold) ? 1 : 0;
    }
    return 0;
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int surface_interior_kokkos(double *x, double cutoff)
  {
    int ix = static_cast<int>((x[0] - k_boxlo0) / k_dx + KOFFSET) - KOFFSET;
    int iy = static_cast<int>((x[1] - k_boxlo1) / k_dy + KOFFSET) - KOFFSET;
    int iz = static_cast<int>((x[2] - k_boxlo2) / k_dz + KOFFSET) - KOFFSET;

    if (!k_cell_inside(ix, iy, iz)) return 0;

    int n = 0;
    double delta;

    delta = x[0] - (k_boxlo0 + ix * k_dx);
    if (delta < cutoff && !k_cell_inside(ix - 1, iy, iz)) {
      d_contact[n].r = delta;
      d_contact[n].delx = delta;
      d_contact[n].dely = d_contact[n].delz = 0.0;
      d_contact[n].radius = 0;
      d_contact[n].iwall = 0;
      n++;
    }

    delta = (k_boxlo0 + (ix + 1) * k_dx) - x[0];
    if (delta < cutoff && !k_cell_inside(ix + 1, iy, iz)) {
      d_contact[n].r = delta;
      d_contact[n].delx = -delta;
      d_contact[n].dely = d_contact[n].delz = 0.0;
      d_contact[n].radius = 0;
      d_contact[n].iwall = 1;
      n++;
    }

    delta = x[1] - (k_boxlo1 + iy * k_dy);
    if (delta < cutoff && !k_cell_inside(ix, iy - 1, iz)) {
      d_contact[n].r = delta;
      d_contact[n].dely = delta;
      d_contact[n].delx = d_contact[n].delz = 0.0;
      d_contact[n].radius = 0;
      d_contact[n].iwall = 2;
      n++;
    }

    delta = (k_boxlo1 + (iy + 1) * k_dy) - x[1];
    if (delta < cutoff && !k_cell_inside(ix, iy + 1, iz)) {
      d_contact[n].r = delta;
      d_contact[n].dely = -delta;
      d_contact[n].delx = d_contact[n].delz = 0.0;
      d_contact[n].radius = 0;
      d_contact[n].iwall = 3;
      n++;
    }

    delta = x[2] - (k_boxlo2 + iz * k_dz);
    if (delta < cutoff && !k_cell_inside(ix, iy, iz - 1)) {
      d_contact[n].r = delta;
      d_contact[n].delz = delta;
      d_contact[n].delx = d_contact[n].dely = 0.0;
      d_contact[n].radius = 0;
      d_contact[n].iwall = 4;
      n++;
    }

    delta = (k_boxlo2 + (iz + 1) * k_dz) - x[2];
    if (delta < cutoff && !k_cell_inside(ix, iy, iz + 1)) {
      d_contact[n].r = delta;
      d_contact[n].delz = -delta;
      d_contact[n].delx = d_contact[n].dely = 0.0;
      d_contact[n].radius = 0;
      d_contact[n].iwall = 5;
      n++;
    }

    return n;
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int surface_exterior_kokkos(double *x, double cutoff)
  {
    int ix = static_cast<int>((x[0] - k_boxlo0) / k_dx + KOFFSET) - KOFFSET;
    int iy = static_cast<int>((x[1] - k_boxlo1) / k_dy + KOFFSET) - KOFFSET;
    int iz = static_cast<int>((x[2] - k_boxlo2) / k_dz + KOFFSET) - KOFFSET;

    if (ix < nxlo_out || ix > nxhi_out ||
        iy < nylo_out || iy > nyhi_out ||
        iz < nzlo_out || iz > nzhi_out) return 0;

    if (k_cell_inside(ix, iy, iz)) return 0;

    double mindist = cutoff;
    double best_delx = 0.0, best_dely = 0.0, best_delz = 0.0;
    int found = 0;
    double delta;

    delta = x[0] - (k_boxlo0 + ix * k_dx);
    if (delta < mindist && k_cell_inside(ix - 1, iy, iz)) {
      mindist = delta;
      best_delx = delta;
      best_dely = best_delz = 0.0;
      found = 1;
    }

    delta = (k_boxlo0 + (ix + 1) * k_dx) - x[0];
    if (delta < mindist && k_cell_inside(ix + 1, iy, iz)) {
      mindist = delta;
      best_delx = -delta;
      best_dely = best_delz = 0.0;
      found = 1;
    }

    delta = x[1] - (k_boxlo1 + iy * k_dy);
    if (delta < mindist && k_cell_inside(ix, iy - 1, iz)) {
      mindist = delta;
      best_dely = delta;
      best_delx = best_delz = 0.0;
      found = 1;
    }

    delta = (k_boxlo1 + (iy + 1) * k_dy) - x[1];
    if (delta < mindist && k_cell_inside(ix, iy + 1, iz)) {
      mindist = delta;
      best_dely = -delta;
      best_delx = best_delz = 0.0;
      found = 1;
    }

    delta = x[2] - (k_boxlo2 + iz * k_dz);
    if (delta < mindist && k_cell_inside(ix, iy, iz - 1)) {
      mindist = delta;
      best_delz = delta;
      best_delx = best_dely = 0.0;
      found = 1;
    }

    delta = (k_boxlo2 + (iz + 1) * k_dz) - x[2];
    if (delta < mindist && k_cell_inside(ix, iy, iz + 1)) {
      mindist = delta;
      best_delz = -delta;
      best_delx = best_dely = 0.0;
      found = 1;
    }

    if (!found) return 0;

    d_contact[0].r = mindist;
    d_contact[0].delx = best_delx;
    d_contact[0].dely = best_dely;
    d_contact[0].delz = best_delz;
    d_contact[0].radius = 0;
    d_contact[0].iwall = 0;
    return 1;
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
