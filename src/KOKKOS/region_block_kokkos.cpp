// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "region_block_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "math_extra.h"
#include "memory_kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
RegBlockKokkos<DeviceType>::RegBlockKokkos(LAMMPS *lmp, int narg, char **arg)
  : RegBlock(lmp, narg, arg)
{
  atomKK = (AtomKokkos*) atom;
  memoryKK->create_kokkos(d_contact,6,"region_block:d_contact");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
RegBlockKokkos<DeviceType>::~RegBlockKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(d_contact);
}

/* ----------------------------------------------------------------------
   generate list of contact points for interior or exterior regions
   if region has variable shape, invoke shape_update() once per timestep
   if region is dynamic:
     before: inverse transform x,y,z (unmove, then unrotate)
     after: forward transform contact point xs,yx,zs (rotate, then move),
            then reset contact delx,dely,delz based on new contact point
            no need to do this if no rotation since delxyz doesn't change
   caller is responsible for wrapping this call with
     modify->clearstep_compute() and modify->addstep_compute() if needed
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_FUNCTION
int RegBlockKokkos<DeviceType>::surface_kokkos(double x, double y, double z, double cutoff)
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

/* ----------------------------------------------------------------------
   contact if 0 <= x < cutoff from one or more inner surfaces of block
   can be one contact for each of 6 faces
   no contact if outside (possible if called from union/intersect)
   delxyz = vector from nearest point on block to x
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int RegBlockKokkos<DeviceType>::surface_interior_kokkos(double *x, double cutoff)
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

/* ----------------------------------------------------------------------
   one contact if 0 <= x < cutoff from outer surface of block
   no contact if inside (possible if called from union/intersect)
   delxyz = vector from nearest point on block to x
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int RegBlockKokkos<DeviceType>::surface_exterior_kokkos(double *x, double cutoff)
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

/* ----------------------------------------------------------------------
   add a single contact at Nth location in contact array
   x = particle position
   xp,yp,zp = region surface point
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void RegBlockKokkos<DeviceType>::add_contact(int n, double *x, double xp, double yp, double zp)
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

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int RegBlockKokkos<DeviceType>::k_inside(double x, double y, double z) const
{
  if (x >= xlo && x <= xhi && y >= ylo && y <= yhi && z >= zlo && z <= zhi)
    return 1;
  return 0;
}

template<class DeviceType>
void RegBlockKokkos<DeviceType>::match_all_kokkos(int groupbit_in, DAT::tdual_int_1d k_match_in)
{
  groupbit = groupbit_in;
  d_match = k_match_in.template view<DeviceType>();

  auto execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  atomKK->sync(execution_space, X_MASK | MASK_MASK);

  x = atomKK->k_x.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagRegBlockMatchAll>(0,nlocal),*this);
  copymode = 0;

  k_match_in.template modify<DeviceType>();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void RegBlockKokkos<DeviceType>::operator()(TagRegBlockMatchAll, const int &i) const {
  if (mask[i] & groupbit) {
    double x_tmp = x(i,0);
    double y_tmp = x(i,1);
    double z_tmp = x(i,2);
    d_match[i] = match_kokkos(x_tmp,y_tmp,z_tmp);
  }
}

/* ----------------------------------------------------------------------
   determine if point x,y,z is a match to region volume
   XOR computes 0 if 2 args are the same, 1 if different
   note that k_inside() returns 1 for points on surface of region
   thus point on surface of exterior region will not match
   if region has variable shape, invoke shape_update() once per timestep
   if region is dynamic, apply inverse transform to x,y,z
     unmove first, then unrotate, so don't have to change rotation point
   caller is responsible for wrapping this call with
     modify->clearstep_compute() and modify->addstep_compute() if needed
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int RegBlockKokkos<DeviceType>::match_kokkos(double x, double y, double z) const
{
  if (dynamic) inverse_transform(x,y,z);
  if (openflag) return 1;
  return !(k_inside(x,y,z) ^ interior);
}

/* ----------------------------------------------------------------------
   transform a point x,y,z in region space to moved space
   rotate first (around original P), then displace
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void RegBlockKokkos<DeviceType>::forward_transform(double &x, double &y, double &z) const
{
  if (rotateflag) rotate(x, y, z, theta);
  if (moveflag) {
    x += dx;
    y += dy;
    z += dz;
  }
}

/* ----------------------------------------------------------------------
   transform a point x,y,z in moved space back to region space
   undisplace first, then unrotate (around original P)
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void RegBlockKokkos<DeviceType>::inverse_transform(double &x, double &y, double &z) const
{
  if (moveflag) {
    x -= dx;
    y -= dy;
    z -= dz;
  }
  if (rotateflag) rotate(x,y,z,-theta);
}

/* ----------------------------------------------------------------------
   rotate x,y,z by angle via right-hand rule around point and runit normal
   sign of angle determines whether rotating forward/backward in time
   return updated x,y,z
   R = vector axis of rotation
   P = point = point to rotate around
   R0 = runit = unit vector for R
   X0 = x,y,z = initial coord of atom
   D = X0 - P = vector from P to X0
   C = (D dot R0) R0 = projection of D onto R, i.e. Dparallel
   A = D - C = vector from R line to X0, i.e. Dperp
   B = R0 cross A = vector perp to A in plane of rotation, same len as A
   A,B define plane of circular rotation around R line
   new x,y,z = P + C + A cos(angle) + B sin(angle)
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void RegBlockKokkos<DeviceType>::rotate(double &x, double &y, double &z, double angle) const
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

/* ----------------------------------------------------------------------
   find nearest point to C on line segment A,B and return it as D
   project (C-A) onto (B-A)
   t = length of that projection, normalized by length of (B-A)
   t <= 0, C is closest to A
   t >= 1, C is closest to B
   else closest point is between A and B
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void RegBlockKokkos<DeviceType>::point_on_line_segment(double *a, double *b, double *c, double *d)
{
  double ba[3], ca[3];

  MathExtra::sub3(b, a, ba);
  MathExtra::sub3(c, a, ca);
  double t = MathExtra::dot3(ca, ba) / MathExtra::dot3(ba, ba);
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


/*------------------------------------------------------------------------
  return distance to closest point on surface I of block region
  store closest point in xc,yc,zc
--------------------------------------------------------------------------*/

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double RegBlockKokkos<DeviceType>::find_closest_point(int i, double *x, double &xc, double &yc, double &zc) 
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



namespace LAMMPS_NS {
template class RegBlockKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class RegBlockKokkos<LMPHostType>;
#endif
}

