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

/* ----------------------------------------------------------------------
   Contributing author: Mitch Murphy (alphataubio at gmail)
------------------------------------------------------------------------- */

#include "region_sphere_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "memory_kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
RegSphereKokkos<DeviceType>::RegSphereKokkos(LAMMPS *lmp, int narg, char **arg)
  : RegSphere(lmp, narg, arg)
{
  atomKK = (AtomKokkos*) atom;
  memoryKK->create_kokkos(d_contact,1,"region_sphere:d_contact");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
RegSphereKokkos<DeviceType>::~RegSphereKokkos()
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
KOKKOS_INLINE_FUNCTION
int RegSphereKokkos<DeviceType>::surface(double x, double y, double z, double cutoff)
{
  int ncontact;
  double xs, ys, zs;
  double xnear[3], xorig[3];

  utils::logmesg(lmp, " *** RegSphereKokkos<DeviceType>::surface\n");

  if (dynamic) {
    xorig[0] = x; xorig[1] = y; xorig[2] = z;
    inverse_transform(x, y, z);
  }

  xnear[0] = x; xnear[1] = y; xnear[2] = z;

  if (!openflag) {
    if (interior)
      ncontact = surface_interior(xnear, cutoff);
    else
      ncontact = surface_exterior(xnear, cutoff);
  } else {
    // one of surface_int/ext() will return 0
    // so no need to worry about offset of contact indices
    ncontact = surface_exterior(xnear, cutoff) + surface_interior(xnear, cutoff);
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
   one contact if 0 <= x < cutoff from inner surface of sphere
   no contact if outside (possible if called from union/intersect)
   delxyz = vector from nearest point on sphere to x
   special case: no contact if x is at center of sphere
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int RegSphereKokkos<DeviceType>::surface_interior(double *x, double cutoff)
{
  double delx = x[0] - xc;
  double dely = x[1] - yc;
  double delz = x[2] - zc;
  double r = sqrt(delx * delx + dely * dely + delz * delz);
  if (r > radius || r == 0.0) return 0;

  double delta = radius - r;
  if (delta < cutoff) {
    d_contact[0].r = delta;
    d_contact[0].delx = delx * (1.0 - radius / r);
    d_contact[0].dely = dely * (1.0 - radius / r);
    d_contact[0].delz = delz * (1.0 - radius / r);
    d_contact[0].radius = -radius;
    d_contact[0].iwall = 0;
    d_contact[0].varflag = 1;
    return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   one contact if 0 <= x < cutoff from outer surface of sphere
   no contact if inside (possible if called from union/intersect)
   delxyz = vector from nearest point on sphere to x
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int RegSphereKokkos<DeviceType>::surface_exterior(double *x, double cutoff)
{
  double delx = x[0] - xc;
  double dely = x[1] - yc;
  double delz = x[2] - zc;
  double r = sqrt(delx * delx + dely * dely + delz * delz);
  if (r < radius) return 0;

  double delta = r - radius;
  if (delta < cutoff) {
    d_contact[0].r = delta;
    d_contact[0].delx = delx * (1.0 - radius / r);
    d_contact[0].dely = dely * (1.0 - radius / r);
    d_contact[0].delz = delz * (1.0 - radius / r);
    d_contact[0].radius = radius;
    d_contact[0].iwall = 0;
    d_contact[0].varflag = 1;
    return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int RegSphereKokkos<DeviceType>::k_inside(double x, double y, double z) const
{
  const double delx = x - xc;
  const double dely = y - yc;
  const double delz = z - zc;
  const double r = sqrt(delx * delx + dely * dely + delz * delz);

  if (r <= radius) return 1;
  return 0;
}

template<class DeviceType>
void RegSphereKokkos<DeviceType>::match_all_kokkos(int groupbit_in, DAT::tdual_int_1d k_match_in)
{

  auto execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  atomKK->sync(execution_space, X_MASK | MASK_MASK);

  auto d_x = atomKK->k_x.template view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();
  auto d_match = k_match_in.template view<DeviceType>();
  auto l_groupbit = groupbit_in;

  copymode = 1;

  // capture lambda reference to KOKKOS_INLINE_FUNCTION match()
  // use KOKKOS_CLASS_LAMBDA instead of KOKKOS_LAMBDA
  // https://github.com/kokkos/kokkos/issues/695

  Kokkos::parallel_for(atom->nlocal, KOKKOS_CLASS_LAMBDA( const int &i ) {
    if (d_mask[i] & l_groupbit) {
      double x_tmp = d_x(i,0);
      double y_tmp = d_x(i,1);
      double z_tmp = d_x(i,2);
      d_match[i] = match(x_tmp,y_tmp,z_tmp);
    }});

  copymode = 0;

  k_match_in.template modify<DeviceType>();
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
int RegSphereKokkos<DeviceType>::match(double x, double y, double z) const
{
  if (dynamic) inverse_transform(x,y,z);
  if (openflag) return 1;
  return !(k_inside(x,y,z) ^ interior);
}

/* ----------------------------------------------------------------------
   transform a point x,y,z in moved space back to region space
   undisplace first, then unrotate (around original P)
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void RegSphereKokkos<DeviceType>::inverse_transform(double &x, double &y, double &z) const
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
void RegSphereKokkos<DeviceType>::rotate(double &x, double &y, double &z, double angle) const
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

namespace LAMMPS_NS {
template class RegSphereKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class RegSphereKokkos<LMPHostType>;
#endif
}

