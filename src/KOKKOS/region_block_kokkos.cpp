/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "region_block_kokkos.h"
#include "atom_kokkos.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

template<class DeviceType>
RegBlockKokkos<DeviceType>::RegBlockKokkos(LAMMPS *lmp, int narg, char **arg) : RegBlock(lmp, narg, arg)
{
  atomKK = (AtomKokkos*) atom;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
RegBlockKokkos<DeviceType>::~RegBlockKokkos()
{

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

  atomKK->sync(Device, X_MASK | MASK_MASK);

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
    d_match[i] = match(x_tmp,y_tmp,z_tmp);
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
int RegBlockKokkos<DeviceType>::match(double x, double y, double z) const
{
  if (dynamic) inverse_transform(x,y,z);
  return !(k_inside(x,y,z) ^ interior);
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

namespace LAMMPS_NS {
template class RegBlockKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class RegBlockKokkos<LMPHostType>;
#endif
}

