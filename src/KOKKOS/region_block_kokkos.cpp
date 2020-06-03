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

template<ExecutionSpace Space>
RegBlockKokkos<Space>::RegBlockKokkos(LAMMPS *lmp, int narg, char **arg) : RegBlock(lmp, narg, arg)
{
  atomKK = (AtomKokkos*) atom;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
RegBlockKokkos<Space>::~RegBlockKokkos()
{

}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
int RegBlockKokkos<Space>::k_inside(KK_FLOAT x, KK_FLOAT y, KK_FLOAT z) const
{
  if (x >= xlo && x <= xhi && y >= ylo && y <= yhi && z >= zlo && z <= zhi)
    return 1;
  return 0;
}

template<ExecutionSpace Space>
void RegBlockKokkos<Space>::match_all_kokkos(int groupbit_in, DAT::tdual_int_1d k_match_in)
{
  groupbit = groupbit_in;
  d_match = DualViewHelper<Space>::view(k_match_in);

  atomKK->sync(Device, X_MASK | MASK_MASK);

  x = DualViewHelper<Space>::view(atomKK->k_x);
  mask = DualViewHelper<Space>::view(atomKK->k_mask);
  int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagRegBlockMatchAll>(0,nlocal),*this);
  copymode = 0;

  DualViewHelper<Space>::modify(k_match_in);
}

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void RegBlockKokkos<Space>::operator()(TagRegBlockMatchAll, const int &i) const {
  if (mask[i] & groupbit) {
    KK_FLOAT x_tmp = x(i,0);
    KK_FLOAT y_tmp = x(i,1);
    KK_FLOAT z_tmp = x(i,2);
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

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
int RegBlockKokkos<Space>::match(KK_FLOAT x, KK_FLOAT y, KK_FLOAT z) const
{
  if (dynamic) inverse_transform(x,y,z);
  return !(k_inside(x,y,z) ^ interior);
}

/* ----------------------------------------------------------------------
   transform a point x,y,z in moved space back to region space
   undisplace first, then unrotate (around original P)
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void RegBlockKokkos<Space>::inverse_transform(KK_FLOAT &x, KK_FLOAT &y, KK_FLOAT &z) const
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

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void RegBlockKokkos<Space>::rotate(KK_FLOAT &x, KK_FLOAT &y, KK_FLOAT &z, KK_FLOAT angle) const
{
  KK_FLOAT a[3],b[3],c[3],d[3],disp[3];

  KK_FLOAT sine = sin(angle);
  KK_FLOAT cosine = cos(angle);
  d[0] = x - point[0];
  d[1] = y - point[1];
  d[2] = z - point[2];
  KK_FLOAT x0dotr = d[0]*runit[0] + d[1]*runit[1] + d[2]*runit[2];
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
template class RegBlockKokkos<Device>;
template class RegBlockKokkos<Host>;
}

