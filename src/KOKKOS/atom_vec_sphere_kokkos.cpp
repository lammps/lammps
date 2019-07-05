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

#include "atom_vec_sphere_kokkos.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm_kokkos.h"
#include "domain.h"
#include "modify.h"
#include "force.h"
#include "fix.h"
#include "fix_adapt.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "memory_kokkos.h"

using namespace LAMMPS_NS;

#define DELTA 10

static const double MY_PI  = 3.14159265358979323846; // pi

/* ---------------------------------------------------------------------- */

AtomVecSphereKokkos::AtomVecSphereKokkos(LAMMPS *lmp) : AtomVecKokkos(lmp)
{
  molecular = 0;

  comm_x_only = 1;
  comm_f_only = 0;
  size_forward = 3;
  size_reverse = 6;
  size_border = 8;
  size_velocity = 6;
  size_data_atom = 7;
  size_data_vel = 7;
  xcol_data = 5;

  atom->sphere_flag = 1;
  atom->radius_flag = atom->rmass_flag = atom->omega_flag =
    atom->torque_flag = 1;

  k_count = DAT::tdual_int_1d("atom::k_count",1);
  atomKK = (AtomKokkos *) atom;
  commKK = (CommKokkos *) comm;

  no_border_vel_flag = 0;
}

/* ---------------------------------------------------------------------- */

void AtomVecSphereKokkos::init()
{
  AtomVec::init();

  // set radvary if particle diameters are time-varying due to fix adapt

  radvary = 0;
  comm_x_only = 1;
  size_forward = 3;

  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"adapt") == 0) {
      FixAdapt *fix = (FixAdapt *) modify->fix[i];
      if (fix->diamflag) {
        radvary = 1;
        comm_x_only = 0;
        size_forward = 5;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecSphereKokkos::grow(int n)
{
  int step = MAX(DELTA,nmax*0.01);
  if (n == 0) nmax += step;
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  atomKK->sync(Device,ALL_MASK);
  atomKK->modified(Device,ALL_MASK);

  memoryKK->grow_kokkos(atomKK->k_tag,atomKK->tag,nmax,"atom:tag");
  memoryKK->grow_kokkos(atomKK->k_type,atomKK->type,nmax,"atom:type");
  memoryKK->grow_kokkos(atomKK->k_mask,atomKK->mask,nmax,"atom:mask");
  memoryKK->grow_kokkos(atomKK->k_image,atomKK->image,nmax,"atom:image");

  memoryKK->grow_kokkos(atomKK->k_x,atomKK->x,nmax,3,"atom:x");
  memoryKK->grow_kokkos(atomKK->k_v,atomKK->v,nmax,3,"atom:v");
  memoryKK->grow_kokkos(atomKK->k_f,atomKK->f,nmax,3,"atom:f");
  memoryKK->grow_kokkos(atomKK->k_radius,atomKK->radius,nmax,"atom:radius");
  memoryKK->grow_kokkos(atomKK->k_rmass,atomKK->rmass,nmax,"atom:rmass");
  memoryKK->grow_kokkos(atomKK->k_omega,atomKK->omega,nmax,3,"atom:omega");
  memoryKK->grow_kokkos(atomKK->k_torque,atomKK->torque,nmax,3,"atom:torque");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);

  grow_reset();
  atomKK->sync(Host,ALL_MASK);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecSphereKokkos::grow_reset()
{
  tag = atomKK->tag;
  d_tag = atomKK->k_tag.d_view;
  h_tag = atomKK->k_tag.h_view;

  type = atomKK->type;
  d_type = atomKK->k_type.d_view;
  h_type = atomKK->k_type.h_view;
  mask = atomKK->mask;
  d_mask = atomKK->k_mask.d_view;
  h_mask = atomKK->k_mask.h_view;
  image = atomKK->image;
  d_image = atomKK->k_image.d_view;
  h_image = atomKK->k_image.h_view;

  x = atomKK->x;
  d_x = atomKK->k_x.d_view;
  h_x = atomKK->k_x.h_view;
  v = atomKK->v;
  d_v = atomKK->k_v.d_view;
  h_v = atomKK->k_v.h_view;
  f = atomKK->f;
  d_f = atomKK->k_f.d_view;
  h_f = atomKK->k_f.h_view;
  radius = atomKK->radius;
  d_radius = atomKK->k_radius.d_view;
  h_radius = atomKK->k_radius.h_view;
  rmass = atomKK->rmass;
  d_rmass = atomKK->k_rmass.d_view;
  h_rmass = atomKK->k_rmass.h_view;
  omega = atomKK->omega;
  d_omega = atomKK->k_omega.d_view;
  h_omega = atomKK->k_omega.h_view;
  torque = atomKK->torque;
  d_torque = atomKK->k_torque.d_view;
  h_torque = atomKK->k_torque.h_view;
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecSphereKokkos::copy(int i, int j, int delflag)
{
  atomKK->sync(Host,X_MASK | V_MASK | TAG_MASK | TYPE_MASK |
            MASK_MASK | IMAGE_MASK | RADIUS_MASK |
            RMASS_MASK | OMEGA_MASK);

  h_tag[j] = h_tag[i];
  h_type[j] = h_type[i];
  h_mask[j] = h_mask[i];
  h_image[j] = h_image[i];
  h_x(j,0) = h_x(i,0);
  h_x(j,1) = h_x(i,1);
  h_x(j,2) = h_x(i,2);
  h_v(j,0) = h_v(i,0);
  h_v(j,1) = h_v(i,1);
  h_v(j,2) = h_v(i,2);

  h_radius[j] = h_radius[i];
  h_rmass[j] = h_rmass[i];
  h_omega(j,0) = h_omega(i,0);
  h_omega(j,1) = h_omega(i,1);
  h_omega(j,2) = h_omega(i,2);

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);

  atomKK->modified(Host,X_MASK | V_MASK | TAG_MASK | TYPE_MASK |
                MASK_MASK | IMAGE_MASK | RADIUS_MASK |
                RMASS_MASK | OMEGA_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG,int TRICLINIC>
struct AtomVecSphereKokkos_PackComm {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_x_array_randomread _x;
  typename ArrayTypes<DeviceType>::t_float_1d _radius,_rmass;
  typename ArrayTypes<DeviceType>::t_xfloat_2d_um _buf;
  typename ArrayTypes<DeviceType>::t_int_2d_const _list;
  const int _iswap;
  X_FLOAT _xprd,_yprd,_zprd,_xy,_xz,_yz;
  X_FLOAT _pbc[6];

  AtomVecSphereKokkos_PackComm(
    const typename DAT::tdual_x_array &x,
    const typename DAT::tdual_float_1d &radius,
    const typename DAT::tdual_float_1d &rmass,
    const typename DAT::tdual_xfloat_2d &buf,
    const typename DAT::tdual_int_2d &list,
    const int & iswap,
    const X_FLOAT &xprd, const X_FLOAT &yprd, const X_FLOAT &zprd,
    const X_FLOAT &xy, const X_FLOAT &xz, const X_FLOAT &yz, const int* const pbc):
    _x(x.view<DeviceType>()),
    _radius(radius.view<DeviceType>()),
    _rmass(rmass.view<DeviceType>()),
    _list(list.view<DeviceType>()),_iswap(iswap),
    _xprd(xprd),_yprd(yprd),_zprd(zprd),
    _xy(xy),_xz(xz),_yz(yz) {
    const size_t elements = 5;
    const size_t maxsend = (buf.view<DeviceType>().extent(0)*buf.view<DeviceType>().extent(1))/elements;
    _buf = typename ArrayTypes<DeviceType>::t_xfloat_2d_um(buf.view<DeviceType>().data(),maxsend,elements);
    _pbc[0] = pbc[0]; _pbc[1] = pbc[1]; _pbc[2] = pbc[2];
    _pbc[3] = pbc[3]; _pbc[4] = pbc[4]; _pbc[5] = pbc[5];
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(_iswap,i);
    if (PBC_FLAG == 0) {
      _buf(i,0) = _x(j,0);
      _buf(i,1) = _x(j,1);
      _buf(i,2) = _x(j,2);
    } else {
      if (TRICLINIC == 0) {
	_buf(i,0) = _x(j,0) + _pbc[0]*_xprd;
	_buf(i,1) = _x(j,1) + _pbc[1]*_yprd;
	_buf(i,2) = _x(j,2) + _pbc[2]*_zprd;
      } else {
	_buf(i,0) = _x(j,0) + _pbc[0]*_xprd + _pbc[5]*_xy + _pbc[4]*_xz;
	_buf(i,1) = _x(j,1) + _pbc[1]*_yprd + _pbc[3]*_yz;
	_buf(i,2) = _x(j,2) + _pbc[2]*_zprd;
      }
    }
    _buf(i,3) = _radius(j);
    _buf(i,4) = _rmass(j);
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecSphereKokkos::pack_comm_kokkos(
  const int &n,
  const DAT::tdual_int_2d &list,
  const int & iswap,
  const DAT::tdual_xfloat_2d &buf,
  const int &pbc_flag,
  const int* const pbc)
{
  // Fallback to AtomVecKokkos if radvary == 0
  if (radvary == 0)
    return AtomVecKokkos::pack_comm_kokkos(n,list,iswap,buf,pbc_flag,pbc);
  // Check whether to always run forward communication on the host
  // Choose correct forward PackComm kernel
  if(commKK->forward_comm_on_host) {
    atomKK->sync(Host,X_MASK|RADIUS_MASK|RMASS_MASK);
    if(pbc_flag) {
      if(domain->triclinic) {
        struct AtomVecSphereKokkos_PackComm<LMPHostType,1,1> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          buf,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecSphereKokkos_PackComm<LMPHostType,1,0> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          buf,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      }
    } else {
      if(domain->triclinic) {
        struct AtomVecSphereKokkos_PackComm<LMPHostType,0,1> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          buf,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecSphereKokkos_PackComm<LMPHostType,0,0> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          buf,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      }
    }
  } else {
    atomKK->sync(Device,X_MASK|RADIUS_MASK|RMASS_MASK);
    if(pbc_flag) {
      if(domain->triclinic) {
        struct AtomVecSphereKokkos_PackComm<LMPDeviceType,1,1> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          buf,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecSphereKokkos_PackComm<LMPDeviceType,1,0> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          buf,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      }
    } else {
      if(domain->triclinic) {
        struct AtomVecSphereKokkos_PackComm<LMPDeviceType,0,1> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          buf,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecSphereKokkos_PackComm<LMPDeviceType,0,0> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          buf,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      }
    }
  }
  return n*size_forward;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int RADVARY,int PBC_FLAG,int TRICLINIC,int DEFORM_VREMAP>
struct AtomVecSphereKokkos_PackCommVel {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_x_array_randomread _x;
  typename ArrayTypes<DeviceType>::t_int_1d _mask;
  typename ArrayTypes<DeviceType>::t_float_1d _radius,_rmass;
  typename ArrayTypes<DeviceType>::t_v_array _v, _omega;
  typename ArrayTypes<DeviceType>::t_xfloat_2d_um _buf;
  typename ArrayTypes<DeviceType>::t_int_2d_const _list;
  const int _iswap;
  X_FLOAT _xprd,_yprd,_zprd,_xy,_xz,_yz;
  X_FLOAT _pbc[6];
  X_FLOAT _h_rate[6];
  const int _deform_vremap;

  AtomVecSphereKokkos_PackCommVel(
    const typename DAT::tdual_x_array &x,
    const typename DAT::tdual_int_1d &mask,
    const typename DAT::tdual_float_1d &radius,
    const typename DAT::tdual_float_1d &rmass,
    const typename DAT::tdual_v_array &v,
    const typename DAT::tdual_v_array &omega,
    const typename DAT::tdual_xfloat_2d &buf,
    const typename DAT::tdual_int_2d &list,
    const int &iswap,
    const X_FLOAT &xprd, const X_FLOAT &yprd, const X_FLOAT &zprd,
    const X_FLOAT &xy, const X_FLOAT &xz, const X_FLOAT &yz, const int* const pbc,
    const double * const h_rate,
    const int &deform_vremap):
    _x(x.view<DeviceType>()),
    _mask(mask.view<DeviceType>()),
    _radius(radius.view<DeviceType>()),
    _rmass(rmass.view<DeviceType>()),
    _v(v.view<DeviceType>()),
    _omega(omega.view<DeviceType>()),
    _list(list.view<DeviceType>()),_iswap(iswap),
    _xprd(xprd),_yprd(yprd),_zprd(zprd),
    _xy(xy),_xz(xz),_yz(yz),
    _deform_vremap(deform_vremap)
  {
    const size_t elements = 9 + 2 * RADVARY;
    const int maxsend = (buf.template view<DeviceType>().extent(0)*buf.template view<DeviceType>().extent(1))/elements;
    _buf = typename ArrayTypes<DeviceType>::t_xfloat_2d_um(buf.view<DeviceType>().data(),maxsend,elements);
    _pbc[0] = pbc[0]; _pbc[1] = pbc[1]; _pbc[2] = pbc[2];
    _pbc[3] = pbc[3]; _pbc[4] = pbc[4]; _pbc[5] = pbc[5];
    _h_rate[0] = h_rate[0]; _h_rate[1] = h_rate[1]; _h_rate[2] = h_rate[2];
    _h_rate[3] = h_rate[3]; _h_rate[4] = h_rate[4]; _h_rate[5] = h_rate[5];
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(_iswap,i);
    if (PBC_FLAG == 0) {
      _buf(i,0) = _x(j,0);
      _buf(i,1) = _x(j,1);
      _buf(i,2) = _x(j,2);
    } else {
      if (TRICLINIC == 0) {
	_buf(i,0) = _x(j,0) + _pbc[0]*_xprd;
	_buf(i,1) = _x(j,1) + _pbc[1]*_yprd;
	_buf(i,2) = _x(j,2) + _pbc[2]*_zprd;
      } else {
	_buf(i,0) = _x(j,0) + _pbc[0]*_xprd + _pbc[5]*_xy + _pbc[4]*_xz;
	_buf(i,1) = _x(j,1) + _pbc[1]*_yprd + _pbc[3]*_yz;
	_buf(i,2) = _x(j,2) + _pbc[2]*_zprd;
      }
    }
    if (DEFORM_VREMAP == 0) {
      _buf(i,3) = _v(j,0);
      _buf(i,4) = _v(j,1);
      _buf(i,5) = _v(j,2);
    } else {
      if (_mask(i) & _deform_vremap) {
        _buf(i,3) = _v(j,0) + _pbc[0]*_h_rate[0] + _pbc[5]*_h_rate[5] + _pbc[4]*_h_rate[4];
        _buf(i,4) = _v(j,1) + _pbc[1]*_h_rate[1] + _pbc[3]*_h_rate[3];
        _buf(i,5) = _v(j,2) + _pbc[2]*_h_rate[2];
      } else {
        _buf(i,3) = _v(j,0);
        _buf(i,4) = _v(j,1);
        _buf(i,5) = _v(j,2);
      }
    }
    _buf(i,6) = _omega(j,0);
    _buf(i,7) = _omega(j,1);
    _buf(i,8) = _omega(j,2);
    if (RADVARY) {
      _buf(i,9) = _radius(j);
      _buf(i,10) = _rmass(j);
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecSphereKokkos::pack_comm_vel_kokkos(
  const int &n,
  const DAT::tdual_int_2d &list,
  const int & iswap,
  const DAT::tdual_xfloat_2d &buf,
  const int &pbc_flag,
  const int* const pbc)
{
  if(commKK->forward_comm_on_host) {
    atomKK->sync(Host,X_MASK|RADIUS_MASK|RMASS_MASK|V_MASK|OMEGA_MASK);
    if(pbc_flag) {
      if(deform_vremap) {
        if(domain->triclinic) {
          if (radvary == 0) {
            struct AtomVecSphereKokkos_PackCommVel<LMPHostType,0,1,1,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecSphereKokkos_PackCommVel<LMPHostType,1,1,1,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          }
        } else {
          if (radvary == 0) {
            struct AtomVecSphereKokkos_PackCommVel<LMPHostType,0,1,0,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecSphereKokkos_PackCommVel<LMPHostType,1,1,0,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          }
        }
      } else {
        if(domain->triclinic) {
          if (radvary == 0) {
            struct AtomVecSphereKokkos_PackCommVel<LMPHostType,0,1,1,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecSphereKokkos_PackCommVel<LMPHostType,1,1,1,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          }
        } else {
          if (radvary == 0) {
            struct AtomVecSphereKokkos_PackCommVel<LMPHostType,0,1,0,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecSphereKokkos_PackCommVel<LMPHostType,1,1,0,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          }
        }
      }
    } else {
      if(domain->triclinic) {
        if (radvary == 0) {
          struct AtomVecSphereKokkos_PackCommVel<LMPHostType,0,0,1,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecSphereKokkos_PackCommVel<LMPHostType,1,0,1,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        }
      } else {
        if (radvary == 0) {
          struct AtomVecSphereKokkos_PackCommVel<LMPHostType,0,0,0,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecSphereKokkos_PackCommVel<LMPHostType,1,0,0,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        }
      }
    }
  } else {
    atomKK->sync(Device,X_MASK|RADIUS_MASK|RMASS_MASK|V_MASK|OMEGA_MASK);
    if(pbc_flag) {
      if(deform_vremap) {
        if(domain->triclinic) {
          if (radvary == 0) {
            struct AtomVecSphereKokkos_PackCommVel<LMPDeviceType,0,1,1,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecSphereKokkos_PackCommVel<LMPDeviceType,1,1,1,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          }
        } else {
          if (radvary == 0) {
            struct AtomVecSphereKokkos_PackCommVel<LMPDeviceType,0,1,0,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecSphereKokkos_PackCommVel<LMPDeviceType,1,1,0,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          }
        }
      } else {
        if(domain->triclinic) {
          if (radvary == 0) {
            struct AtomVecSphereKokkos_PackCommVel<LMPDeviceType,0,1,1,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecSphereKokkos_PackCommVel<LMPDeviceType,1,1,1,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          }
        } else {
          if (radvary == 0) {
            struct AtomVecSphereKokkos_PackCommVel<LMPDeviceType,0,1,0,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecSphereKokkos_PackCommVel<LMPDeviceType,1,1,0,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          }
        }
      }
    } else {
      if(domain->triclinic) {
        if (radvary == 0) {
          struct AtomVecSphereKokkos_PackCommVel<LMPDeviceType,0,0,1,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecSphereKokkos_PackCommVel<LMPDeviceType,1,0,1,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        }
      } else {
        if (radvary == 0) {
          struct AtomVecSphereKokkos_PackCommVel<LMPDeviceType,0,0,0,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecSphereKokkos_PackCommVel<LMPDeviceType,1,0,0,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        }
      }
    }
  }
  return n*(size_forward+size_velocity);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG,int TRICLINIC>
struct AtomVecSphereKokkos_PackCommSelf {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_x_array_randomread _x;
  typename ArrayTypes<DeviceType>::t_x_array _xw;
  typename ArrayTypes<DeviceType>::t_float_1d _radius,_rmass;
  int _nfirst;
  typename ArrayTypes<DeviceType>::t_int_2d_const _list;
  const int _iswap;
  X_FLOAT _xprd,_yprd,_zprd,_xy,_xz,_yz;
  X_FLOAT _pbc[6];

  AtomVecSphereKokkos_PackCommSelf(
    const typename DAT::tdual_x_array &x,
    const typename DAT::tdual_float_1d &radius,
    const typename DAT::tdual_float_1d &rmass,
    const int &nfirst,
    const typename DAT::tdual_int_2d &list,
    const int & iswap,
    const X_FLOAT &xprd, const X_FLOAT &yprd, const X_FLOAT &zprd,
    const X_FLOAT &xy, const X_FLOAT &xz, const X_FLOAT &yz, const int* const pbc):
    _x(x.view<DeviceType>()),_xw(x.view<DeviceType>()),
    _radius(radius.view<DeviceType>()),
    _rmass(rmass.view<DeviceType>()),
    _nfirst(nfirst),_list(list.view<DeviceType>()),_iswap(iswap),
    _xprd(xprd),_yprd(yprd),_zprd(zprd),
    _xy(xy),_xz(xz),_yz(yz) {
    _pbc[0] = pbc[0]; _pbc[1] = pbc[1]; _pbc[2] = pbc[2];
    _pbc[3] = pbc[3]; _pbc[4] = pbc[4]; _pbc[5] = pbc[5];
  };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(_iswap,i);
    if (PBC_FLAG == 0) {
      _xw(i+_nfirst,0) = _x(j,0);
      _xw(i+_nfirst,1) = _x(j,1);
      _xw(i+_nfirst,2) = _x(j,2);
    } else {
      if (TRICLINIC == 0) {
	_xw(i+_nfirst,0) = _x(j,0) + _pbc[0]*_xprd;
	_xw(i+_nfirst,1) = _x(j,1) + _pbc[1]*_yprd;
	_xw(i+_nfirst,2) = _x(j,2) + _pbc[2]*_zprd;
      } else {
	_xw(i+_nfirst,0) = _x(j,0) + _pbc[0]*_xprd + _pbc[5]*_xy + _pbc[4]*_xz;
	_xw(i+_nfirst,1) = _x(j,1) + _pbc[1]*_yprd + _pbc[3]*_yz;
	_xw(i+_nfirst,2) = _x(j,2) + _pbc[2]*_zprd;
      }
    }
    _radius(i+_nfirst) = _radius(j);
    _rmass(i+_nfirst) = _rmass(j);
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecSphereKokkos::pack_comm_self(
  const int &n, const DAT::tdual_int_2d &list, const int &iswap,
  const int nfirst, const int &pbc_flag, const int* const pbc) {
  // Fallback to AtomVecKokkos if radvary == 0
  if (radvary == 0)
    return AtomVecKokkos::pack_comm_self(n,list,iswap,nfirst,pbc_flag,pbc);
  if(commKK->forward_comm_on_host) {
    atomKK->sync(Host,X_MASK|RADIUS_MASK|RMASS_MASK);
    atomKK->modified(Host,X_MASK|RADIUS_MASK|RMASS_MASK);
    if(pbc_flag) {
      if(domain->triclinic) {
	struct AtomVecSphereKokkos_PackCommSelf<LMPHostType,1,1> f(
          atomKK->k_x,
	  atomKK->k_radius,atomKK->k_rmass,
	  nfirst,list,iswap,
	  domain->xprd,domain->yprd,domain->zprd,
	  domain->xy,domain->xz,domain->yz,pbc);
	Kokkos::parallel_for(n,f);
      } else {
	struct AtomVecSphereKokkos_PackCommSelf<LMPHostType,1,0> f(
          atomKK->k_x,
	  atomKK->k_radius,atomKK->k_rmass,
	  nfirst,list,iswap,
	  domain->xprd,domain->yprd,domain->zprd,
	  domain->xy,domain->xz,domain->yz,pbc);
	Kokkos::parallel_for(n,f);
      }
    } else {
      if(domain->triclinic) {
	struct AtomVecSphereKokkos_PackCommSelf<LMPHostType,0,1> f(
          atomKK->k_x,
	  atomKK->k_radius,atomKK->k_rmass,
	  nfirst,list,iswap,
	  domain->xprd,domain->yprd,domain->zprd,
	  domain->xy,domain->xz,domain->yz,pbc);
	Kokkos::parallel_for(n,f);
      } else {
	struct AtomVecSphereKokkos_PackCommSelf<LMPHostType,0,0> f(
          atomKK->k_x,
	  atomKK->k_radius,atomKK->k_rmass,
	  nfirst,list,iswap,
	  domain->xprd,domain->yprd,domain->zprd,
	  domain->xy,domain->xz,domain->yz,pbc);
	Kokkos::parallel_for(n,f);
      }
    }
  } else {
    atomKK->sync(Device,X_MASK|RADIUS_MASK|RMASS_MASK);
    atomKK->modified(Device,X_MASK|RADIUS_MASK|RMASS_MASK);
    if(pbc_flag) {
      if(domain->triclinic) {
	struct AtomVecSphereKokkos_PackCommSelf<LMPDeviceType,1,1> f(
          atomKK->k_x,
	  atomKK->k_radius,atomKK->k_rmass,
	  nfirst,list,iswap,
	  domain->xprd,domain->yprd,domain->zprd,
	  domain->xy,domain->xz,domain->yz,pbc);
	Kokkos::parallel_for(n,f);
      } else {
	struct AtomVecSphereKokkos_PackCommSelf<LMPDeviceType,1,0> f(
          atomKK->k_x,
	  atomKK->k_radius,atomKK->k_rmass,
	  nfirst,list,iswap,
	  domain->xprd,domain->yprd,domain->zprd,
	  domain->xy,domain->xz,domain->yz,pbc);
	Kokkos::parallel_for(n,f);
      }
    } else {
      if(domain->triclinic) {
	struct AtomVecSphereKokkos_PackCommSelf<LMPDeviceType,0,1> f(
          atomKK->k_x,
	  atomKK->k_radius,atomKK->k_rmass,
	  nfirst,list,iswap,
	  domain->xprd,domain->yprd,domain->zprd,
	  domain->xy,domain->xz,domain->yz,pbc);
	Kokkos::parallel_for(n,f);
      } else {
	struct AtomVecSphereKokkos_PackCommSelf<LMPDeviceType,0,0> f(
          atomKK->k_x,
	  atomKK->k_radius,atomKK->k_rmass,
	  nfirst,list,iswap,
	  domain->xprd,domain->yprd,domain->zprd,
	  domain->xy,domain->xz,domain->yz,pbc);
	Kokkos::parallel_for(n,f);
      }
    }
  }
  return n*size_forward;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecSphereKokkos_UnpackComm {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_x_array _x;
  typename ArrayTypes<DeviceType>::t_float_1d _radius,_rmass;
  typename ArrayTypes<DeviceType>::t_xfloat_2d_const_um _buf;
  int _first;

  AtomVecSphereKokkos_UnpackComm(
    const typename DAT::tdual_x_array &x,
    const typename DAT::tdual_float_1d &radius,
    const typename DAT::tdual_float_1d &rmass,
    const typename DAT::tdual_xfloat_2d &buf,
    const int& first):
    _x(x.view<DeviceType>()),
    _radius(radius.view<DeviceType>()),
    _rmass(rmass.view<DeviceType>()),
    _first(first)
  {
    const size_t elements = 5;
    const size_t maxsend = (buf.view<DeviceType>().extent(0)*buf.view<DeviceType>().extent(1))/elements;
    _buf = typename ArrayTypes<DeviceType>::t_xfloat_2d_const_um(buf.view<DeviceType>().data(),maxsend,elements);
  };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    _x(i+_first,0) = _buf(i,0);
    _x(i+_first,1) = _buf(i,1);
    _x(i+_first,2) = _buf(i,2);
    _radius(i+_first) = _buf(i,3);
    _rmass(i+_first) = _buf(i,4);
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecSphereKokkos::unpack_comm_kokkos(
  const int &n, const int &first,
  const DAT::tdual_xfloat_2d &buf ) {
  // Fallback to AtomVecKokkos if radvary == 0
  if (radvary == 0) {
    AtomVecKokkos::unpack_comm_kokkos(n,first,buf);
    return;
  }
  if(commKK->forward_comm_on_host) {
    atomKK->modified(Host,X_MASK|RADIUS_MASK|RMASS_MASK);
    struct AtomVecSphereKokkos_UnpackComm<LMPHostType> f(
      atomKK->k_x,
      atomKK->k_radius,atomKK->k_rmass,
      buf,first);
    Kokkos::parallel_for(n,f);
  } else {
    atomKK->modified(Device,X_MASK|RADIUS_MASK|RMASS_MASK);
    struct AtomVecSphereKokkos_UnpackComm<LMPDeviceType> f(
      atomKK->k_x,
      atomKK->k_radius,atomKK->k_rmass,
      buf,first);
    Kokkos::parallel_for(n,f);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int RADVARY>
struct AtomVecSphereKokkos_UnpackCommVel {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_x_array _x;
  typename ArrayTypes<DeviceType>::t_float_1d _radius,_rmass;
  typename ArrayTypes<DeviceType>::t_v_array _v, _omega;
  typename ArrayTypes<DeviceType>::t_xfloat_2d_const _buf;
  int _first;

  AtomVecSphereKokkos_UnpackCommVel(
    const typename DAT::tdual_x_array &x,
    const typename DAT::tdual_float_1d &radius,
    const typename DAT::tdual_float_1d &rmass,
    const typename DAT::tdual_v_array &v,
    const typename DAT::tdual_v_array &omega,
    const typename DAT::tdual_xfloat_2d &buf,
    const int& first):
    _x(x.view<DeviceType>()),
    _radius(radius.view<DeviceType>()),
    _rmass(rmass.view<DeviceType>()),
    _v(v.view<DeviceType>()),
    _omega(omega.view<DeviceType>()),
    _first(first)
  {
    const size_t elements = 9 + 2 * RADVARY;
    const int maxsend = (buf.template view<DeviceType>().extent(0)*buf.template view<DeviceType>().extent(1))/elements;
    buffer_view<DeviceType>(_buf,buf,maxsend,elements);
  };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    _x(i+_first,0) = _buf(i,0);
    _x(i+_first,1) = _buf(i,1);
    _x(i+_first,2) = _buf(i,2);
    _v(i+_first,0) = _buf(i,3);
    _v(i+_first,1) = _buf(i,4);
    _v(i+_first,2) = _buf(i,5);
    _omega(i+_first,0) = _buf(i,6);
    _omega(i+_first,1) = _buf(i,7);
    _omega(i+_first,2) = _buf(i,8);
    if (RADVARY) {
      _radius(i+_first) = _buf(i,9);
      _rmass(i+_first) = _buf(i,10);
    }
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecSphereKokkos::unpack_comm_vel_kokkos(
  const int &n, const int &first,
  const DAT::tdual_xfloat_2d &buf ) {
  if(commKK->forward_comm_on_host) {
    atomKK->modified(Host,X_MASK|RADIUS_MASK|RMASS_MASK|V_MASK|OMEGA_MASK);
    if (radvary == 0) {
      struct AtomVecSphereKokkos_UnpackCommVel<LMPHostType,0> f(
        atomKK->k_x,
        atomKK->k_radius,atomKK->k_rmass,
        atomKK->k_v,atomKK->k_omega,
        buf,first);
      Kokkos::parallel_for(n,f);
    } else {
      struct AtomVecSphereKokkos_UnpackCommVel<LMPHostType,1> f(
        atomKK->k_x,
        atomKK->k_radius,atomKK->k_rmass,
        atomKK->k_v,atomKK->k_omega,
        buf,first);
      Kokkos::parallel_for(n,f);
    }
  } else {
    atomKK->modified(Device,X_MASK|RADIUS_MASK|RMASS_MASK|V_MASK|OMEGA_MASK);
    if (radvary == 0) {
      struct AtomVecSphereKokkos_UnpackCommVel<LMPDeviceType,0> f(
        atomKK->k_x,
        atomKK->k_radius,atomKK->k_rmass,
        atomKK->k_v,atomKK->k_omega,
        buf,first);
      Kokkos::parallel_for(n,f);
    } else {
      struct AtomVecSphereKokkos_UnpackCommVel<LMPDeviceType,1> f(
        atomKK->k_x,
        atomKK->k_radius,atomKK->k_rmass,
        atomKK->k_v,atomKK->k_omega,
        buf,first);
      Kokkos::parallel_for(n,f);
    }
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecSphereKokkos::pack_comm(int n, int *list, double *buf,
			           int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  if (radvary == 0) {
    // Not sure if we need to call sync for X here
    atomKK->sync(Host,X_MASK);
    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = h_x(j,0);
        buf[m++] = h_x(j,1);
        buf[m++] = h_x(j,2);
      }
    } else {
      if (domain->triclinic == 0) {
        dx = pbc[0]*domain->xprd;
        dy = pbc[1]*domain->yprd;
        dz = pbc[2]*domain->zprd;
      } else {
        dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
        dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
        dz = pbc[2]*domain->zprd;
      }
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = h_x(j,0) + dx;
        buf[m++] = h_x(j,1) + dy;
        buf[m++] = h_x(j,2) + dz;
      }
    }
  } else {
    atomKK->sync(Host,X_MASK|RADIUS_MASK|RMASS_MASK);
    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = h_x(j,0);
        buf[m++] = h_x(j,1);
        buf[m++] = h_x(j,2);
        buf[m++] = h_radius[j];
        buf[m++] = h_rmass[j];
      }
    } else {
      if (domain->triclinic == 0) {
        dx = pbc[0]*domain->xprd;
        dy = pbc[1]*domain->yprd;
        dz = pbc[2]*domain->zprd;
      } else {
        dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
        dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
        dz = pbc[2]*domain->zprd;
      }
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = h_x(j,0) + dx;
        buf[m++] = h_x(j,1) + dy;
        buf[m++] = h_x(j,2) + dz;
        buf[m++] = h_radius[j];
        buf[m++] = h_rmass[j];
      }
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSphereKokkos::pack_comm_vel(int n, int *list, double *buf,
				       int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  if (radvary == 0) {
    atomKK->sync(Host,X_MASK|V_MASK|OMEGA_MASK);
    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = h_x(j,0);
        buf[m++] = h_x(j,1);
        buf[m++] = h_x(j,2);
        buf[m++] = h_v(j,0);
        buf[m++] = h_v(j,1);
        buf[m++] = h_v(j,2);
        buf[m++] = h_omega(j,0);
        buf[m++] = h_omega(j,1);
        buf[m++] = h_omega(j,2);
      }
    } else {
      if (domain->triclinic == 0) {
        dx = pbc[0]*domain->xprd;
        dy = pbc[1]*domain->yprd;
        dz = pbc[2]*domain->zprd;
      } else {
        dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
        dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
        dz = pbc[2]*domain->zprd;
      }
      if (!deform_vremap) {
        for (i = 0; i < n; i++) {
          j = list[i];
          buf[m++] = h_x(j,0) + dx;
          buf[m++] = h_x(j,1) + dy;
          buf[m++] = h_x(j,2) + dz;
          buf[m++] = h_v(j,0);
          buf[m++] = h_v(j,1);
          buf[m++] = h_v(j,2);
          buf[m++] = h_omega(j,0);
          buf[m++] = h_omega(j,1);
          buf[m++] = h_omega(j,2);
        }
      } else {
        dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
        dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
        dvz = pbc[2]*h_rate[2];
        for (i = 0; i < n; i++) {
          j = list[i];
          buf[m++] = h_x(j,0) + dx;
          buf[m++] = h_x(j,1) + dy;
          buf[m++] = h_x(j,2) + dz;
         if (mask[i] & deform_groupbit) {
            buf[m++] = h_v(j,0) + dvx;
            buf[m++] = h_v(j,1) + dvy;
            buf[m++] = h_v(j,2) + dvz;
          } else {
            buf[m++] = h_v(j,0);
            buf[m++] = h_v(j,1);
            buf[m++] = h_v(j,2);
          }
          buf[m++] = h_omega(j,0);
          buf[m++] = h_omega(j,1);
          buf[m++] = h_omega(j,2);
        }
      }
    }
  } else {
    atomKK->sync(Host,X_MASK|RADIUS_MASK|RMASS_MASK|V_MASK|OMEGA_MASK);
    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = h_x(j,0);
        buf[m++] = h_x(j,1);
        buf[m++] = h_x(j,2);
        buf[m++] = h_radius[j];
        buf[m++] = h_rmass[j];
        buf[m++] = h_v(j,0);
        buf[m++] = h_v(j,1);
        buf[m++] = h_v(j,2);
        buf[m++] = h_omega(j,0);
        buf[m++] = h_omega(j,1);
        buf[m++] = h_omega(j,2);
      }
    } else {
      if (domain->triclinic == 0) {
        dx = pbc[0]*domain->xprd;
        dy = pbc[1]*domain->yprd;
        dz = pbc[2]*domain->zprd;
      } else {
        dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
        dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
        dz = pbc[2]*domain->zprd;
      }
      if (!deform_vremap) {
        for (i = 0; i < n; i++) {
          j = list[i];
          buf[m++] = h_x(j,0) + dx;
          buf[m++] = h_x(j,1) + dy;
          buf[m++] = h_x(j,2) + dz;
          buf[m++] = h_radius[j];
          buf[m++] = h_rmass[j];
          buf[m++] = h_v(j,0);
          buf[m++] = h_v(j,1);
          buf[m++] = h_v(j,2);
          buf[m++] = h_omega(j,0);
          buf[m++] = h_omega(j,1);
          buf[m++] = h_omega(j,2);
        }
      } else {
        dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
        dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
        dvz = pbc[2]*h_rate[2];
        for (i = 0; i < n; i++) {
          j = list[i];
          buf[m++] = h_x(j,0) + dx;
          buf[m++] = h_x(j,1) + dy;
          buf[m++] = h_x(j,2) + dz;
          buf[m++] = h_radius[j];
          buf[m++] = h_rmass[j];
          if (mask[i] & deform_groupbit) {
            buf[m++] = h_v(j,0) + dvx;
            buf[m++] = h_v(j,1) + dvy;
            buf[m++] = h_v(j,2) + dvz;
          } else {
            buf[m++] = h_v(j,0);
            buf[m++] = h_v(j,1);
            buf[m++] = h_v(j,2);
          }
          buf[m++] = h_omega(j,0);
          buf[m++] = h_omega(j,1);
          buf[m++] = h_omega(j,2);
        }
      }
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSphereKokkos::pack_comm_hybrid(int n, int *list, double *buf)
{
  if (radvary == 0) return 0;

  atomKK->sync(Host,RADIUS_MASK|RMASS_MASK);

  int m = 0;
  for (int i = 0; i < n; i++) {
    const int j = list[i];
    buf[m++] = h_radius[j];
    buf[m++] = h_rmass[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecSphereKokkos::unpack_comm(int n, int first, double *buf)
{
  if (radvary == 0) {
    int m = 0;
    const int last = first + n;
    for (int i = first; i < last; i++) {
      h_x(i,0) = buf[m++];
      h_x(i,1) = buf[m++];
      h_x(i,2) = buf[m++];
    }
    atomKK->modified(Host,X_MASK);
  } else {
    int m = 0;
    const int last = first + n;
    for (int i = first; i < last; i++) {
      h_x(i,0) = buf[m++];
      h_x(i,1) = buf[m++];
      h_x(i,2) = buf[m++];
      h_radius[i] = buf[m++];
      h_rmass[i] = buf[m++];
    }
    atomKK->modified(Host,X_MASK|RADIUS_MASK|RMASS_MASK);
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecSphereKokkos::unpack_comm_vel(int n, int first, double *buf)
{
  if (radvary == 0) {
    int m = 0;
    const int last = first + n;
    for (int i = first; i < last; i++) {
      h_x(i,0) = buf[m++];
      h_x(i,1) = buf[m++];
      h_x(i,2) = buf[m++];
      h_v(i,0) = buf[m++];
      h_v(i,1) = buf[m++];
      h_v(i,2) = buf[m++];
      h_omega(i,0) = buf[m++];
      h_omega(i,1) = buf[m++];
      h_omega(i,2) = buf[m++];
    }
    atomKK->modified(Host,X_MASK|V_MASK|OMEGA_MASK);
  } else {
    int m = 0;
    const int last = first + n;
    for (int i = first; i < last; i++) {
      h_x(i,0) = buf[m++];
      h_x(i,1) = buf[m++];
      h_x(i,2) = buf[m++];
      h_radius[i] = buf[m++];
      h_rmass[i] = buf[m++];
      h_v(i,0) = buf[m++];
      h_v(i,1) = buf[m++];
      h_v(i,2) = buf[m++];
      h_omega(i,0) = buf[m++];
      h_omega(i,1) = buf[m++];
      h_omega(i,2) = buf[m++];
    }
    atomKK->modified(Host,X_MASK|RADIUS_MASK|RMASS_MASK|V_MASK|OMEGA_MASK);
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecSphereKokkos::unpack_comm_hybrid(int n, int first, double *buf)
{
  if (radvary == 0) return 0;

  int m = 0;
  const int last = first + n;
  for (int i = first; i < last; i++) {
    h_radius[i] = buf[m++];
    h_rmass[i] = buf[m++];
  }
  atomKK->modified(Host,RADIUS_MASK|RMASS_MASK);
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSphereKokkos::pack_reverse(int n, int first, double *buf)
{
  if(n > 0)
    atomKK->sync(Host,F_MASK|TORQUE_MASK);

  int m = 0;
  const int last = first + n;
  for (int i = first; i < last; i++) {
    buf[m++] = h_f(i,0);
    buf[m++] = h_f(i,1);
    buf[m++] = h_f(i,2);
    buf[m++] = h_torque(i,0);
    buf[m++] = h_torque(i,1);
    buf[m++] = h_torque(i,2);
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSphereKokkos::pack_reverse_hybrid(int n, int first, double *buf)
{
  if(n > 0)
    atomKK->sync(Host,TORQUE_MASK);

  int m = 0;
  const int last = first + n;
  for (int i = first; i < last; i++) {
    buf[m++] = h_torque(i,0);
    buf[m++] = h_torque(i,1);
    buf[m++] = h_torque(i,2);
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecSphereKokkos::unpack_reverse(int n, int *list, double *buf)
{
  if(n > 0) {
    atomKK->modified(Host,F_MASK|TORQUE_MASK);
  }

  int m = 0;
  for (int i = 0; i < n; i++) {
    const int j = list[i];
    h_f(j,0) += buf[m++];
    h_f(j,1) += buf[m++];
    h_f(j,2) += buf[m++];
    h_torque(j,0) += buf[m++];
    h_torque(j,1) += buf[m++];
    h_torque(j,2) += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecSphereKokkos::unpack_reverse_hybrid(int n, int *list, double *buf)
{
  if(n > 0) {
    atomKK->modified(Host,TORQUE_MASK);
  }

  int m = 0;
  for (int i = 0; i < n; i++) {
    const int j = list[i];
    h_torque(j,0) += buf[m++];
    h_torque(j,1) += buf[m++];
    h_torque(j,2) += buf[m++];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG>
struct AtomVecSphereKokkos_PackBorder {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_xfloat_2d_um _buf;
  const typename ArrayTypes<DeviceType>::t_int_2d_const _list;
  const int _iswap;
  const typename ArrayTypes<DeviceType>::t_x_array_randomread _x;
  const typename ArrayTypes<DeviceType>::t_tagint_1d _tag;
  const typename ArrayTypes<DeviceType>::t_int_1d _type;
  const typename ArrayTypes<DeviceType>::t_int_1d _mask;
  typename ArrayTypes<DeviceType>::t_float_1d _radius,_rmass;
  X_FLOAT _dx,_dy,_dz;

  AtomVecSphereKokkos_PackBorder(
    const typename ArrayTypes<DeviceType>::t_xfloat_2d &buf,
    const typename ArrayTypes<DeviceType>::t_int_2d_const &list,
    const int &iswap,
    const typename ArrayTypes<DeviceType>::t_x_array &x,
    const typename ArrayTypes<DeviceType>::t_tagint_1d &tag,
    const typename ArrayTypes<DeviceType>::t_int_1d &type,
    const typename ArrayTypes<DeviceType>::t_int_1d &mask,
    const typename ArrayTypes<DeviceType>::t_float_1d &radius,
    const typename ArrayTypes<DeviceType>::t_float_1d &rmass,
    const X_FLOAT &dx, const X_FLOAT &dy, const X_FLOAT &dz):
    _list(list),_iswap(iswap),
    _x(x),_tag(tag),_type(type),_mask(mask),
    _radius(radius),
    _rmass(rmass),
    _dx(dx),_dy(dy),_dz(dz)
  {
    const size_t elements = 8;
    const int maxsend = (buf.extent(0)*buf.extent(1))/elements;
    _buf = typename ArrayTypes<DeviceType>::t_xfloat_2d_um(buf.data(),maxsend,elements);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(_iswap,i);
    if (PBC_FLAG == 0) {
      _buf(i,0) = _x(j,0);
      _buf(i,1) = _x(j,1);
      _buf(i,2) = _x(j,2);
    } else {
      _buf(i,0) = _x(j,0) + _dx;
      _buf(i,1) = _x(j,1) + _dy;
      _buf(i,2) = _x(j,2) + _dz;
    }
    _buf(i,3) = d_ubuf(_tag(j)).d;
    _buf(i,4) = d_ubuf(_type(j)).d;
    _buf(i,5) = d_ubuf(_mask(j)).d;
    _buf(i,6) = _radius(j);
    _buf(i,7) = _rmass(j);
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecSphereKokkos::pack_border_kokkos(
  int n, DAT::tdual_int_2d k_sendlist, DAT::tdual_xfloat_2d buf,int iswap,
  int pbc_flag, int *pbc, ExecutionSpace space)
{
  X_FLOAT dx,dy,dz;

  // This was in atom_vec_dpd_kokkos but doesn't appear in any other atom_vec
  atomKK->sync(space,ALL_MASK);

  if (pbc_flag != 0) {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    if(space==Host) {
      AtomVecSphereKokkos_PackBorder<LMPHostType,1> f(
        buf.view<LMPHostType>(), k_sendlist.view<LMPHostType>(),
        iswap,h_x,h_tag,h_type,h_mask,
        h_radius,h_rmass,
        dx,dy,dz);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecSphereKokkos_PackBorder<LMPDeviceType,1> f(
        buf.view<LMPDeviceType>(), k_sendlist.view<LMPDeviceType>(),
        iswap,d_x,d_tag,d_type,d_mask,
        d_radius,d_rmass,
        dx,dy,dz);
      Kokkos::parallel_for(n,f);
    }
  } else {
    dx = dy = dz = 0;
    if(space==Host) {
      AtomVecSphereKokkos_PackBorder<LMPHostType,0> f(
        buf.view<LMPHostType>(), k_sendlist.view<LMPHostType>(),
        iswap,h_x,h_tag,h_type,h_mask,
        h_radius,h_rmass,
        dx,dy,dz);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecSphereKokkos_PackBorder<LMPDeviceType,0> f(
        buf.view<LMPDeviceType>(), k_sendlist.view<LMPDeviceType>(),
        iswap,d_x,d_tag,d_type,d_mask,
        d_radius,d_rmass,
        dx,dy,dz);
      Kokkos::parallel_for(n,f);
    }
  }
  return n*size_border;
}

/* ---------------------------------------------------------------------- */

int AtomVecSphereKokkos::pack_border(
  int n, int *list, double *buf,
  int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  atomKK->sync(Host,ALL_MASK);

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = h_x(j,0);
      buf[m++] = h_x(j,1);
      buf[m++] = h_x(j,2);
      buf[m++] = ubuf(h_tag[j]).d;
      buf[m++] = ubuf(h_type[j]).d;
      buf[m++] = ubuf(h_mask[j]).d;
      buf[m++] = h_radius[j];
      buf[m++] = h_rmass[j];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = h_x(j,0) + dx;
      buf[m++] = h_x(j,1) + dy;
      buf[m++] = h_x(j,2) + dz;
      buf[m++] = ubuf(h_tag[j]).d;
      buf[m++] = ubuf(h_type[j]).d;
      buf[m++] = ubuf(h_mask[j]).d;
      buf[m++] = h_radius[j];
      buf[m++] = h_rmass[j];
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG,int DEFORM_VREMAP>
struct AtomVecSphereKokkos_PackBorderVel {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_xfloat_2d_um _buf;
  const typename ArrayTypes<DeviceType>::t_int_2d_const _list;
  const int _iswap;
  const typename ArrayTypes<DeviceType>::t_x_array_randomread _x;
  const typename ArrayTypes<DeviceType>::t_tagint_1d _tag;
  const typename ArrayTypes<DeviceType>::t_int_1d _type;
  const typename ArrayTypes<DeviceType>::t_int_1d _mask;
  typename ArrayTypes<DeviceType>::t_float_1d _radius,_rmass;
  typename ArrayTypes<DeviceType>::t_v_array _v, _omega;
  X_FLOAT _dx,_dy,_dz, _dvx, _dvy, _dvz;
  const int _deform_groupbit;

  AtomVecSphereKokkos_PackBorderVel(
    const typename ArrayTypes<DeviceType>::t_xfloat_2d &buf,
    const typename ArrayTypes<DeviceType>::t_int_2d_const &list,
    const int &iswap,
    const typename ArrayTypes<DeviceType>::t_x_array &x,
    const typename ArrayTypes<DeviceType>::t_tagint_1d &tag,
    const typename ArrayTypes<DeviceType>::t_int_1d &type,
    const typename ArrayTypes<DeviceType>::t_int_1d &mask,
    const typename ArrayTypes<DeviceType>::t_float_1d &radius,
    const typename ArrayTypes<DeviceType>::t_float_1d &rmass,
    const typename ArrayTypes<DeviceType>::t_v_array &v,
    const typename ArrayTypes<DeviceType>::t_v_array &omega,
    const X_FLOAT &dx, const X_FLOAT &dy, const X_FLOAT &dz,
    const X_FLOAT &dvx, const X_FLOAT &dvy, const X_FLOAT &dvz,
    const int &deform_groupbit):
    _buf(buf),_list(list),_iswap(iswap),
    _x(x),_tag(tag),_type(type),_mask(mask),
    _radius(radius),
    _rmass(rmass),
    _v(v), _omega(omega),
    _dx(dx),_dy(dy),_dz(dz),
    _dvx(dvx),_dvy(dvy),_dvz(dvz),
    _deform_groupbit(deform_groupbit)
  {
    const size_t elements = 14;
    const int maxsend = (buf.extent(0)*buf.extent(1))/elements;
    _buf = typename ArrayTypes<DeviceType>::t_xfloat_2d_um(buf.data(),maxsend,elements);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(_iswap,i);
    if (PBC_FLAG == 0) {
      _buf(i,0) = _x(j,0);
      _buf(i,1) = _x(j,1);
      _buf(i,2) = _x(j,2);
    } else {
      _buf(i,0) = _x(j,0) + _dx;
      _buf(i,1) = _x(j,1) + _dy;
      _buf(i,2) = _x(j,2) + _dz;
    }
    _buf(i,3) = d_ubuf(_tag(j)).d;
    _buf(i,4) = d_ubuf(_type(j)).d;
    _buf(i,5) = d_ubuf(_mask(j)).d;
    _buf(i,6) = _radius(j);
    _buf(i,7) = _rmass(j);
    if (DEFORM_VREMAP) {
      if (_mask(i) & _deform_groupbit) {
        _buf(i,8) = _v(j,0) + _dvx;
        _buf(i,9) = _v(j,1) + _dvy;
        _buf(i,10) = _v(j,2) + _dvz;
      }
    }
    else {
      _buf(i,8) = _v(j,0);
      _buf(i,9) = _v(j,1);
      _buf(i,10) = _v(j,2);
    }
    _buf(i,11) = _omega(j,0);
    _buf(i,12) = _omega(j,1);
    _buf(i,13) = _omega(j,2);
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecSphereKokkos::pack_border_vel_kokkos(
  int n, DAT::tdual_int_2d k_sendlist, DAT::tdual_xfloat_2d buf,int iswap,
  int pbc_flag, int *pbc, ExecutionSpace space)
{
  X_FLOAT dx=0,dy=0,dz=0;
  X_FLOAT dvx=0,dvy=0,dvz=0;

  // This was in atom_vec_dpd_kokkos but doesn't appear in any other atom_vec
  atomKK->sync(space,ALL_MASK);

  if (pbc_flag != 0) {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    if (!deform_vremap) {
      if(space==Host) {
        AtomVecSphereKokkos_PackBorderVel<LMPHostType,1,0> f(
          buf.view<LMPHostType>(), k_sendlist.view<LMPHostType>(),
          iswap,h_x,h_tag,h_type,h_mask,
          h_radius,h_rmass,
          h_v, h_omega,
          dx,dy,dz,dvx,dvy,dvz,
          deform_groupbit);
        Kokkos::parallel_for(n,f);
      } else {
        AtomVecSphereKokkos_PackBorderVel<LMPDeviceType,1,0> f(
          buf.view<LMPDeviceType>(), k_sendlist.view<LMPDeviceType>(),
          iswap,d_x,d_tag,d_type,d_mask,
          d_radius,d_rmass,
          d_v, d_omega,
          dx,dy,dz,dvx,dvy,dvz,
          deform_groupbit);
        Kokkos::parallel_for(n,f);
      }
    }
    else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      if(space==Host) {
        AtomVecSphereKokkos_PackBorderVel<LMPHostType,1,1> f(
          buf.view<LMPHostType>(), k_sendlist.view<LMPHostType>(),
          iswap,h_x,h_tag,h_type,h_mask,
          h_radius,h_rmass,
          h_v, h_omega,
          dx,dy,dz,dvx,dvy,dvz,
          deform_groupbit);
        Kokkos::parallel_for(n,f);
      } else {
        AtomVecSphereKokkos_PackBorderVel<LMPDeviceType,1,1> f(
          buf.view<LMPDeviceType>(), k_sendlist.view<LMPDeviceType>(),
          iswap,d_x,d_tag,d_type,d_mask,
          d_radius,d_rmass,
          d_v, d_omega,
          dx,dy,dz,dvx,dvy,dvz,
          deform_groupbit);
        Kokkos::parallel_for(n,f);
      }
    }
  } else {
    if(space==Host) {
      AtomVecSphereKokkos_PackBorderVel<LMPHostType,0,0> f(
        buf.view<LMPHostType>(), k_sendlist.view<LMPHostType>(),
        iswap,h_x,h_tag,h_type,h_mask,
        h_radius,h_rmass,
        h_v, h_omega,
        dx,dy,dz,dvx,dvy,dvz,
        deform_groupbit);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecSphereKokkos_PackBorderVel<LMPDeviceType,0,0> f(
        buf.view<LMPDeviceType>(), k_sendlist.view<LMPDeviceType>(),
        iswap,d_x,d_tag,d_type,d_mask,
        d_radius,d_rmass,
        d_v, d_omega,
        dx,dy,dz,dvx,dvy,dvz,
        deform_groupbit);
      Kokkos::parallel_for(n,f);
    }
  }

  return n*(size_border + size_velocity);
}

/* ---------------------------------------------------------------------- */

int AtomVecSphereKokkos::pack_border_vel(int n, int *list, double *buf,
                                         int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  atomKK->sync(Host,ALL_MASK);

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = h_x(j,0);
      buf[m++] = h_x(j,1);
      buf[m++] = h_x(j,2);
      buf[m++] = ubuf(h_tag[j]).d;
      buf[m++] = ubuf(h_type[j]).d;
      buf[m++] = ubuf(h_mask[j]).d;
      buf[m++] = h_radius[j];
      buf[m++] = h_rmass[j];
      buf[m++] = h_v(j,0);
      buf[m++] = h_v(j,1);
      buf[m++] = h_v(j,2);
      buf[m++] = h_omega(j,0);
      buf[m++] = h_omega(j,1);
      buf[m++] = h_omega(j,2);
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = h_x(j,0) + dx;
        buf[m++] = h_x(j,1) + dy;
        buf[m++] = h_x(j,2) + dz;
        buf[m++] = ubuf(h_tag[j]).d;
        buf[m++] = ubuf(h_type[j]).d;
        buf[m++] = ubuf(h_mask[j]).d;
        buf[m++] = h_radius[j];
        buf[m++] = h_rmass[j];
        buf[m++] = h_v(j,0);
        buf[m++] = h_v(j,1);
        buf[m++] = h_v(j,2);
        buf[m++] = h_omega(j,0);
        buf[m++] = h_omega(j,1);
        buf[m++] = h_omega(j,2);
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = h_x(j,0) + dx;
        buf[m++] = h_x(j,1) + dy;
        buf[m++] = h_x(j,2) + dz;
        buf[m++] = ubuf(h_tag[j]).d;
        buf[m++] = ubuf(h_type[j]).d;
        buf[m++] = ubuf(h_mask[j]).d;
        buf[m++] = h_radius[j];
        buf[m++] = h_rmass[j];
        if (mask[i] & deform_groupbit) {
          buf[m++] = h_v(j,0) + dvx;
          buf[m++] = h_v(j,1) + dvy;
          buf[m++] = h_v(j,2) + dvz;
        } else {
          buf[m++] = h_v(j,0);
          buf[m++] = h_v(j,1);
          buf[m++] = h_v(j,2);
        }
        buf[m++] = h_omega(j,0);
        buf[m++] = h_omega(j,1);
        buf[m++] = h_omega(j,2);
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSphereKokkos::pack_border_hybrid(int n, int *list, double *buf)
{
  atomKK->sync(Host,RADIUS_MASK|RMASS_MASK);

  int m = 0;
  for (int i = 0; i < n; i++) {
    const int j = list[i];
    buf[m++] = h_radius[j];
    buf[m++] = h_rmass[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecSphereKokkos_UnpackBorder {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_xfloat_2d_const_um _buf;
  typename ArrayTypes<DeviceType>::t_x_array _x;
  typename ArrayTypes<DeviceType>::t_tagint_1d _tag;
  typename ArrayTypes<DeviceType>::t_int_1d _type;
  typename ArrayTypes<DeviceType>::t_int_1d _mask;
  typename ArrayTypes<DeviceType>::t_float_1d _radius,_rmass;
  int _first;

  AtomVecSphereKokkos_UnpackBorder(
    const typename ArrayTypes<DeviceType>::t_xfloat_2d &buf,
    const typename ArrayTypes<DeviceType>::t_x_array &x,
    const typename ArrayTypes<DeviceType>::t_tagint_1d &tag,
    const typename ArrayTypes<DeviceType>::t_int_1d &type,
    const typename ArrayTypes<DeviceType>::t_int_1d &mask,
    const typename ArrayTypes<DeviceType>::t_float_1d &radius,
    const typename ArrayTypes<DeviceType>::t_float_1d &rmass,
    const int& first):
    _buf(buf),_x(x),_tag(tag),_type(type),_mask(mask),
    _radius(radius),
    _rmass(rmass),
    _first(first)
  {
    const size_t elements = 8;
    const int maxsend = (buf.extent(0)*buf.extent(1))/elements;
    _buf = typename ArrayTypes<DeviceType>::t_xfloat_2d_const_um(buf.data(),maxsend,elements);
  };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    _x(i+_first,0) = _buf(i,0);
    _x(i+_first,1) = _buf(i,1);
    _x(i+_first,2) = _buf(i,2);
    _tag(i+_first) = static_cast<tagint> (d_ubuf(_buf(i,3)).i);
    _type(i+_first) = static_cast<int>  (d_ubuf(_buf(i,4)).i);
    _mask(i+_first) = static_cast<int>  (d_ubuf(_buf(i,5)).i);
    _radius(i+_first) = _buf(i,6);
    _rmass(i+_first) = _buf(i,7);
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecSphereKokkos::unpack_border_kokkos(const int &n, const int &first,
					       const DAT::tdual_xfloat_2d &buf,ExecutionSpace space) {
  while (first+n >= nmax) grow(0);
  if(space==Host) {
    struct AtomVecSphereKokkos_UnpackBorder<LMPHostType> f(buf.view<LMPHostType>(),
      h_x,h_tag,h_type,h_mask,
      h_radius,h_rmass,
      first);
    Kokkos::parallel_for(n,f);
  } else {
    struct AtomVecSphereKokkos_UnpackBorder<LMPDeviceType> f(buf.view<LMPDeviceType>(),
      d_x,d_tag,d_type,d_mask,
      d_radius,d_rmass,
      first);
    Kokkos::parallel_for(n,f);
  }

  atomKK->modified(space,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|
	         RADIUS_MASK|RMASS_MASK);
}

/* ---------------------------------------------------------------------- */

void AtomVecSphereKokkos::unpack_border(int n, int first, double *buf)
{
  int m = 0;
  const int last = first + n;
  for (int i = first; i < last; i++) {
    if (i == nmax) grow(0);
    h_x(i,0) = buf[m++];
    h_x(i,1) = buf[m++];
    h_x(i,2) = buf[m++];
    h_tag[i] = (tagint) ubuf(buf[m++]).i;
    h_type[i] = (int) ubuf(buf[m++]).i;
    h_mask[i] = (int) ubuf(buf[m++]).i;
    h_radius[i] = buf[m++];
    h_rmass[i] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);

  atomKK->modified(Host,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|RADIUS_MASK|RMASS_MASK);
}


/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecSphereKokkos_UnpackBorderVel {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_xfloat_2d_const_um _buf;
  typename ArrayTypes<DeviceType>::t_x_array _x;
  typename ArrayTypes<DeviceType>::t_tagint_1d _tag;
  typename ArrayTypes<DeviceType>::t_int_1d _type;
  typename ArrayTypes<DeviceType>::t_int_1d _mask;
  typename ArrayTypes<DeviceType>::t_float_1d _radius,_rmass;
  typename ArrayTypes<DeviceType>::t_v_array _v;
  typename ArrayTypes<DeviceType>::t_v_array _omega;
  int _first;

  AtomVecSphereKokkos_UnpackBorderVel(
    const typename ArrayTypes<DeviceType>::t_xfloat_2d_const &buf,
    const typename ArrayTypes<DeviceType>::t_x_array &x,
    const typename ArrayTypes<DeviceType>::t_tagint_1d &tag,
    const typename ArrayTypes<DeviceType>::t_int_1d &type,
    const typename ArrayTypes<DeviceType>::t_int_1d &mask,
    const typename ArrayTypes<DeviceType>::t_float_1d &radius,
    const typename ArrayTypes<DeviceType>::t_float_1d &rmass,
    const typename ArrayTypes<DeviceType>::t_v_array &v,
    const typename ArrayTypes<DeviceType>::t_v_array &omega,
    const int& first):
    _buf(buf),_x(x),_tag(tag),_type(type),_mask(mask),
    _radius(radius),
    _rmass(rmass),
    _v(v), _omega(omega),
    _first(first)
  {
    const size_t elements = 14;
    const int maxsend = (buf.extent(0)*buf.extent(1))/elements;
    _buf = typename ArrayTypes<DeviceType>::t_xfloat_2d_const_um(buf.data(),maxsend,elements);
  };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    _x(i+_first,0) = _buf(i,0);
    _x(i+_first,1) = _buf(i,1);
    _x(i+_first,2) = _buf(i,2);
    _tag(i+_first) = static_cast<tagint> (d_ubuf(_buf(i,3)).i);
    _type(i+_first) = static_cast<int>  (d_ubuf(_buf(i,4)).i);
    _mask(i+_first) = static_cast<int>  (d_ubuf(_buf(i,5)).i);
    _radius(i+_first) = _buf(i,6);
    _rmass(i+_first) = _buf(i,7);
    _v(i+_first,0) = _buf(i,8);
    _v(i+_first,1) = _buf(i,9);
    _v(i+_first,2) = _buf(i,10);
    _omega(i+_first,0) = _buf(i,11);
    _omega(i+_first,1) = _buf(i,12);
    _omega(i+_first,2) = _buf(i,13);
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecSphereKokkos::unpack_border_vel_kokkos(
  const int &n, const int &first,
  const DAT::tdual_xfloat_2d &buf,ExecutionSpace space) {
  while (first+n >= nmax) grow(0);
  if(space==Host) {
    struct AtomVecSphereKokkos_UnpackBorderVel<LMPHostType> f(buf.view<LMPHostType>(),
      h_x,h_tag,h_type,h_mask,
      h_radius,h_rmass,
      h_v, h_omega,
      first);
    Kokkos::parallel_for(n,f);
  } else {
    struct AtomVecSphereKokkos_UnpackBorderVel<LMPDeviceType> f(buf.view<LMPDeviceType>(),
      d_x,d_tag,d_type,d_mask,
      d_radius,d_rmass,
      d_v, d_omega,
      first);
    Kokkos::parallel_for(n,f);
  }

  atomKK->modified(space,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|
	         RADIUS_MASK|RMASS_MASK|V_MASK|OMEGA_MASK);
}

/* ---------------------------------------------------------------------- */

void AtomVecSphereKokkos::unpack_border_vel(int n, int first, double *buf)
{
  int m = 0;
  const int last = first + n;
  for (int i = first; i < last; i++) {
    if (i == nmax) grow(0);
    h_x(i,0) = buf[m++];
    h_x(i,1) = buf[m++];
    h_x(i,2) = buf[m++];
    h_tag[i] = (tagint) ubuf(buf[m++]).i;
    h_type[i] = (int) ubuf(buf[m++]).i;
    h_mask[i] = (int) ubuf(buf[m++]).i;
    h_radius[i] = buf[m++];
    h_rmass[i] = buf[m++];
    h_v(i,0) = buf[m++];
    h_v(i,1) = buf[m++];
    h_v(i,2) = buf[m++];
    h_omega(i,0) = buf[m++];
    h_omega(i,1) = buf[m++];
    h_omega(i,2) = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);

  atomKK->modified(Host,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|RADIUS_MASK|RMASS_MASK|V_MASK|OMEGA_MASK);
}

/* ---------------------------------------------------------------------- */

int AtomVecSphereKokkos::unpack_border_hybrid(int n, int first, double *buf)
{
  int m = 0;
  const int last = first + n;
  for (int i = first; i < last; i++) {
    h_radius[i] = buf[m++];
    h_rmass[i] = buf[m++];
  }
  atomKK->modified(Host,RADIUS_MASK|RMASS_MASK);
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecSphereKokkos_PackExchangeFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_x_array_randomread _x;
  typename AT::t_v_array_randomread _v;
  typename AT::t_tagint_1d_randomread _tag;
  typename AT::t_int_1d_randomread _type;
  typename AT::t_int_1d_randomread _mask;
  typename AT::t_imageint_1d_randomread _image;
  typename AT::t_float_1d_randomread _radius,_rmass;
  typename AT::t_v_array_randomread _omega;
  typename AT::t_x_array _xw;
  typename AT::t_v_array _vw;
  typename AT::t_tagint_1d _tagw;
  typename AT::t_int_1d _typew;
  typename AT::t_int_1d _maskw;
  typename AT::t_imageint_1d _imagew;
  typename AT::t_float_1d _radiusw,_rmassw;
  typename AT::t_v_array _omegaw;
  typename AT::t_xfloat_2d_um _buf;
  typename AT::t_int_1d_const _sendlist;
  typename AT::t_int_1d_const _copylist;
  int _nlocal,_dim;
  X_FLOAT _lo,_hi;

  AtomVecSphereKokkos_PackExchangeFunctor(
    const AtomKokkos* atom,
    const typename AT::tdual_xfloat_2d buf,
    typename AT::tdual_int_1d sendlist,
    typename AT::tdual_int_1d copylist,int nlocal, int dim,X_FLOAT lo, X_FLOAT hi):
    _x(atom->k_x.view<DeviceType>()),
    _v(atom->k_v.view<DeviceType>()),
    _tag(atom->k_tag.view<DeviceType>()),
    _type(atom->k_type.view<DeviceType>()),
    _mask(atom->k_mask.view<DeviceType>()),
    _image(atom->k_image.view<DeviceType>()),
    _radius(atom->k_radius.view<DeviceType>()),
    _rmass(atom->k_rmass.view<DeviceType>()),
    _omega(atom->k_omega.view<DeviceType>()),
    _xw(atom->k_x.view<DeviceType>()),
    _vw(atom->k_v.view<DeviceType>()),
    _tagw(atom->k_tag.view<DeviceType>()),
    _typew(atom->k_type.view<DeviceType>()),
    _maskw(atom->k_mask.view<DeviceType>()),
    _imagew(atom->k_image.view<DeviceType>()),
    _radiusw(atom->k_radius.view<DeviceType>()),
    _rmassw(atom->k_rmass.view<DeviceType>()),
    _omegaw(atom->k_omega.view<DeviceType>()),
    _sendlist(sendlist.template view<DeviceType>()),
    _copylist(copylist.template view<DeviceType>()),
    _nlocal(nlocal),_dim(dim),
    _lo(lo),_hi(hi)
  {
    const size_t elements = 16;
    const int maxsend = (buf.template view<DeviceType>().extent(0)*buf.template view<DeviceType>().extent(1))/elements;

    _buf = typename AT::t_xfloat_2d_um(buf.template view<DeviceType>().data(),maxsend,elements);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &mysend) const {
    const int i = _sendlist(mysend);
    _buf(mysend,0) = 16;
    _buf(mysend,1) = _x(i,0);
    _buf(mysend,2) = _x(i,1);
    _buf(mysend,3) = _x(i,2);
    _buf(mysend,4) = _v(i,0);
    _buf(mysend,5) = _v(i,1);
    _buf(mysend,6) = _v(i,2);
    _buf(mysend,7) = d_ubuf(_tag[i]).d;
    _buf(mysend,8) = d_ubuf(_type[i]).d;
    _buf(mysend,9) = d_ubuf(_mask[i]).d;
    _buf(mysend,10) = d_ubuf(_image[i]).d;
    _buf(mysend,11) = _radius[i];
    _buf(mysend,12) = _rmass[i];
    _buf(mysend,13) = _omega(i,0);
    _buf(mysend,14) = _omega(i,1);
    _buf(mysend,15) = _omega(i,2);
    const int j = _copylist(mysend);

    if (j>-1) {
      _xw(i,0) = _x(j,0);
      _xw(i,1) = _x(j,1);
      _xw(i,2) = _x(j,2);
      _vw(i,0) = _v(j,0);
      _vw(i,1) = _v(j,1);
      _vw(i,2) = _v(j,2);
      _tagw[i] = _tag(j);
      _typew[i] = _type(j);
      _maskw[i] = _mask(j);
      _imagew[i] = _image(j);
      _radiusw[i] = _radius(j);
      _rmassw[i] = _rmass(j);
      _omegaw(i,0) = _omega(j,0);
      _omegaw(i,1) = _omega(j,1);
      _omegaw(i,2) = _omega(j,2);
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecSphereKokkos::pack_exchange_kokkos(
  const int &nsend,
  DAT::tdual_xfloat_2d &k_buf,
  DAT::tdual_int_1d k_sendlist,
  DAT::tdual_int_1d k_copylist,
  ExecutionSpace space,int dim,X_FLOAT lo,X_FLOAT hi)
{
  if(nsend > (int) (k_buf.view<LMPHostType>().extent(0)*k_buf.view<LMPHostType>().extent(1))/16) {
    int newsize = nsend*17/k_buf.view<LMPHostType>().extent(1)+1;
    k_buf.resize(newsize,k_buf.view<LMPHostType>().extent(1));
  }
  atomKK->sync(space,X_MASK | V_MASK | TAG_MASK | TYPE_MASK |
             MASK_MASK | IMAGE_MASK| RADIUS_MASK | RMASS_MASK |
             OMEGA_MASK);

  if(space == Host) {
    AtomVecSphereKokkos_PackExchangeFunctor<LMPHostType> f(atomKK,k_buf,k_sendlist,k_copylist,atom->nlocal,dim,lo,hi);
    Kokkos::parallel_for(nsend,f);
  } else {
    AtomVecSphereKokkos_PackExchangeFunctor<LMPDeviceType> f(atomKK,k_buf,k_sendlist,k_copylist,atom->nlocal,dim,lo,hi);
    Kokkos::parallel_for(nsend,f);
  }
  return nsend*16;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecSphereKokkos::pack_exchange(int i, double *buf)
{
  atomKK->sync(Host,X_MASK | V_MASK | TAG_MASK | TYPE_MASK |
            MASK_MASK | IMAGE_MASK| RADIUS_MASK | RMASS_MASK |
            OMEGA_MASK);


  int m = 1;
  buf[m++] = h_x(i,0);
  buf[m++] = h_x(i,1);
  buf[m++] = h_x(i,2);
  buf[m++] = h_v(i,0);
  buf[m++] = h_v(i,1);
  buf[m++] = h_v(i,2);
  buf[m++] = ubuf(h_tag[i]).d;
  buf[m++] = ubuf(h_type[i]).d;
  buf[m++] = ubuf(h_mask[i]).d;
  buf[m++] = ubuf(h_image[i]).d;

  buf[m++] = h_radius[i];
  buf[m++] = h_rmass[i];
  buf[m++] = h_omega(i,0);
  buf[m++] = h_omega(i,1);
  buf[m++] = h_omega(i,2);

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecSphereKokkos_UnpackExchangeFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_x_array _x;
  typename AT::t_v_array _v;
  typename AT::t_tagint_1d _tag;
  typename AT::t_int_1d _type;
  typename AT::t_int_1d _mask;
  typename AT::t_imageint_1d _image;
  typename AT::t_float_1d _radius;
  typename AT::t_float_1d _rmass;
  typename AT::t_v_array _omega;
  typename AT::t_xfloat_2d_um _buf;
  typename AT::t_int_1d _nlocal;
  int _dim;
  X_FLOAT _lo,_hi;

  AtomVecSphereKokkos_UnpackExchangeFunctor(
    const AtomKokkos* atom,
    const typename AT::tdual_xfloat_2d buf,
    typename AT::tdual_int_1d nlocal,
    int dim, X_FLOAT lo, X_FLOAT hi):
    _x(atom->k_x.view<DeviceType>()),
    _v(atom->k_v.view<DeviceType>()),
    _tag(atom->k_tag.view<DeviceType>()),
    _type(atom->k_type.view<DeviceType>()),
    _mask(atom->k_mask.view<DeviceType>()),
    _image(atom->k_image.view<DeviceType>()),
    _radius(atom->k_radius.view<DeviceType>()),
    _rmass(atom->k_rmass.view<DeviceType>()),
    _omega(atom->k_omega.view<DeviceType>()),
    _nlocal(nlocal.template view<DeviceType>()),_dim(dim),
    _lo(lo),_hi(hi)
  {
    const size_t elements = 16;
    const int maxsendlist = (buf.template view<DeviceType>().extent(0)*buf.template view<DeviceType>().extent(1))/elements;

    buffer_view<DeviceType>(_buf,buf,maxsendlist,elements);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &myrecv) const {
    X_FLOAT x = _buf(myrecv,_dim+1);
    if (x >= _lo && x < _hi) {
      int i = Kokkos::atomic_fetch_add(&_nlocal(0),1);
      _x(i,0) = _buf(myrecv,1);
      _x(i,1) = _buf(myrecv,2);
      _x(i,2) = _buf(myrecv,3);
      _v(i,0) = _buf(myrecv,4);
      _v(i,1) = _buf(myrecv,5);
      _v(i,2) = _buf(myrecv,6);
      _tag[i] = (tagint) d_ubuf(_buf(myrecv,7)).i;
      _type[i] = (int) d_ubuf(_buf(myrecv,8)).i;
      _mask[i] = (int) d_ubuf(_buf(myrecv,9)).i;
      _image[i] = (imageint) d_ubuf(_buf(myrecv,10)).i;
      _radius[i] = _buf(myrecv,11);
      _rmass[i] = _buf(myrecv,12);
      _omega(i,0) = _buf(myrecv,13);
      _omega(i,1) = _buf(myrecv,14);
      _omega(i,2) = _buf(myrecv,15);
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecSphereKokkos::unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf,int nrecv,int nlocal,int dim,X_FLOAT lo,X_FLOAT hi,ExecutionSpace space) {
  if(space == Host) {
    k_count.h_view(0) = nlocal;
    AtomVecSphereKokkos_UnpackExchangeFunctor<LMPHostType> f(atomKK,k_buf,k_count,dim,lo,hi);
    Kokkos::parallel_for(nrecv/16,f);
  } else {
    k_count.h_view(0) = nlocal;
    k_count.modify<LMPHostType>();
    k_count.sync<LMPDeviceType>();
    AtomVecSphereKokkos_UnpackExchangeFunctor<LMPDeviceType> f(atomKK,k_buf,k_count,dim,lo,hi);
    Kokkos::parallel_for(nrecv/16,f);
    k_count.modify<LMPDeviceType>();
    k_count.sync<LMPHostType>();
  }

  atomKK->modified(space,X_MASK | V_MASK | TAG_MASK | TYPE_MASK |
                 MASK_MASK | IMAGE_MASK| RADIUS_MASK | RMASS_MASK |
                 OMEGA_MASK);

  return k_count.h_view(0);
}

/* ---------------------------------------------------------------------- */

int AtomVecSphereKokkos::unpack_exchange(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  int m = 1;
  h_x(nlocal,0) = buf[m++];
  h_x(nlocal,1) = buf[m++];
  h_x(nlocal,2) = buf[m++];
  h_v(nlocal,0) = buf[m++];
  h_v(nlocal,1) = buf[m++];
  h_v(nlocal,2) = buf[m++];
  h_tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  h_type[nlocal] = (int) ubuf(buf[m++]).i;
  h_mask[nlocal] = (int) ubuf(buf[m++]).i;
  h_image[nlocal] = (imageint) ubuf(buf[m++]).i;

  h_radius[nlocal] = buf[m++];
  h_rmass[nlocal] = buf[m++];
  h_omega(nlocal,0) = buf[m++];
  h_omega(nlocal,1) = buf[m++];
  h_omega(nlocal,2) = buf[m++];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->
        unpack_exchange(nlocal,&buf[m]);

  atomKK->modified(Host,X_MASK | V_MASK | TAG_MASK | TYPE_MASK |
           MASK_MASK | IMAGE_MASK | RADIUS_MASK | RMASS_MASK |
           OMEGA_MASK);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecSphereKokkos::size_restart()
{
  int i;

  int nlocal = atom->nlocal;
  int n = 16 * nlocal;

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
        n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVecSphereKokkos::pack_restart(int i, double *buf)
{
  atomKK->sync(Host,X_MASK | TAG_MASK | TYPE_MASK |
            MASK_MASK | IMAGE_MASK | V_MASK |
            RADIUS_MASK | RMASS_MASK | OMEGA_MASK);

  int m = 1;
  buf[m++] = h_x(i,0);
  buf[m++] = h_x(i,1);
  buf[m++] = h_x(i,2);
  buf[m++] = ubuf(h_tag[i]).d;
  buf[m++] = ubuf(h_type[i]).d;
  buf[m++] = ubuf(h_mask[i]).d;
  buf[m++] = ubuf(h_image[i]).d;
  buf[m++] = h_v(i,0);
  buf[m++] = h_v(i,1);
  buf[m++] = h_v(i,2);

  buf[m++] = h_radius[i];
  buf[m++] = h_rmass[i];
  buf[m++] = h_omega(i,0);
  buf[m++] = h_omega(i,1);
  buf[m++] = h_omega(i,2);

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecSphereKokkos::unpack_restart(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra,nmax,atom->nextra_store,"atom:extra");
  }

  int m = 1;
  h_x(nlocal,0) = buf[m++];
  h_x(nlocal,1) = buf[m++];
  h_x(nlocal,2) = buf[m++];
  h_tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  h_type[nlocal] = (int) ubuf(buf[m++]).i;
  h_mask[nlocal] = (int) ubuf(buf[m++]).i;
  h_image[nlocal] = (imageint) ubuf(buf[m++]).i;
  h_v(nlocal,0) = buf[m++];
  h_v(nlocal,1) = buf[m++];
  h_v(nlocal,2) = buf[m++];

  h_radius[nlocal] = buf[m++];
  h_rmass[nlocal] = buf[m++];
  h_omega(nlocal,0) = buf[m++];
  h_omega(nlocal,1) = buf[m++];
  h_omega(nlocal,2) = buf[m++];

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atomKK->modified(Host,X_MASK | TAG_MASK | TYPE_MASK |
                MASK_MASK | IMAGE_MASK | V_MASK |
	        RADIUS_MASK | RMASS_MASK | OMEGA_MASK);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecSphereKokkos::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
  }

  h_tag[nlocal] = 0;
  h_type[nlocal] = itype;
  h_x(nlocal,0) = coord[0];
  h_x(nlocal,1) = coord[1];
  h_x(nlocal,2) = coord[2];
  h_mask[nlocal] = 1;
  h_image[nlocal] = ((imageint) IMGMAX << IMG2BITS) |
    ((imageint) IMGMAX << IMGBITS) | IMGMAX;
  h_v(nlocal,0) = 0.0;
  h_v(nlocal,1) = 0.0;
  h_v(nlocal,2) = 0.0;

  h_radius[nlocal] = 0.5;
  h_rmass[nlocal] = 4.0*MY_PI/3.0 * h_radius[nlocal]*h_radius[nlocal]*h_radius[nlocal];
  h_omega(nlocal,0) = 0.0;
  h_omega(nlocal,1) = 0.0;
  h_omega(nlocal,2) = 0.0;

  atomKK->modified(Host,ALL_MASK);

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecSphereKokkos::data_atom(double *coord, imageint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = ATOTAGINT(values[0]);
  type[nlocal] = atoi(values[1]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  radius[nlocal] = 0.5 * atof(values[2]);
  if (radius[nlocal] < 0.0)
    error->one(FLERR,"Invalid radius in Atoms section of data file");

  double density = atof(values[3]);
  if (density <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  if (radius[nlocal] == 0.0) rmass[nlocal] = density;
  else
    rmass[nlocal] = 4.0*MY_PI/3.0 *
      radius[nlocal]*radius[nlocal]*radius[nlocal] * density;

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  omega[nlocal][0] = 0.0;
  omega[nlocal][1] = 0.0;
  omega[nlocal][2] = 0.0;

  atomKK->modified(Host,ALL_MASK);

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecSphereKokkos::data_atom_hybrid(int nlocal, char **values)
{
  radius[nlocal] = 0.5 * atof(values[0]);
  if (radius[nlocal] < 0.0)
    error->one(FLERR,"Invalid radius in Atoms section of data file");

  double density = atof(values[1]);
  if (density <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  if (radius[nlocal] == 0.0) rmass[nlocal] = density;
  else
    rmass[nlocal] = 4.0*MY_PI/3.0 *
      radius[nlocal]*radius[nlocal]*radius[nlocal] * density;


  atomKK->modified(Host,RADIUS_MASK|RMASS_MASK);

  return 2;
}

/* ----------------------------------------------------------------------
   unpack one line from Velocities section of data file
------------------------------------------------------------------------- */

void AtomVecSphereKokkos::data_vel(int m, char **values)
{
  atomKK->sync(Host,V_MASK|OMEGA_MASK);
  h_v(m,0) = atof(values[0]);
  h_v(m,1) = atof(values[1]);
  h_v(m,2) = atof(values[2]);
  h_omega(m,0) = atof(values[3]);
  h_omega(m,1) = atof(values[4]);
  h_omega(m,2) = atof(values[5]);
  atomKK->modified(Host,V_MASK|OMEGA_MASK);
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Velocities section of data file
------------------------------------------------------------------------- */

int AtomVecSphereKokkos::data_vel_hybrid(int m, char **values)
{
  atomKK->sync(Host,OMEGA_MASK);
  omega[m][0] = atof(values[0]);
  omega[m][1] = atof(values[1]);
  omega[m][2] = atof(values[2]);
  atomKK->modified(Host,OMEGA_MASK);
  return 3;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecSphereKokkos::pack_data(double **buf)
{
  atomKK->sync(Host,TAG_MASK|TYPE_MASK|RADIUS_MASK|RMASS_MASK|X_MASK|IMAGE_MASK);

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(h_tag[i]).d;
    buf[i][1] = ubuf(h_type[i]).d;
    buf[i][2] = 2.0*h_radius[i];
    if (h_radius[i] == 0.0) buf[i][3] = h_rmass[i];
    else
      buf[i][3] = h_rmass[i] / (4.0*MY_PI/3.0 * h_radius[i]*h_radius[i]*h_radius[i]);
    buf[i][4] = h_x(i,0);
    buf[i][5] = h_x(i,1);
    buf[i][6] = h_x(i,2);
    buf[i][7] = ubuf((h_image[i] & IMGMASK) - IMGMAX).d;
    buf[i][8] = ubuf((h_image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][9] = ubuf((h_image[i] >> IMG2BITS) - IMGMAX).d;
  }
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecSphereKokkos::pack_data_hybrid(int i, double *buf)
{
  atomKK->sync(Host,RADIUS_MASK|RMASS_MASK);

  buf[0] = 2.0*h_radius[i];
  if (h_radius[i] == 0.0) buf[1] = h_rmass[i];
  else buf[1] = h_rmass[i] / (4.0*MY_PI/3.0 * h_radius[i]*h_radius[i]*h_radius[i]);
  return 2;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecSphereKokkos::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT
            " %d %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %d %d %d\n",
            (tagint) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,
            buf[i][2],buf[i][3],
            buf[i][4],buf[i][5],buf[i][6],
            (int) ubuf(buf[i][7]).i,(int) ubuf(buf[i][8]).i,
            (int) ubuf(buf[i][9]).i);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecSphereKokkos::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e %-1.16e",buf[0],buf[1]);
  return 2;
}

/* ----------------------------------------------------------------------
   pack velocity info for data file
------------------------------------------------------------------------- */

void AtomVecSphereKokkos::pack_vel(double **buf)
{
  atomKK->sync(Host,TAG_MASK|V_MASK|OMEGA_MASK);

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(h_tag[i]).d;
    buf[i][1] = h_v(i,0);
    buf[i][2] = h_v(i,1);
    buf[i][3] = h_v(i,2);
    buf[i][4] = h_omega(i,0);
    buf[i][5] = h_omega(i,1);
    buf[i][6] = h_omega(i,2);
  }
}

/* ----------------------------------------------------------------------
   pack hybrid velocity info for data file
------------------------------------------------------------------------- */

int AtomVecSphereKokkos::pack_vel_hybrid(int i, double *buf)
{
  atomKK->sync(Host,OMEGA_MASK);

  buf[0] = h_omega(i,0);
  buf[1] = h_omega(i,1);
  buf[2] = h_omega(i,2);
  return 3;
}

/* ----------------------------------------------------------------------
   write velocity info to data file
------------------------------------------------------------------------- */

void AtomVecSphereKokkos::write_vel(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT
            " %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e\n",
            (tagint) ubuf(buf[i][0]).i,buf[i][1],buf[i][2],buf[i][3],
            buf[i][4],buf[i][5],buf[i][6]);
}

/* ----------------------------------------------------------------------
   write hybrid velocity info to data file
------------------------------------------------------------------------- */

int AtomVecSphereKokkos::write_vel_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e %-1.16e %-1.16e",buf[0],buf[1],buf[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecSphereKokkos::memory_usage()
{
  bigint bytes = 0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);

  if (atom->memcheck("radius")) bytes += memory->usage(radius,nmax);
  if (atom->memcheck("rmass")) bytes += memory->usage(rmass,nmax);
  if (atom->memcheck("omega")) bytes += memory->usage(omega,nmax,3);
  if (atom->memcheck("torque"))
    bytes += memory->usage(torque,nmax*comm->nthreads,3);

  return bytes;
}

/* ---------------------------------------------------------------------- */

void AtomVecSphereKokkos::sync(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.sync<LMPDeviceType>();
    if (mask & V_MASK) atomKK->k_v.sync<LMPDeviceType>();
    if (mask & F_MASK) atomKK->k_f.sync<LMPDeviceType>();
    if (mask & TAG_MASK) atomKK->k_tag.sync<LMPDeviceType>();
    if (mask & TYPE_MASK) atomKK->k_type.sync<LMPDeviceType>();
    if (mask & MASK_MASK) atomKK->k_mask.sync<LMPDeviceType>();
    if (mask & IMAGE_MASK) atomKK->k_image.sync<LMPDeviceType>();
    if (mask & RADIUS_MASK) atomKK->k_radius.sync<LMPDeviceType>();
    if (mask & RMASS_MASK) atomKK->k_rmass.sync<LMPDeviceType>();
    if (mask & OMEGA_MASK) atomKK->k_omega.sync<LMPDeviceType>();
    if (mask & TORQUE_MASK) atomKK->k_torque.sync<LMPDeviceType>();
  } else {
    if (mask & X_MASK) atomKK->k_x.sync<LMPHostType>();
    if (mask & V_MASK) atomKK->k_v.sync<LMPHostType>();
    if (mask & F_MASK) atomKK->k_f.sync<LMPHostType>();
    if (mask & TAG_MASK) atomKK->k_tag.sync<LMPHostType>();
    if (mask & TYPE_MASK) atomKK->k_type.sync<LMPHostType>();
    if (mask & MASK_MASK) atomKK->k_mask.sync<LMPHostType>();
    if (mask & IMAGE_MASK) atomKK->k_image.sync<LMPHostType>();
    if (mask & RADIUS_MASK) atomKK->k_radius.sync<LMPHostType>();
    if (mask & RMASS_MASK) atomKK->k_rmass.sync<LMPHostType>();
    if (mask & OMEGA_MASK) atomKK->k_omega.sync<LMPHostType>();
    if (mask & TORQUE_MASK) atomKK->k_torque.sync<LMPHostType>();
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecSphereKokkos::sync_overlapping_device(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if ((mask & X_MASK) && atomKK->k_x.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_x_array>(atomKK->k_x,space);
    if ((mask & V_MASK) && atomKK->k_v.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_v_array>(atomKK->k_v,space);
    if ((mask & F_MASK) && atomKK->k_f.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_f_array>(atomKK->k_f,space);
    if ((mask & TAG_MASK) && atomKK->k_tag.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_tagint_1d>(atomKK->k_tag,space);
    if ((mask & TYPE_MASK) && atomKK->k_type.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_int_1d>(atomKK->k_type,space);
    if ((mask & MASK_MASK) && atomKK->k_mask.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_int_1d>(atomKK->k_mask,space);
    if ((mask & IMAGE_MASK) && atomKK->k_image.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_imageint_1d>(atomKK->k_image,space);
    if ((mask & RADIUS_MASK) && atomKK->k_radius.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_float_1d>(atomKK->k_radius,space);
    if ((mask & RMASS_MASK) && atomKK->k_rmass.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_float_1d>(atomKK->k_rmass,space);
    if ((mask & OMEGA_MASK) && atomKK->k_omega.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_v_array>(atomKK->k_omega,space);
    if ((mask & TORQUE_MASK) && atomKK->k_torque.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_f_array>(atomKK->k_torque,space);
  } else {
    if ((mask & X_MASK) && atomKK->k_x.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_x_array>(atomKK->k_x,space);
    if ((mask & V_MASK) && atomKK->k_v.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_v_array>(atomKK->k_v,space);
    if ((mask & F_MASK) && atomKK->k_f.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_f_array>(atomKK->k_f,space);
    if ((mask & TAG_MASK) && atomKK->k_tag.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_tagint_1d>(atomKK->k_tag,space);
    if ((mask & TYPE_MASK) && atomKK->k_type.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_int_1d>(atomKK->k_type,space);
    if ((mask & MASK_MASK) && atomKK->k_mask.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_int_1d>(atomKK->k_mask,space);
    if ((mask & IMAGE_MASK) && atomKK->k_image.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_imageint_1d>(atomKK->k_image,space);
    if ((mask & RADIUS_MASK) && atomKK->k_radius.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_float_1d>(atomKK->k_radius,space);
    if ((mask & RMASS_MASK) && atomKK->k_rmass.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_float_1d>(atomKK->k_rmass,space);
    if ((mask & OMEGA_MASK) && atomKK->k_omega.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_v_array>(atomKK->k_omega,space);
    if ((mask & TORQUE_MASK) && atomKK->k_torque.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_f_array>(atomKK->k_torque,space);
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecSphereKokkos::modified(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.modify<LMPDeviceType>();
    if (mask & V_MASK) atomKK->k_v.modify<LMPDeviceType>();
    if (mask & F_MASK) atomKK->k_f.modify<LMPDeviceType>();
    if (mask & TAG_MASK) atomKK->k_tag.modify<LMPDeviceType>();
    if (mask & TYPE_MASK) atomKK->k_type.modify<LMPDeviceType>();
    if (mask & MASK_MASK) atomKK->k_mask.modify<LMPDeviceType>();
    if (mask & IMAGE_MASK) atomKK->k_image.modify<LMPDeviceType>();
    if (mask & RADIUS_MASK) atomKK->k_radius.modify<LMPDeviceType>();
    if (mask & RMASS_MASK) atomKK->k_rmass.modify<LMPDeviceType>();
    if (mask & OMEGA_MASK) atomKK->k_omega.modify<LMPDeviceType>();
    if (mask & TORQUE_MASK) atomKK->k_torque.modify<LMPDeviceType>();
  } else {
    if (mask & X_MASK) atomKK->k_x.modify<LMPHostType>();
    if (mask & V_MASK) atomKK->k_v.modify<LMPHostType>();
    if (mask & F_MASK) atomKK->k_f.modify<LMPHostType>();
    if (mask & TAG_MASK) atomKK->k_tag.modify<LMPHostType>();
    if (mask & TYPE_MASK) atomKK->k_type.modify<LMPHostType>();
    if (mask & MASK_MASK) atomKK->k_mask.modify<LMPHostType>();
    if (mask & IMAGE_MASK) atomKK->k_image.modify<LMPHostType>();
    if (mask & RADIUS_MASK) atomKK->k_radius.modify<LMPHostType>();
    if (mask & RMASS_MASK) atomKK->k_rmass.modify<LMPHostType>();
    if (mask & OMEGA_MASK) atomKK->k_omega.modify<LMPHostType>();
    if (mask & TORQUE_MASK) atomKK->k_torque.modify<LMPHostType>();
  }
}
