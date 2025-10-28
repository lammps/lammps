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

#include "atom_vec_sphere_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm_kokkos.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "modify.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

AtomVecSphereKokkos::AtomVecSphereKokkos(LAMMPS *lmp) : AtomVec(lmp),
AtomVecKokkos(lmp), AtomVecSphere(lmp)
{
  no_border_vel_flag = 0;
  unpack_exchange_indices_flag = 1;
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecSphereKokkos::grow(int n)
{
  auto DELTA = LMP_KOKKOS_AV_DELTA;
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

  memoryKK->grow_kokkos(atomKK->k_x,atomKK->x,nmax,"atom:x");
  memoryKK->grow_kokkos(atomKK->k_v,atomKK->v,nmax,"atom:v");
  memoryKK->grow_kokkos(atomKK->k_f,atomKK->f,nmax,"atom:f");
  memoryKK->grow_kokkos(atomKK->k_radius,atomKK->radius,nmax,"atom:radius");
  memoryKK->grow_kokkos(atomKK->k_rmass,atomKK->rmass,nmax,"atom:rmass");
  memoryKK->grow_kokkos(atomKK->k_omega,atomKK->omega,nmax,"atom:omega");
  memoryKK->grow_kokkos(atomKK->k_torque,atomKK->torque,nmax,"atom:torque");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);

  grow_pointers();
  atomKK->sync(Host,ALL_MASK);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecSphereKokkos::grow_pointers()
{
  tag = atomKK->tag;
  d_tag = atomKK->k_tag.view_device();
  h_tag = atomKK->k_tag.view_host();

  type = atomKK->type;
  d_type = atomKK->k_type.view_device();
  h_type = atomKK->k_type.view_host();
  mask = atomKK->mask;
  d_mask = atomKK->k_mask.view_device();
  h_mask = atomKK->k_mask.view_host();
  image = atomKK->image;
  d_image = atomKK->k_image.view_device();
  h_image = atomKK->k_image.view_host();

  x = atomKK->x;
  d_x = atomKK->k_x.view_device();
  h_x = atomKK->k_x.view_hostkk();
  v = atomKK->v;
  d_v = atomKK->k_v.view_device();
  h_v = atomKK->k_v.view_hostkk();
  f = atomKK->f;
  d_f = atomKK->k_f.view_device();
  h_f = atomKK->k_f.view_hostkk();
  radius = atomKK->radius;
  d_radius = atomKK->k_radius.view_device();
  h_radius = atomKK->k_radius.view_hostkk();
  rmass = atomKK->rmass;
  d_rmass = atomKK->k_rmass.view_device();
  h_rmass = atomKK->k_rmass.view_hostkk();
  omega = atomKK->omega;
  d_omega = atomKK->k_omega.view_device();
  h_omega = atomKK->k_omega.view_hostkk();
  torque = atomKK->torque;
  d_torque = atomKK->k_torque.view_device();
  h_torque = atomKK->k_torque.view_hostkk();
}

/* ----------------------------------------------------------------------
   sort atom arrays on device
------------------------------------------------------------------------- */

void AtomVecSphereKokkos::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  atomKK->sync(Device, ALL_MASK & ~F_MASK & ~TORQUE_MASK);

  Sorter.sort(LMPDeviceType(), d_tag);
  Sorter.sort(LMPDeviceType(), d_type);
  Sorter.sort(LMPDeviceType(), d_mask);
  Sorter.sort(LMPDeviceType(), d_image);
  Sorter.sort(LMPDeviceType(), d_x);
  Sorter.sort(LMPDeviceType(), d_v);
  Sorter.sort(LMPDeviceType(), d_radius);
  Sorter.sort(LMPDeviceType(), d_rmass);
  Sorter.sort(LMPDeviceType(), d_omega);

  atomKK->modified(Device, ALL_MASK & ~F_MASK & ~TORQUE_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG,int TRICLINIC>
struct AtomVecSphereKokkos_PackComm {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr_randomread _x;
  typename AT::t_kkfloat_1d _radius,_rmass;
  typename AT::t_double_2d_lr_um _buf;
  typename AT::t_int_1d_const _list;
  double _xprd,_yprd,_zprd,_xy,_xz,_yz;
  double _pbc[6];

  AtomVecSphereKokkos_PackComm(
    const typename DAT::ttransform_kkfloat_1d_3_lr &x,
    const typename DAT::ttransform_kkfloat_1d &radius,
    const typename DAT::ttransform_kkfloat_1d &rmass,
    const typename DAT::tdual_double_2d_lr &buf,
    const typename DAT::tdual_int_1d &list,
    const double &xprd, const double &yprd, const double &zprd,
    const double &xy, const double &xz, const double &yz, const int* const pbc):
    _x(x.view<DeviceType>()),
    _radius(radius.view<DeviceType>()),
    _rmass(rmass.view<DeviceType>()),
    _list(list.view<DeviceType>()),
    _xprd(xprd),_yprd(yprd),_zprd(zprd),
    _xy(xy),_xz(xz),_yz(yz) {
    const size_t elements = 5;
    const size_t maxsend = (buf.view<DeviceType>().extent(0)*buf.view<DeviceType>().extent(1))/elements;
    _buf = typename AT::t_double_2d_lr_um(buf.view<DeviceType>().data(),maxsend,elements);
    _pbc[0] = pbc[0]; _pbc[1] = pbc[1]; _pbc[2] = pbc[2];
    _pbc[3] = pbc[3]; _pbc[4] = pbc[4]; _pbc[5] = pbc[5];
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(i);
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
  const DAT::tdual_int_1d &list,
  const DAT::tdual_double_2d_lr &buf,
  const int &pbc_flag,
  const int* const pbc)
{
  // Fallback to AtomVecKokkos if radvary == 0
  if (radvary == 0)
    return AtomVecKokkos::pack_comm_kokkos(n,list,buf,pbc_flag,pbc);
  // Check whether to always run forward communication on the host
  // Choose correct forward PackComm kernel
  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->sync(HostKK,X_MASK|RADIUS_MASK|RMASS_MASK);
    if (pbc_flag) {
      if (domain->triclinic) {
        struct AtomVecSphereKokkos_PackComm<LMPHostType,1,1> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          buf,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecSphereKokkos_PackComm<LMPHostType,1,0> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          buf,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      }
    } else {
      if (domain->triclinic) {
        struct AtomVecSphereKokkos_PackComm<LMPHostType,0,1> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          buf,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecSphereKokkos_PackComm<LMPHostType,0,0> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          buf,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      }
    }
  } else {
    atomKK->sync(Device,X_MASK|RADIUS_MASK|RMASS_MASK);
    if (pbc_flag) {
      if (domain->triclinic) {
        struct AtomVecSphereKokkos_PackComm<LMPDeviceType,1,1> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          buf,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecSphereKokkos_PackComm<LMPDeviceType,1,0> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          buf,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      }
    } else {
      if (domain->triclinic) {
        struct AtomVecSphereKokkos_PackComm<LMPDeviceType,0,1> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          buf,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecSphereKokkos_PackComm<LMPDeviceType,0,0> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          buf,list,
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
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr_randomread _x;
  typename AT::t_int_1d _mask;
  typename AT::t_kkfloat_1d _radius,_rmass;
  typename AT::t_kkfloat_1d_3 _v, _omega;
  typename AT::t_double_2d_lr_um _buf;
  typename AT::t_int_1d_const _list;
  double _xprd,_yprd,_zprd,_xy,_xz,_yz;
  double _pbc[6];
  double _h_rate[6];
  const int _deform_vremap;

  AtomVecSphereKokkos_PackCommVel(
    const typename DAT::ttransform_kkfloat_1d_3_lr &x,
    const typename DAT::tdual_int_1d &mask,
    const typename DAT::ttransform_kkfloat_1d &radius,
    const typename DAT::ttransform_kkfloat_1d &rmass,
    const typename DAT::ttransform_kkfloat_1d_3 &v,
    const typename DAT::ttransform_kkfloat_1d_3 &omega,
    const typename DAT::tdual_double_2d_lr &buf,
    const typename DAT::tdual_int_1d &list,
    const double &xprd, const double &yprd, const double &zprd,
    const double &xy, const double &xz, const double &yz, const int* const pbc,
    const double * const h_rate,
    const int &deform_vremap):
    _x(x.view<DeviceType>()),
    _mask(mask.view<DeviceType>()),
    _radius(radius.view<DeviceType>()),
    _rmass(rmass.view<DeviceType>()),
    _v(v.view<DeviceType>()),
    _omega(omega.view<DeviceType>()),
    _list(list.view<DeviceType>()),
    _xprd(xprd),_yprd(yprd),_zprd(zprd),
    _xy(xy),_xz(xz),_yz(yz),
    _deform_vremap(deform_vremap)
  {
    const size_t elements = 9 + 2 * RADVARY;
    const int maxsend = (buf.template view<DeviceType>().extent(0)*buf.template view<DeviceType>().extent(1))/elements;
    _buf = typename AT::t_double_2d_lr_um(buf.view<DeviceType>().data(),maxsend,elements);
    _pbc[0] = pbc[0]; _pbc[1] = pbc[1]; _pbc[2] = pbc[2];
    _pbc[3] = pbc[3]; _pbc[4] = pbc[4]; _pbc[5] = pbc[5];
    _h_rate[0] = h_rate[0]; _h_rate[1] = h_rate[1]; _h_rate[2] = h_rate[2];
    _h_rate[3] = h_rate[3]; _h_rate[4] = h_rate[4]; _h_rate[5] = h_rate[5];
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(i);
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
  const DAT::tdual_int_1d &list,
  const DAT::tdual_double_2d_lr &buf,
  const int &pbc_flag,
  const int* const pbc)
{
  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->sync(HostKK,X_MASK|RADIUS_MASK|RMASS_MASK|V_MASK|OMEGA_MASK);
    if (pbc_flag) {
      if (deform_vremap) {
        if (domain->triclinic) {
          if (radvary == 0) {
            struct AtomVecSphereKokkos_PackCommVel<LMPHostType,0,1,1,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecSphereKokkos_PackCommVel<LMPHostType,1,1,1,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,
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
              buf,list,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecSphereKokkos_PackCommVel<LMPHostType,1,1,0,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          }
        }
      } else {
        if (domain->triclinic) {
          if (radvary == 0) {
            struct AtomVecSphereKokkos_PackCommVel<LMPHostType,0,1,1,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecSphereKokkos_PackCommVel<LMPHostType,1,1,1,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,
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
              buf,list,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecSphereKokkos_PackCommVel<LMPHostType,1,1,0,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          }
        }
      }
    } else {
      if (domain->triclinic) {
        if (radvary == 0) {
          struct AtomVecSphereKokkos_PackCommVel<LMPHostType,0,0,1,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecSphereKokkos_PackCommVel<LMPHostType,1,0,1,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,
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
            buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecSphereKokkos_PackCommVel<LMPHostType,1,0,0,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        }
      }
    }
  } else {
    atomKK->sync(Device,X_MASK|RADIUS_MASK|RMASS_MASK|V_MASK|OMEGA_MASK);
    if (pbc_flag) {
      if (deform_vremap) {
        if (domain->triclinic) {
          if (radvary == 0) {
            struct AtomVecSphereKokkos_PackCommVel<LMPDeviceType,0,1,1,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecSphereKokkos_PackCommVel<LMPDeviceType,1,1,1,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,
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
              buf,list,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecSphereKokkos_PackCommVel<LMPDeviceType,1,1,0,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          }
        }
      } else {
        if (domain->triclinic) {
          if (radvary == 0) {
            struct AtomVecSphereKokkos_PackCommVel<LMPDeviceType,0,1,1,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecSphereKokkos_PackCommVel<LMPDeviceType,1,1,1,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,
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
              buf,list,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecSphereKokkos_PackCommVel<LMPDeviceType,1,1,0,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          }
        }
      }
    } else {
      if (domain->triclinic) {
        if (radvary == 0) {
          struct AtomVecSphereKokkos_PackCommVel<LMPDeviceType,0,0,1,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecSphereKokkos_PackCommVel<LMPDeviceType,1,0,1,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,
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
            buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecSphereKokkos_PackCommVel<LMPDeviceType,1,0,0,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,
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
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr_randomread _x;
  typename AT::t_kkfloat_1d_3_lr _xw;
  typename AT::t_kkfloat_1d _radius,_rmass;
  int _nfirst;
  typename AT::t_int_1d_const _list;
  double _xprd,_yprd,_zprd,_xy,_xz,_yz;
  double _pbc[6];

  AtomVecSphereKokkos_PackCommSelf(
    const typename DAT::ttransform_kkfloat_1d_3_lr &x,
    const typename DAT::ttransform_kkfloat_1d &radius,
    const typename DAT::ttransform_kkfloat_1d &rmass,
    const int &nfirst,
    const typename DAT::tdual_int_1d &list,
    const double &xprd, const double &yprd, const double &zprd,
    const double &xy, const double &xz, const double &yz, const int* const pbc):
    _x(x.view<DeviceType>()),_xw(x.view<DeviceType>()),
    _radius(radius.view<DeviceType>()),
    _rmass(rmass.view<DeviceType>()),
    _nfirst(nfirst),_list(list.view<DeviceType>()),
    _xprd(xprd),_yprd(yprd),_zprd(zprd),
    _xy(xy),_xz(xz),_yz(yz) {
    _pbc[0] = pbc[0]; _pbc[1] = pbc[1]; _pbc[2] = pbc[2];
    _pbc[3] = pbc[3]; _pbc[4] = pbc[4]; _pbc[5] = pbc[5];
  };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(i);
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
  const int &n, const DAT::tdual_int_1d &list,
  const int nfirst, const int &pbc_flag, const int* const pbc) {
  // Fallback to AtomVecKokkos if radvary == 0
  if (radvary == 0)
    return AtomVecKokkos::pack_comm_self(n,list,nfirst,pbc_flag,pbc);
  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->sync(HostKK,X_MASK|RADIUS_MASK|RMASS_MASK);
    atomKK->modified(HostKK,X_MASK|RADIUS_MASK|RMASS_MASK);
    if (pbc_flag) {
      if (domain->triclinic) {
        struct AtomVecSphereKokkos_PackCommSelf<LMPHostType,1,1> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          nfirst,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecSphereKokkos_PackCommSelf<LMPHostType,1,0> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          nfirst,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      }
    } else {
      if (domain->triclinic) {
        struct AtomVecSphereKokkos_PackCommSelf<LMPHostType,0,1> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          nfirst,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecSphereKokkos_PackCommSelf<LMPHostType,0,0> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          nfirst,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      }
    }
  } else {
    atomKK->sync(Device,X_MASK|RADIUS_MASK|RMASS_MASK);
    atomKK->modified(Device,X_MASK|RADIUS_MASK|RMASS_MASK);
    if (pbc_flag) {
      if (domain->triclinic) {
        struct AtomVecSphereKokkos_PackCommSelf<LMPDeviceType,1,1> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          nfirst,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecSphereKokkos_PackCommSelf<LMPDeviceType,1,0> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          nfirst,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      }
    } else {
      if (domain->triclinic) {
        struct AtomVecSphereKokkos_PackCommSelf<LMPDeviceType,0,1> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          nfirst,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecSphereKokkos_PackCommSelf<LMPDeviceType,0,0> f(
          atomKK->k_x,
          atomKK->k_radius,atomKK->k_rmass,
          nfirst,list,
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
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr _x;
  typename AT::t_kkfloat_1d _radius,_rmass;
  typename AT::t_double_2d_lr_const_um _buf;
  int _first;

  AtomVecSphereKokkos_UnpackComm(
    const typename DAT::ttransform_kkfloat_1d_3_lr &x,
    const typename DAT::ttransform_kkfloat_1d &radius,
    const typename DAT::ttransform_kkfloat_1d &rmass,
    const typename DAT::tdual_double_2d_lr &buf,
    const int& first):
    _x(x.view<DeviceType>()),
    _radius(radius.view<DeviceType>()),
    _rmass(rmass.view<DeviceType>()),
    _first(first)
  {
    const size_t elements = 5;
    const size_t maxsend = (buf.view<DeviceType>().extent(0)*buf.view<DeviceType>().extent(1))/elements;
    _buf = typename AT::t_double_2d_lr_const_um(buf.view<DeviceType>().data(),maxsend,elements);
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
  const DAT::tdual_double_2d_lr &buf) {
  // Fallback to AtomVecKokkos if radvary == 0
  if (radvary == 0) {
    AtomVecKokkos::unpack_comm_kokkos(n,first,buf);
    return;
  }
  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->modified(HostKK,X_MASK|RADIUS_MASK|RMASS_MASK);
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
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr _x;
  typename AT::t_kkfloat_1d _radius,_rmass;
  typename AT::t_kkfloat_1d_3 _v, _omega;
  typename AT::t_double_2d_lr_const _buf;
  int _first;

  AtomVecSphereKokkos_UnpackCommVel(
    const typename DAT::ttransform_kkfloat_1d_3_lr &x,
    const typename DAT::ttransform_kkfloat_1d &radius,
    const typename DAT::ttransform_kkfloat_1d &rmass,
    const typename DAT::ttransform_kkfloat_1d_3 &v,
    const typename DAT::ttransform_kkfloat_1d_3 &omega,
    const typename DAT::tdual_double_2d_lr &buf,
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
  const DAT::tdual_double_2d_lr &buf) {
  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->modified(HostKK,X_MASK|RADIUS_MASK|RMASS_MASK|V_MASK|OMEGA_MASK);
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

template<class DeviceType,int PBC_FLAG>
struct AtomVecSphereKokkos_PackBorder {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_double_2d_lr_um _buf;
  const typename AT::t_int_1d_const _list;
  const typename AT::t_kkfloat_1d_3_lr_randomread _x;
  const typename AT::t_tagint_1d _tag;
  const typename AT::t_int_1d _type;
  const typename AT::t_int_1d _mask;
  typename AT::t_kkfloat_1d _radius,_rmass;
  double _dx,_dy,_dz;

  AtomVecSphereKokkos_PackBorder(
    const typename AT::t_double_2d_lr &buf,
    const typename AT::t_int_1d_const &list,
    const typename AT::t_kkfloat_1d_3_lr &x,
    const typename AT::t_tagint_1d &tag,
    const typename AT::t_int_1d &type,
    const typename AT::t_int_1d &mask,
    const typename AT::t_kkfloat_1d &radius,
    const typename AT::t_kkfloat_1d &rmass,
    const double &dx, const double &dy, const double &dz):
    _buf(buf),_list(list),
    _x(x),_tag(tag),_type(type),_mask(mask),
    _radius(radius),
    _rmass(rmass),
    _dx(dx),_dy(dy),_dz(dz)
  {
    const size_t elements = 8;
    const int maxsend = (buf.extent(0)*buf.extent(1))/elements;
    _buf = typename AT::t_double_2d_lr_um(buf.data(),maxsend,elements);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(i);
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
  int n, DAT::tdual_int_1d k_sendlist, DAT::tdual_double_2d_lr buf,
  int pbc_flag, int *pbc, ExecutionSpace space)
{
  double dx,dy,dz;

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
    if (space==Host) {
      AtomVecSphereKokkos_PackBorder<LMPHostType,1> f(
        buf.view_host(), k_sendlist.view_host(),
        h_x,h_tag,h_type,h_mask,
        h_radius,h_rmass,
        dx,dy,dz);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecSphereKokkos_PackBorder<LMPDeviceType,1> f(
        buf.view_device(), k_sendlist.view_device(),
        d_x,d_tag,d_type,d_mask,
        d_radius,d_rmass,
        dx,dy,dz);
      Kokkos::parallel_for(n,f);
    }
  } else {
    dx = dy = dz = 0;
    if (space==Host) {
      AtomVecSphereKokkos_PackBorder<LMPHostType,0> f(
        buf.view_host(), k_sendlist.view_host(),
        h_x,h_tag,h_type,h_mask,
        h_radius,h_rmass,
        dx,dy,dz);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecSphereKokkos_PackBorder<LMPDeviceType,0> f(
        buf.view_device(), k_sendlist.view_device(),
        d_x,d_tag,d_type,d_mask,
        d_radius,d_rmass,
        dx,dy,dz);
      Kokkos::parallel_for(n,f);
    }
  }
  return n*size_border;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG,int DEFORM_VREMAP>
struct AtomVecSphereKokkos_PackBorderVel {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_double_2d_lr_um _buf;
  const typename AT::t_int_1d_const _list;
  const typename AT::t_kkfloat_1d_3_lr_randomread _x;
  const typename AT::t_tagint_1d _tag;
  const typename AT::t_int_1d _type;
  const typename AT::t_int_1d _mask;
  typename AT::t_kkfloat_1d _radius,_rmass;
  typename AT::t_kkfloat_1d_3 _v, _omega;
  double _dx,_dy,_dz, _dvx, _dvy, _dvz;
  const int _deform_groupbit;

  AtomVecSphereKokkos_PackBorderVel(
    const typename AT::t_double_2d_lr &buf,
    const typename AT::t_int_1d_const &list,
    const typename AT::t_kkfloat_1d_3_lr &x,
    const typename AT::t_tagint_1d &tag,
    const typename AT::t_int_1d &type,
    const typename AT::t_int_1d &mask,
    const typename AT::t_kkfloat_1d &radius,
    const typename AT::t_kkfloat_1d &rmass,
    const typename AT::t_kkfloat_1d_3 &v,
    const typename AT::t_kkfloat_1d_3 &omega,
    const double &dx, const double &dy, const double &dz,
    const double &dvx, const double &dvy, const double &dvz,
    const int &deform_groupbit):
    _buf(buf),_list(list),
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
    _buf = typename AT::t_double_2d_lr_um(buf.data(),maxsend,elements);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(i);
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
  int n, DAT::tdual_int_1d k_sendlist, DAT::tdual_double_2d_lr buf,
  int pbc_flag, int *pbc, ExecutionSpace space)
{
  double dx=0,dy=0,dz=0;
  double dvx=0,dvy=0,dvz=0;

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
      if (space==Host) {
        AtomVecSphereKokkos_PackBorderVel<LMPHostType,1,0> f(
          buf.view_host(), k_sendlist.view_host(),
          h_x,h_tag,h_type,h_mask,
          h_radius,h_rmass,
          h_v, h_omega,
          dx,dy,dz,dvx,dvy,dvz,
          deform_groupbit);
        Kokkos::parallel_for(n,f);
      } else {
        AtomVecSphereKokkos_PackBorderVel<LMPDeviceType,1,0> f(
          buf.view_device(), k_sendlist.view_device(),
          d_x,d_tag,d_type,d_mask,
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
      if (space==Host) {
        AtomVecSphereKokkos_PackBorderVel<LMPHostType,1,1> f(
          buf.view_host(), k_sendlist.view_host(),
          h_x,h_tag,h_type,h_mask,
          h_radius,h_rmass,
          h_v, h_omega,
          dx,dy,dz,dvx,dvy,dvz,
          deform_groupbit);
        Kokkos::parallel_for(n,f);
      } else {
        AtomVecSphereKokkos_PackBorderVel<LMPDeviceType,1,1> f(
          buf.view_device(), k_sendlist.view_device(),
          d_x,d_tag,d_type,d_mask,
          d_radius,d_rmass,
          d_v, d_omega,
          dx,dy,dz,dvx,dvy,dvz,
          deform_groupbit);
        Kokkos::parallel_for(n,f);
      }
    }
  } else {
    if (space==Host) {
      AtomVecSphereKokkos_PackBorderVel<LMPHostType,0,0> f(
        buf.view_host(), k_sendlist.view_host(),
        h_x,h_tag,h_type,h_mask,
        h_radius,h_rmass,
        h_v, h_omega,
        dx,dy,dz,dvx,dvy,dvz,
        deform_groupbit);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecSphereKokkos_PackBorderVel<LMPDeviceType,0,0> f(
        buf.view_device(), k_sendlist.view_device(),
        d_x,d_tag,d_type,d_mask,
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

template<class DeviceType>
struct AtomVecSphereKokkos_UnpackBorder {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_double_2d_lr_const_um _buf;
  typename AT::t_kkfloat_1d_3_lr _x;
  typename AT::t_tagint_1d _tag;
  typename AT::t_int_1d _type;
  typename AT::t_int_1d _mask;
  typename AT::t_kkfloat_1d _radius,_rmass;
  int _first;

  AtomVecSphereKokkos_UnpackBorder(
    const typename AT::t_double_2d_lr &buf,
    const typename AT::t_kkfloat_1d_3_lr &x,
    const typename AT::t_tagint_1d &tag,
    const typename AT::t_int_1d &type,
    const typename AT::t_int_1d &mask,
    const typename AT::t_kkfloat_1d &radius,
    const typename AT::t_kkfloat_1d &rmass,
    const int& first):
    _buf(buf),_x(x),_tag(tag),_type(type),_mask(mask),
    _radius(radius),
    _rmass(rmass),
    _first(first)
  {
    const size_t elements = 8;
    const int maxsend = (buf.extent(0)*buf.extent(1))/elements;
    _buf = typename AT::t_double_2d_lr_const_um(buf.data(),maxsend,elements);
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
                                               const DAT::tdual_double_2d_lr &buf,ExecutionSpace space) {
  while (first+n >= nmax) grow(0);
  if (space==Host) {
    struct AtomVecSphereKokkos_UnpackBorder<LMPHostType> f(buf.view_host(),
      h_x,h_tag,h_type,h_mask,
      h_radius,h_rmass,
      first);
    Kokkos::parallel_for(n,f);
  } else {
    struct AtomVecSphereKokkos_UnpackBorder<LMPDeviceType> f(buf.view_device(),
      d_x,d_tag,d_type,d_mask,
      d_radius,d_rmass,
      first);
    Kokkos::parallel_for(n,f);
  }

  atomKK->modified(space,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|
                 RADIUS_MASK|RMASS_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecSphereKokkos_UnpackBorderVel {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_double_2d_lr_const_um _buf;
  typename AT::t_kkfloat_1d_3_lr _x;
  typename AT::t_tagint_1d _tag;
  typename AT::t_int_1d _type;
  typename AT::t_int_1d _mask;
  typename AT::t_kkfloat_1d _radius,_rmass;
  typename AT::t_kkfloat_1d_3 _v;
  typename AT::t_kkfloat_1d_3 _omega;
  int _first;

  AtomVecSphereKokkos_UnpackBorderVel(
    const typename AT::t_double_2d_lr_const &buf,
    const typename AT::t_kkfloat_1d_3_lr &x,
    const typename AT::t_tagint_1d &tag,
    const typename AT::t_int_1d &type,
    const typename AT::t_int_1d &mask,
    const typename AT::t_kkfloat_1d &radius,
    const typename AT::t_kkfloat_1d &rmass,
    const typename AT::t_kkfloat_1d_3 &v,
    const typename AT::t_kkfloat_1d_3 &omega,
    const int& first):
    _buf(buf),_x(x),_tag(tag),_type(type),_mask(mask),
    _radius(radius),
    _rmass(rmass),
    _v(v), _omega(omega),
    _first(first)
  {
    const size_t elements = 14;
    const int maxsend = (buf.extent(0)*buf.extent(1))/elements;
    _buf = typename AT::t_double_2d_lr_const_um(buf.data(),maxsend,elements);
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
  const DAT::tdual_double_2d_lr &buf,ExecutionSpace space) {
  while (first+n >= nmax) grow(0);
  if (space==Host) {
    struct AtomVecSphereKokkos_UnpackBorderVel<LMPHostType> f(buf.view_host(),
      h_x,h_tag,h_type,h_mask,
      h_radius,h_rmass,
      h_v, h_omega,
      first);
    Kokkos::parallel_for(n,f);
  } else {
    struct AtomVecSphereKokkos_UnpackBorderVel<LMPDeviceType> f(buf.view_device(),
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

template<class DeviceType>
struct AtomVecSphereKokkos_PackExchangeFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_kkfloat_1d_3_lr_randomread _x;
  typename AT::t_kkfloat_1d_3_randomread _v;
  typename AT::t_tagint_1d_randomread _tag;
  typename AT::t_int_1d_randomread _type;
  typename AT::t_int_1d_randomread _mask;
  typename AT::t_imageint_1d_randomread _image;
  typename AT::t_kkfloat_1d_randomread _radius,_rmass;
  typename AT::t_kkfloat_1d_3_randomread _omega;
  typename AT::t_kkfloat_1d_3_lr _xw;
  typename AT::t_kkfloat_1d_3 _vw;
  typename AT::t_tagint_1d _tagw;
  typename AT::t_int_1d _typew;
  typename AT::t_int_1d _maskw;
  typename AT::t_imageint_1d _imagew;
  typename AT::t_kkfloat_1d _radiusw,_rmassw;
  typename AT::t_kkfloat_1d_3 _omegaw;
  typename AT::t_double_2d_lr_um _buf;
  typename AT::t_int_1d_const _sendlist;
  typename AT::t_int_1d_const _copylist;
  int _size_exchange;

  AtomVecSphereKokkos_PackExchangeFunctor(
    const AtomKokkos* atom,
    const DAT::tdual_double_2d_lr buf,
    DAT::tdual_int_1d sendlist,
    DAT::tdual_int_1d copylist):
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
    _size_exchange(atom->avecKK->size_exchange) {
    const int maxsend = (buf.template view<DeviceType>().extent(0)*buf.template view<DeviceType>().extent(1))/_size_exchange;

    _buf = typename AT::t_double_2d_lr_um(buf.template view<DeviceType>().data(),maxsend,_size_exchange);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &mysend) const {
    const int i = _sendlist(mysend);
    _buf(mysend,0) = _size_exchange;
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
  DAT::tdual_double_2d_lr &k_buf,
  DAT::tdual_int_1d k_sendlist,
  DAT::tdual_int_1d k_copylist,
  ExecutionSpace space)
{
  size_exchange = 16;

  if (nsend > (int) (k_buf.view_host().extent(0)*k_buf.view_host().extent(1))/size_exchange) {
    int newsize = nsend*17/k_buf.view_host().extent(1)+1;
    k_buf.resize(newsize,k_buf.view_host().extent(1));
  }
  atomKK->sync(space,X_MASK | V_MASK | TAG_MASK | TYPE_MASK |
             MASK_MASK | IMAGE_MASK| RADIUS_MASK | RMASS_MASK |
             OMEGA_MASK);

  if (space == HostKK) {
    AtomVecSphereKokkos_PackExchangeFunctor<LMPHostType> f(atomKK,k_buf,k_sendlist,k_copylist);
    Kokkos::parallel_for(nsend,f);
  } else {
    AtomVecSphereKokkos_PackExchangeFunctor<LMPDeviceType> f(atomKK,k_buf,k_sendlist,k_copylist);
    Kokkos::parallel_for(nsend,f);
  }
  return nsend*size_exchange;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int OUTPUT_INDICES>
struct AtomVecSphereKokkos_UnpackExchangeFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_kkfloat_1d_3_lr _x;
  typename AT::t_kkfloat_1d_3 _v;
  typename AT::t_tagint_1d _tag;
  typename AT::t_int_1d _type;
  typename AT::t_int_1d _mask;
  typename AT::t_imageint_1d _image;
  typename AT::t_kkfloat_1d _radius;
  typename AT::t_kkfloat_1d _rmass;
  typename AT::t_kkfloat_1d_3 _omega;
  typename AT::t_double_2d_lr_um _buf;
  typename AT::t_int_1d _nlocal;
  typename AT::t_int_1d _indices;
  int _dim;
  double _lo,_hi;
  int _size_exchange;

  AtomVecSphereKokkos_UnpackExchangeFunctor(
    const AtomKokkos* atom,
    const DAT::tdual_double_2d_lr buf,
    DAT::tdual_int_1d nlocal,
    DAT::tdual_int_1d indices,
    int dim, double lo, double hi):
      _x(atom->k_x.view<DeviceType>()),
      _v(atom->k_v.view<DeviceType>()),
      _tag(atom->k_tag.view<DeviceType>()),
      _type(atom->k_type.view<DeviceType>()),
      _mask(atom->k_mask.view<DeviceType>()),
      _image(atom->k_image.view<DeviceType>()),
      _radius(atom->k_radius.view<DeviceType>()),
      _rmass(atom->k_rmass.view<DeviceType>()),
      _omega(atom->k_omega.view<DeviceType>()),
      _nlocal(nlocal.template view<DeviceType>()),
      _indices(indices.template view<DeviceType>()),
      _dim(dim),_lo(lo),_hi(hi),_size_exchange(atom->avecKK->size_exchange) {
    const size_t size_exchange = 16;
    const int maxsendlist = (buf.template view<DeviceType>().extent(0)*buf.template view<DeviceType>().extent(1))/size_exchange;

    buffer_view<DeviceType>(_buf,buf,maxsendlist,size_exchange);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &myrecv) const {
    double x = _buf(myrecv,_dim+1);
    int i = -1;
    if (x >= _lo && x < _hi) {
      i = Kokkos::atomic_fetch_add(&_nlocal(0),1);
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
    if (OUTPUT_INDICES)
      _indices(myrecv) = i;
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecSphereKokkos::unpack_exchange_kokkos(DAT::tdual_double_2d_lr &k_buf, int nrecv, int nlocal,
                                                int dim, double lo, double hi, ExecutionSpace space,
                                                DAT::tdual_int_1d &k_indices)
{
  while (nlocal + nrecv/size_exchange >= nmax) grow(0);

  if (space == HostKK) {
    k_count.view_host()(0) = nlocal;
    if (k_indices.view_host().data()) {
      AtomVecSphereKokkos_UnpackExchangeFunctor<LMPHostType,1> f(atomKK,k_buf,k_count,k_indices,dim,lo,hi);
      Kokkos::parallel_for(nrecv/size_exchange,f);
    } else {
      AtomVecSphereKokkos_UnpackExchangeFunctor<LMPHostType,0> f(atomKK,k_buf,k_count,k_indices,dim,lo,hi);
      Kokkos::parallel_for(nrecv/size_exchange,f);
    }
  } else {
    k_count.view_host()(0) = nlocal;
    k_count.modify_host();
    k_count.sync_device();
    if (k_indices.view_host().data()) {
      AtomVecSphereKokkos_UnpackExchangeFunctor<LMPDeviceType,1> f(atomKK,k_buf,k_count,k_indices,dim,lo,hi);
      Kokkos::parallel_for(nrecv/size_exchange,f);
    } else {
      AtomVecSphereKokkos_UnpackExchangeFunctor<LMPDeviceType,0> f(atomKK,k_buf,k_count,k_indices,dim,lo,hi);
      Kokkos::parallel_for(nrecv/size_exchange,f);
    }
    k_count.modify_device();
    k_count.sync_host();
  }

  atomKK->modified(space,X_MASK | V_MASK | TAG_MASK | TYPE_MASK |
                 MASK_MASK | IMAGE_MASK| RADIUS_MASK | RMASS_MASK |
                 OMEGA_MASK);

  return k_count.view_host()(0);
}

/* ---------------------------------------------------------------------- */

void AtomVecSphereKokkos::sync(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.sync_device();
    if (mask & V_MASK) atomKK->k_v.sync_device();
    if (mask & F_MASK) atomKK->k_f.sync_device();
    if (mask & TAG_MASK) atomKK->k_tag.sync_device();
    if (mask & TYPE_MASK) atomKK->k_type.sync_device();
    if (mask & MASK_MASK) atomKK->k_mask.sync_device();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_device();
    if (mask & RADIUS_MASK) atomKK->k_radius.sync_device();
    if (mask & RMASS_MASK) atomKK->k_rmass.sync_device();
    if (mask & OMEGA_MASK) atomKK->k_omega.sync_device();
    if (mask & TORQUE_MASK) atomKK->k_torque.sync_device();
  } else if (space == Host) {
    if (mask & X_MASK) atomKK->k_x.sync_host();
    if (mask & V_MASK) atomKK->k_v.sync_host();
    if (mask & F_MASK) atomKK->k_f.sync_host();
    if (mask & TAG_MASK) atomKK->k_tag.sync_host();
    if (mask & TYPE_MASK) atomKK->k_type.sync_host();
    if (mask & MASK_MASK) atomKK->k_mask.sync_host();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_host();
    if (mask & RADIUS_MASK) atomKK->k_radius.sync_host();
    if (mask & RMASS_MASK) atomKK->k_rmass.sync_host();
    if (mask & OMEGA_MASK) atomKK->k_omega.sync_host();
    if (mask & TORQUE_MASK) atomKK->k_torque.sync_host();
  } else if (space == HostKK) {
    if (mask & X_MASK) atomKK->k_x.sync_hostkk();
    if (mask & V_MASK) atomKK->k_v.sync_hostkk();
    if (mask & F_MASK) atomKK->k_f.sync_hostkk();
    if (mask & TAG_MASK) atomKK->k_tag.sync_host();
    if (mask & TYPE_MASK) atomKK->k_type.sync_host();
    if (mask & MASK_MASK) atomKK->k_mask.sync_host();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_host();
    if (mask & RADIUS_MASK) atomKK->k_radius.sync_hostkk();
    if (mask & RMASS_MASK) atomKK->k_rmass.sync_hostkk();
    if (mask & OMEGA_MASK) atomKK->k_omega.sync_hostkk();
    if (mask & TORQUE_MASK) atomKK->k_torque.sync_hostkk();
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecSphereKokkos::sync_pinned(ExecutionSpace space, unsigned int mask, int async_flag)
{
  if (space == Device) {
    if ((mask & X_MASK) && atomKK->k_x.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3_lr>(atomKK->k_x,space,async_flag);
    if ((mask & V_MASK) && atomKK->k_v.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3>(atomKK->k_v,space,async_flag);
    if ((mask & F_MASK) && atomKK->k_f.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkacc_1d_3>(atomKK->k_f,space,async_flag);
    if ((mask & TAG_MASK) && atomKK->k_tag.need_sync_device())
      perform_pinned_copy<DAT::tdual_tagint_1d>(atomKK->k_tag,space,async_flag);
    if ((mask & TYPE_MASK) && atomKK->k_type.need_sync_device())
      perform_pinned_copy<DAT::tdual_int_1d>(atomKK->k_type,space,async_flag);
    if ((mask & MASK_MASK) && atomKK->k_mask.need_sync_device())
      perform_pinned_copy<DAT::tdual_int_1d>(atomKK->k_mask,space,async_flag);
    if ((mask & IMAGE_MASK) && atomKK->k_image.need_sync_device())
      perform_pinned_copy<DAT::tdual_imageint_1d>(atomKK->k_image,space,async_flag);
    if ((mask & RADIUS_MASK) && atomKK->k_radius.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_radius,space,async_flag);
    if ((mask & RMASS_MASK) && atomKK->k_rmass.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_rmass,space,async_flag);
    if ((mask & OMEGA_MASK) && atomKK->k_omega.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3>(atomKK->k_omega,space,async_flag);
    if ((mask & TORQUE_MASK) && atomKK->k_torque.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3>(atomKK->k_torque,space,async_flag);
  } else {
    if ((mask & X_MASK) && atomKK->k_x.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3_lr>(atomKK->k_x,space,async_flag);
    if ((mask & V_MASK) && atomKK->k_v.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3>(atomKK->k_v,space,async_flag);
    if ((mask & F_MASK) && atomKK->k_f.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkacc_1d_3>(atomKK->k_f,space,async_flag);
    if ((mask & TAG_MASK) && atomKK->k_tag.need_sync_host())
      perform_pinned_copy<DAT::tdual_tagint_1d>(atomKK->k_tag,space,async_flag);
    if ((mask & TYPE_MASK) && atomKK->k_type.need_sync_host())
      perform_pinned_copy<DAT::tdual_int_1d>(atomKK->k_type,space,async_flag);
    if ((mask & MASK_MASK) && atomKK->k_mask.need_sync_host())
      perform_pinned_copy<DAT::tdual_int_1d>(atomKK->k_mask,space,async_flag);
    if ((mask & IMAGE_MASK) && atomKK->k_image.need_sync_host())
      perform_pinned_copy<DAT::tdual_imageint_1d>(atomKK->k_image,space,async_flag);
    if ((mask & RADIUS_MASK) && atomKK->k_radius.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_radius,space,async_flag);
    if ((mask & RMASS_MASK) && atomKK->k_rmass.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_rmass,space,async_flag);
    if ((mask & OMEGA_MASK) && atomKK->k_omega.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3>(atomKK->k_omega,space,async_flag);
    if ((mask & TORQUE_MASK) && atomKK->k_torque.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3>(atomKK->k_torque,space,async_flag);
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecSphereKokkos::modified(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.modify_device();
    if (mask & V_MASK) atomKK->k_v.modify_device();
    if (mask & F_MASK) atomKK->k_f.modify_device();
    if (mask & TAG_MASK) atomKK->k_tag.modify_device();
    if (mask & TYPE_MASK) atomKK->k_type.modify_device();
    if (mask & MASK_MASK) atomKK->k_mask.modify_device();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_device();
    if (mask & RADIUS_MASK) atomKK->k_radius.modify_device();
    if (mask & RMASS_MASK) atomKK->k_rmass.modify_device();
    if (mask & OMEGA_MASK) atomKK->k_omega.modify_device();
    if (mask & TORQUE_MASK) atomKK->k_torque.modify_device();
  } else if (space == Host) {
    if (mask & X_MASK) atomKK->k_x.modify_host();
    if (mask & V_MASK) atomKK->k_v.modify_host();
    if (mask & F_MASK) atomKK->k_f.modify_host();
    if (mask & TAG_MASK) atomKK->k_tag.modify_host();
    if (mask & TYPE_MASK) atomKK->k_type.modify_host();
    if (mask & MASK_MASK) atomKK->k_mask.modify_host();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_host();
    if (mask & RADIUS_MASK) atomKK->k_radius.modify_host();
    if (mask & RMASS_MASK) atomKK->k_rmass.modify_host();
    if (mask & OMEGA_MASK) atomKK->k_omega.modify_host();
    if (mask & TORQUE_MASK) atomKK->k_torque.modify_host();
  } else if (space == HostKK) {
    if (mask & X_MASK) atomKK->k_x.modify_hostkk();
    if (mask & V_MASK) atomKK->k_v.modify_hostkk();
    if (mask & F_MASK) atomKK->k_f.modify_hostkk();
    if (mask & TAG_MASK) atomKK->k_tag.modify_host();
    if (mask & TYPE_MASK) atomKK->k_type.modify_host();
    if (mask & MASK_MASK) atomKK->k_mask.modify_host();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_host();
    if (mask & RADIUS_MASK) atomKK->k_radius.modify_hostkk();
    if (mask & RMASS_MASK) atomKK->k_rmass.modify_hostkk();
    if (mask & OMEGA_MASK) atomKK->k_omega.modify_hostkk();
    if (mask & TORQUE_MASK) atomKK->k_torque.modify_hostkk();
  }
}
