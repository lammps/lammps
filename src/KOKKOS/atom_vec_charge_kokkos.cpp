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

#include "atom_vec_charge_kokkos.h"
#include "atom_kokkos.h"
#include "comm_kokkos.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "atom_masks.h"
#include "memory_kokkos.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

#define DELTA 10

/* ---------------------------------------------------------------------- */

AtomVecChargeKokkos::AtomVecChargeKokkos(LAMMPS *lmp) : AtomVecKokkos(lmp)
{
  molecular = 0;
  mass_type = 1;

  comm_x_only = comm_f_only = 1;
  size_forward = 3;
  size_reverse = 3;
  size_border = 7;
  size_velocity = 3;
  size_data_atom = 6;
  size_data_vel = 4;
  xcol_data = 4;

  atom->q_flag = 1;

  k_count = DAT::tdual_int_1d("atom::k_count",1);
  atomKK = (AtomKokkos *) atom;
  commKK = (CommKokkos *) comm;

}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecChargeKokkos::grow(int n)
{
  int step = MAX(DELTA,nmax*0.01);
  if (n == 0) nmax += step;
  else nmax = n;
  atomKK->nmax = nmax;
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

  memoryKK->grow_kokkos(atomKK->k_q,atomKK->q,nmax,"atom:q");

  grow_pointers();
  atomKK->sync(Host,ALL_MASK);

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecChargeKokkos::grow_pointers()
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

  q = atomKK->q;
  d_q = atomKK->k_q.d_view;
  h_q = atomKK->k_q.h_view;

}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecChargeKokkos::copy(int i, int j, int delflag)
{
  h_tag[j] = h_tag[i];
  h_type[j] = h_type[i];
  mask[j] = mask[i];
  h_image[j] = h_image[i];
  h_x(j,0) = h_x(i,0);
  h_x(j,1) = h_x(i,1);
  h_x(j,2) = h_x(i,2);
  h_v(j,0) = h_v(i,0);
  h_v(j,1) = h_v(i,1);
  h_v(j,2) = h_v(i,2);

  h_q[j] = h_q[i];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG,int TRICLINIC>
struct AtomVecChargeKokkos_PackComm {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_x_array_randomread _x;
  typename ArrayTypes<DeviceType>::t_xfloat_2d_um _buf;
  typename ArrayTypes<DeviceType>::t_int_2d_const _list;
  const int _iswap;
  X_FLOAT _xprd,_yprd,_zprd,_xy,_xz,_yz;
  X_FLOAT _pbc[6];

  AtomVecChargeKokkos_PackComm(
      const typename DAT::tdual_x_array &x,
      const typename DAT::tdual_xfloat_2d &buf,
      const typename DAT::tdual_int_2d &list,
      const int & iswap,
      const X_FLOAT &xprd, const X_FLOAT &yprd, const X_FLOAT &zprd,
      const X_FLOAT &xy, const X_FLOAT &xz, const X_FLOAT &yz, const int* const pbc):
      _x(x.view<DeviceType>()),_list(list.view<DeviceType>()),_iswap(iswap),
      _xprd(xprd),_yprd(yprd),_zprd(zprd),
      _xy(xy),_xz(xz),_yz(yz) {
        const size_t maxsend = (buf.view<DeviceType>().extent(0)*buf.view<DeviceType>().extent(1))/3;
        const size_t elements = 3;
        buffer_view<DeviceType>(_buf,buf,maxsend,elements);
        _pbc[0] = pbc[0]; _pbc[1] = pbc[1]; _pbc[2] = pbc[2];
        _pbc[3] = pbc[3]; _pbc[4] = pbc[4]; _pbc[5] = pbc[5];
  };

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
  }
};

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG>
struct AtomVecChargeKokkos_PackBorder {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_xfloat_2d _buf;
  const typename ArrayTypes<DeviceType>::t_int_2d_const _list;
  const int _iswap;
  const typename ArrayTypes<DeviceType>::t_x_array_randomread _x;
  const typename ArrayTypes<DeviceType>::t_tagint_1d _tag;
  const typename ArrayTypes<DeviceType>::t_int_1d _type;
  const typename ArrayTypes<DeviceType>::t_int_1d _mask;
  const typename ArrayTypes<DeviceType>::t_float_1d _q;
  X_FLOAT _dx,_dy,_dz;

  AtomVecChargeKokkos_PackBorder(
      const typename ArrayTypes<DeviceType>::t_xfloat_2d &buf,
      const typename ArrayTypes<DeviceType>::t_int_2d_const &list,
      const int & iswap,
      const typename ArrayTypes<DeviceType>::t_x_array &x,
      const typename ArrayTypes<DeviceType>::t_tagint_1d &tag,
      const typename ArrayTypes<DeviceType>::t_int_1d &type,
      const typename ArrayTypes<DeviceType>::t_int_1d &mask,
      const typename ArrayTypes<DeviceType>::t_float_1d &q,
      const X_FLOAT &dx, const X_FLOAT &dy, const X_FLOAT &dz):
  _buf(buf),_list(list),_iswap(iswap),
    _x(x),_tag(tag),_type(type),_mask(mask),_q(q),
    _dx(dx),_dy(dy),_dz(dz) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
      const int j = _list(_iswap,i);
      if (PBC_FLAG == 0) {
          _buf(i,0) = _x(j,0);
          _buf(i,1) = _x(j,1);
          _buf(i,2) = _x(j,2);
          _buf(i,3) = d_ubuf(_tag(j)).d;
          _buf(i,4) = d_ubuf(_type(j)).d;
          _buf(i,5) = d_ubuf(_mask(j)).d;
          _buf(i,6) = _q(j);
      } else {
          _buf(i,0) = _x(j,0) + _dx;
          _buf(i,1) = _x(j,1) + _dy;
          _buf(i,2) = _x(j,2) + _dz;
          _buf(i,3) = d_ubuf(_tag(j)).d;
          _buf(i,4) = d_ubuf(_type(j)).d;
          _buf(i,5) = d_ubuf(_mask(j)).d;
          _buf(i,6) = _q(j);
      }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecChargeKokkos::pack_border_kokkos(int n, DAT::tdual_int_2d k_sendlist, DAT::tdual_xfloat_2d buf,int iswap,
                               int pbc_flag, int *pbc, ExecutionSpace space)
{
  X_FLOAT dx,dy,dz;

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
      AtomVecChargeKokkos_PackBorder<LMPHostType,1> f(
        buf.view<LMPHostType>(), k_sendlist.view<LMPHostType>(),
        iswap,h_x,h_tag,h_type,h_mask,h_q,dx,dy,dz);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecChargeKokkos_PackBorder<LMPDeviceType,1> f(
        buf.view<LMPDeviceType>(), k_sendlist.view<LMPDeviceType>(),
        iswap,d_x,d_tag,d_type,d_mask,d_q,dx,dy,dz);
      Kokkos::parallel_for(n,f);
    }

  } else {
    dx = dy = dz = 0;
    if(space==Host) {
      AtomVecChargeKokkos_PackBorder<LMPHostType,0> f(
        buf.view<LMPHostType>(), k_sendlist.view<LMPHostType>(),
        iswap,h_x,h_tag,h_type,h_mask,h_q,dx,dy,dz);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecChargeKokkos_PackBorder<LMPDeviceType,0> f(
        buf.view<LMPDeviceType>(), k_sendlist.view<LMPDeviceType>(),
        iswap,d_x,d_tag,d_type,d_mask,d_q,dx,dy,dz);
      Kokkos::parallel_for(n,f);
    }
  }
  return n*size_border;
}

/* ---------------------------------------------------------------------- */

int AtomVecChargeKokkos::pack_border(int n, int *list, double *buf,
                               int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = h_x(j,0);
      buf[m++] = h_x(j,1);
      buf[m++] = h_x(j,2);
      buf[m++] = ubuf(h_tag(j)).d;
      buf[m++] = ubuf(h_type(j)).d;
      buf[m++] = ubuf(h_mask(j)).d;
      buf[m++] = h_q(j);
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
      buf[m++] = ubuf(h_tag(j)).d;
      buf[m++] = ubuf(h_type(j)).d;
      buf[m++] = ubuf(h_mask(j)).d;
      buf[m++] = h_q(j);
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecChargeKokkos::pack_border_vel(int n, int *list, double *buf,
                                   int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = h_x(j,0);
      buf[m++] = h_x(j,1);
      buf[m++] = h_x(j,2);
      buf[m++] = ubuf(h_tag(j)).d;
      buf[m++] = ubuf(h_type(j)).d;
      buf[m++] = ubuf(h_mask(j)).d;
      buf[m++] = h_q[j];
      buf[m++] = h_v(j,0);
      buf[m++] = h_v(j,1);
      buf[m++] = h_v(j,2);
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
        buf[m++] = ubuf(h_tag(j)).d;
        buf[m++] = ubuf(h_type(j)).d;
        buf[m++] = ubuf(h_mask(j)).d;
        buf[m++] = h_q[j];
        buf[m++] = h_v(j,0);
        buf[m++] = h_v(j,1);
        buf[m++] = h_v(j,2);
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
        buf[m++] = ubuf(h_tag(j)).d;
        buf[m++] = ubuf(h_type(j)).d;
        buf[m++] = ubuf(h_mask(j)).d;
        buf[m++] = h_q[j];
        if (mask[i] & deform_groupbit) {
          buf[m++] = h_v(j,0) + dvx;
          buf[m++] = h_v(j,1) + dvy;
          buf[m++] = h_v(j,2) + dvz;
        } else {
          buf[m++] = h_v(j,0);
          buf[m++] = h_v(j,1);
          buf[m++] = h_v(j,2);
        }
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecChargeKokkos::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = h_q[j];
  }
  return m;
}


/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecChargeKokkos_UnpackBorder {
  typedef DeviceType device_type;

  const typename ArrayTypes<DeviceType>::t_xfloat_2d_const _buf;
  typename ArrayTypes<DeviceType>::t_x_array _x;
  typename ArrayTypes<DeviceType>::t_tagint_1d _tag;
  typename ArrayTypes<DeviceType>::t_int_1d _type;
  typename ArrayTypes<DeviceType>::t_int_1d _mask;
  typename ArrayTypes<DeviceType>::t_float_1d _q;
  int _first;


  AtomVecChargeKokkos_UnpackBorder(
      const typename ArrayTypes<DeviceType>::t_xfloat_2d_const &buf,
      typename ArrayTypes<DeviceType>::t_x_array &x,
      typename ArrayTypes<DeviceType>::t_tagint_1d &tag,
      typename ArrayTypes<DeviceType>::t_int_1d &type,
      typename ArrayTypes<DeviceType>::t_int_1d &mask,
      typename ArrayTypes<DeviceType>::t_float_1d &q,
      const int& first):
    _buf(buf),_x(x),_tag(tag),_type(type),_mask(mask),_q(q),_first(first){
  };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
      _x(i+_first,0) = _buf(i,0);
      _x(i+_first,1) = _buf(i,1);
      _x(i+_first,2) = _buf(i,2);
      _tag(i+_first) = (tagint) d_ubuf(_buf(i,3)).i;
      _type(i+_first) = (int) d_ubuf(_buf(i,4)).i;
      _mask(i+_first) = (int) d_ubuf(_buf(i,5)).i;
      _q(i+_first) = _buf(i,6);
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecChargeKokkos::unpack_border_kokkos(const int &n, const int &first,
                     const DAT::tdual_xfloat_2d &buf,ExecutionSpace space) {
  if (first+n >= nmax) {
    grow(first+n+100);
  }
  if(space==Host) {
    struct AtomVecChargeKokkos_UnpackBorder<LMPHostType>
      f(buf.view<LMPHostType>(),h_x,h_tag,h_type,h_mask,h_q,first);
    Kokkos::parallel_for(n,f);
  } else {
    struct AtomVecChargeKokkos_UnpackBorder<LMPDeviceType>
      f(buf.view<LMPDeviceType>(),d_x,d_tag,d_type,d_mask,d_q,first);
    Kokkos::parallel_for(n,f);
  }
  atomKK->modified(space,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|Q_MASK);
}

/* ---------------------------------------------------------------------- */

void AtomVecChargeKokkos::unpack_border(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  for (i = first; i < last; i++) {
    if (i == nmax) {
      grow(0);
    }
    atomKK->modified(Host,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|Q_MASK);
    h_x(i,0) = buf[m++];
    h_x(i,1) = buf[m++];
    h_x(i,2) = buf[m++];
    h_tag(i) =  (tagint)  ubuf(buf[m++]).i;
    h_type(i) = (int) ubuf(buf[m++]).i;
    h_mask(i) = (int) ubuf(buf[m++]).i;
    h_q[i] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecChargeKokkos::unpack_border_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    atomKK->modified(Host,X_MASK|V_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|Q_MASK);
    h_x(i,0) = buf[m++];
    h_x(i,1) = buf[m++];
    h_x(i,2) = buf[m++];
    h_tag(i) =  (tagint)  ubuf(buf[m++]).i;
    h_type(i) = (int) ubuf(buf[m++]).i;
    h_mask(i) = (int) ubuf(buf[m++]).i;
    h_q[i] = buf[m++];
    h_v(i,0) = buf[m++];
    h_v(i,1) = buf[m++];
    h_v(i,2) = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

int AtomVecChargeKokkos::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    h_q[i] = buf[m++];
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecChargeKokkos_PackExchangeFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_x_array_randomread _x;
  typename AT::t_v_array_randomread _v;
  typename AT::t_tagint_1d_randomread _tag;
  typename AT::t_int_1d_randomread _type;
  typename AT::t_int_1d_randomread _mask;
  typename AT::t_imageint_1d_randomread _image;
  typename AT::t_float_1d_randomread _q;
  typename AT::t_x_array _xw;
  typename AT::t_v_array _vw;
  typename AT::t_tagint_1d _tagw;
  typename AT::t_int_1d _typew;
  typename AT::t_int_1d _maskw;
  typename AT::t_imageint_1d _imagew;
  typename AT::t_float_1d _qw;

  typename AT::t_xfloat_2d_um _buf;
  typename AT::t_int_1d_const _sendlist;
  typename AT::t_int_1d_const _copylist;
  int _nlocal,_dim;
  X_FLOAT _lo,_hi;

  AtomVecChargeKokkos_PackExchangeFunctor(
      const AtomKokkos* atom,
      const typename AT::tdual_xfloat_2d buf,
      typename AT::tdual_int_1d sendlist,
      typename AT::tdual_int_1d copylist,int nlocal, int dim,
                X_FLOAT lo, X_FLOAT hi):
    _x(atom->k_x.view<DeviceType>()),
    _v(atom->k_v.view<DeviceType>()),
    _tag(atom->k_tag.view<DeviceType>()),
    _type(atom->k_type.view<DeviceType>()),
    _mask(atom->k_mask.view<DeviceType>()),
    _image(atom->k_image.view<DeviceType>()),
    _q(atom->k_q.view<DeviceType>()),
    _xw(atom->k_x.view<DeviceType>()),
    _vw(atom->k_v.view<DeviceType>()),
    _tagw(atom->k_tag.view<DeviceType>()),
    _typew(atom->k_type.view<DeviceType>()),
    _maskw(atom->k_mask.view<DeviceType>()),
    _imagew(atom->k_image.view<DeviceType>()),
    _qw(atom->k_q.view<DeviceType>()),
    _sendlist(sendlist.template view<DeviceType>()),
    _copylist(copylist.template view<DeviceType>()),
    _nlocal(nlocal),_dim(dim),
    _lo(lo),_hi(hi){
    const size_t elements = 12;
    const int maxsendlist = (buf.template view<DeviceType>().extent(0)*
                             buf.template view<DeviceType>().extent(1))/elements;

    buffer_view<DeviceType>(_buf,buf,maxsendlist,elements);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &mysend) const {
    const int i = _sendlist(mysend);
    _buf(mysend,0) = 12;
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
    _buf(mysend,11) = _q[i];
    const int j = _copylist(mysend);

    if(j>-1) {
    _xw(i,0) = _x(j,0);
    _xw(i,1) = _x(j,1);
    _xw(i,2) = _x(j,2);
    _vw(i,0) = _v(j,0);
    _vw(i,1) = _v(j,1);
    _vw(i,2) = _v(j,2);
    _tagw(i) = _tag(j);
    _typew(i) = _type(j);
    _maskw(i) = _mask(j);
    _imagew(i) = _image(j);
    _qw(i) = _q(j);
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecChargeKokkos::pack_exchange_kokkos(const int &nsend,DAT::tdual_xfloat_2d &k_buf,
                                              DAT::tdual_int_1d k_sendlist,
                                              DAT::tdual_int_1d k_copylist,
                                              ExecutionSpace space,int dim,
                                              X_FLOAT lo,X_FLOAT hi )
{
  if(nsend > (int) (k_buf.view<LMPHostType>().extent(0)*k_buf.view<LMPHostType>().extent(1))/12) {
    int newsize = nsend*12/k_buf.view<LMPHostType>().extent(1)+1;
    k_buf.resize(newsize,k_buf.view<LMPHostType>().extent(1));
  }
  if(space == Host) {
    AtomVecChargeKokkos_PackExchangeFunctor<LMPHostType>
      f(atomKK,k_buf,k_sendlist,k_copylist,atom->nlocal,dim,lo,hi);
    Kokkos::parallel_for(nsend,f);
    return nsend*12;
  } else {
    AtomVecChargeKokkos_PackExchangeFunctor<LMPDeviceType>
      f(atomKK,k_buf,k_sendlist,k_copylist,atom->nlocal,dim,lo,hi);
    Kokkos::parallel_for(nsend,f);
    return nsend*12;
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecChargeKokkos::pack_exchange(int i, double *buf)
{
  int m = 1;
  buf[m++] = h_x(i,0);
  buf[m++] = h_x(i,1);
  buf[m++] = h_x(i,2);
  buf[m++] = h_v(i,0);
  buf[m++] = h_v(i,1);
  buf[m++] = h_v(i,2);
  buf[m++] = ubuf(h_tag(i)).d;
  buf[m++] = ubuf(h_type(i)).d;
  buf[m++] = ubuf(h_mask(i)).d;
  buf[m++] = ubuf(h_image(i)).d;
  buf[m++] = h_q[i];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecChargeKokkos_UnpackExchangeFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_x_array _x;
  typename AT::t_v_array _v;
  typename AT::t_tagint_1d _tag;
  typename AT::t_int_1d _type;
  typename AT::t_int_1d _mask;
  typename AT::t_imageint_1d _image;
  typename AT::t_float_1d _q;
  typename AT::t_xfloat_2d_um _buf;
  typename AT::t_int_1d _nlocal;
  int _dim;
  X_FLOAT _lo,_hi;

  AtomVecChargeKokkos_UnpackExchangeFunctor(
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
    _q(atom->k_q.view<DeviceType>()),
    _nlocal(nlocal.template view<DeviceType>()),_dim(dim),
    _lo(lo),_hi(hi){
    const size_t elements = 12;
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
      _q[i] = _buf(myrecv,11);
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecChargeKokkos::unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf,int nrecv,
                                                int nlocal,int dim,X_FLOAT lo,X_FLOAT hi,
                                                ExecutionSpace space) {
  if(space == Host) {
    k_count.h_view(0) = nlocal;
    AtomVecChargeKokkos_UnpackExchangeFunctor<LMPHostType> f(atomKK,k_buf,k_count,dim,lo,hi);
    Kokkos::parallel_for(nrecv/12,f);
    return k_count.h_view(0);
  } else {
    k_count.h_view(0) = nlocal;
    k_count.modify<LMPHostType>();
    k_count.sync<LMPDeviceType>();
    AtomVecChargeKokkos_UnpackExchangeFunctor<LMPDeviceType>
      f(atomKK,k_buf,k_count,dim,lo,hi);
    Kokkos::parallel_for(nrecv/12,f);
    k_count.modify<LMPDeviceType>();
    k_count.sync<LMPHostType>();

    return k_count.h_view(0);
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecChargeKokkos::unpack_exchange(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);
  atomKK->modified(Host,X_MASK | V_MASK | TAG_MASK | TYPE_MASK |
           MASK_MASK | IMAGE_MASK | Q_MASK);

  int m = 1;
  h_x(nlocal,0) = buf[m++];
  h_x(nlocal,1) = buf[m++];
  h_x(nlocal,2) = buf[m++];
  h_v(nlocal,0) = buf[m++];
  h_v(nlocal,1) = buf[m++];
  h_v(nlocal,2) = buf[m++];
  h_tag(nlocal) = (tagint) ubuf(buf[m++]).i;
  h_type(nlocal) = (int) ubuf(buf[m++]).i;
  h_mask(nlocal) = (int) ubuf(buf[m++]).i;
  h_image(nlocal) = (imageint) ubuf(buf[m++]).i;
  h_q[nlocal] = buf[m++];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->
        unpack_exchange(nlocal,&buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecChargeKokkos::size_restart()
{
  int i;

  int nlocal = atom->nlocal;
  int n = 12 * nlocal;

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

int AtomVecChargeKokkos::pack_restart(int i, double *buf)
{
  atomKK->sync(Host,X_MASK | V_MASK | TAG_MASK | TYPE_MASK |
            MASK_MASK | IMAGE_MASK | Q_MASK);

  int m = 1;
  buf[m++] = h_x(i,0);
  buf[m++] = h_x(i,1);
  buf[m++] = h_x(i,2);
  buf[m++] = ubuf(h_tag(i)).d;
  buf[m++] = ubuf(h_type(i)).d;
  buf[m++] = ubuf(h_mask(i)).d;
  buf[m++] = ubuf(h_image(i)).d;
  buf[m++] = h_v(i,0);
  buf[m++] = h_v(i,1);
  buf[m++] = h_v(i,2);

  buf[m++] = h_q[i];

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecChargeKokkos::unpack_restart(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra,nmax,atom->nextra_store,"atom:extra");
  }

  atomKK->modified(Host,X_MASK | V_MASK | TAG_MASK | TYPE_MASK |
           MASK_MASK | IMAGE_MASK | Q_MASK);

  int m = 1;
  h_x(nlocal,0) = buf[m++];
  h_x(nlocal,1) = buf[m++];
  h_x(nlocal,2) = buf[m++];
  h_tag(nlocal) = (tagint) ubuf(buf[m++]).i;
  h_type(nlocal) = (int) ubuf(buf[m++]).i;
  h_mask(nlocal) = (int) ubuf(buf[m++]).i;
  h_image(nlocal) = (imageint) ubuf(buf[m++]).i;
  h_v(nlocal,0) = buf[m++];
  h_v(nlocal,1) = buf[m++];
  h_v(nlocal,2) = buf[m++];

  h_q[nlocal] = buf[m++];

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecChargeKokkos::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    atomKK->modified(Host,ALL_MASK);
    grow(0);
  }
  atomKK->sync(Host,ALL_MASK);
  atomKK->modified(Host,ALL_MASK);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  h_x(nlocal,0) = coord[0];
  h_x(nlocal,1) = coord[1];
  h_x(nlocal,2) = coord[2];
  h_mask[nlocal] = 1;
  h_image[nlocal] = ((imageint) IMGMAX << IMG2BITS) |
    ((imageint) IMGMAX << IMGBITS) | IMGMAX;
  h_v(nlocal,0) = 0.0;
  h_v(nlocal,1) = 0.0;
  h_v(nlocal,2) = 0.0;

  h_q[nlocal] = 0.0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecChargeKokkos::data_atom(double *coord, imageint imagetmp,
                                    char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  h_tag[nlocal] = utils::inumeric(FLERR,values[0],true,lmp);
  h_type[nlocal] = utils::inumeric(FLERR,values[1],true,lmp);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  h_q[nlocal] = utils::numeric(FLERR,values[2],true,lmp);

  h_x(nlocal,0) = coord[0];
  h_x(nlocal,1) = coord[1];
  h_x(nlocal,2) = coord[2];

  h_image[nlocal] = imagetmp;

  h_mask[nlocal] = 1;
  h_v(nlocal,0) = 0.0;
  h_v(nlocal,1) = 0.0;
  h_v(nlocal,2) = 0.0;

  atomKK->modified(Host,ALL_MASK);

  atom->nlocal++;
}
/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecChargeKokkos::data_atom_hybrid(int nlocal, char **values)
{
  h_q[nlocal] = utils::numeric(FLERR,values[0],true,lmp);

  return 1;
}
/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecChargeKokkos::pack_data(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = h_tag[i];
    buf[i][1] = h_type[i];
    buf[i][2] = h_q[i];
    buf[i][3] = h_x(i,0);
    buf[i][4] = h_x(i,1);
    buf[i][5] = h_x(i,2);
    buf[i][6] = (h_image[i] & IMGMASK) - IMGMAX;
    buf[i][7] = (h_image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    buf[i][8] = (h_image[i] >> IMG2BITS) - IMGMAX;
  }
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecChargeKokkos::pack_data_hybrid(int i, double *buf)
{
  buf[0] = h_q[i];
  return 1;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecChargeKokkos::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,"%d %d %-1.16e %-1.16e %-1.16e %-1.16e %d %d %d\n",
            (int) buf[i][0],(int) buf[i][1],buf[i][2],buf[i][3],buf[i][4],buf[i][5],
            (int) buf[i][6],(int) buf[i][7],(int) buf[i][8]);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecChargeKokkos::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e",buf[0]);
  return 1;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecChargeKokkos::memory_usage()
{
  bigint bytes = 0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*commKK->nthreads,3);

  if (atom->memcheck("q")) bytes += memory->usage(q,nmax);

  return bytes;
}

/* ---------------------------------------------------------------------- */

void AtomVecChargeKokkos::sync(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.sync<LMPDeviceType>();
    if (mask & V_MASK) atomKK->k_v.sync<LMPDeviceType>();
    if (mask & F_MASK) atomKK->k_f.sync<LMPDeviceType>();
    if (mask & TAG_MASK) atomKK->k_tag.sync<LMPDeviceType>();
    if (mask & TYPE_MASK) atomKK->k_type.sync<LMPDeviceType>();
    if (mask & MASK_MASK) atomKK->k_mask.sync<LMPDeviceType>();
    if (mask & IMAGE_MASK) atomKK->k_image.sync<LMPDeviceType>();
    if (mask & Q_MASK) atomKK->k_q.sync<LMPDeviceType>();
  } else {
    if (mask & X_MASK) atomKK->k_x.sync<LMPHostType>();
    if (mask & V_MASK) atomKK->k_v.sync<LMPHostType>();
    if (mask & F_MASK) atomKK->k_f.sync<LMPHostType>();
    if (mask & TAG_MASK) atomKK->k_tag.sync<LMPHostType>();
    if (mask & TYPE_MASK) atomKK->k_type.sync<LMPHostType>();
    if (mask & MASK_MASK) atomKK->k_mask.sync<LMPHostType>();
    if (mask & IMAGE_MASK) atomKK->k_image.sync<LMPHostType>();
    if (mask & Q_MASK) atomKK->k_q.sync<LMPHostType>();
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecChargeKokkos::modified(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.modify<LMPDeviceType>();
    if (mask & V_MASK) atomKK->k_v.modify<LMPDeviceType>();
    if (mask & F_MASK) atomKK->k_f.modify<LMPDeviceType>();
    if (mask & TAG_MASK) atomKK->k_tag.modify<LMPDeviceType>();
    if (mask & TYPE_MASK) atomKK->k_type.modify<LMPDeviceType>();
    if (mask & MASK_MASK) atomKK->k_mask.modify<LMPDeviceType>();
    if (mask & IMAGE_MASK) atomKK->k_image.modify<LMPDeviceType>();
    if (mask & Q_MASK) atomKK->k_q.modify<LMPDeviceType>();
  } else {
    if (mask & X_MASK) atomKK->k_x.modify<LMPHostType>();
    if (mask & V_MASK) atomKK->k_v.modify<LMPHostType>();
    if (mask & F_MASK) atomKK->k_f.modify<LMPHostType>();
    if (mask & TAG_MASK) atomKK->k_tag.modify<LMPHostType>();
    if (mask & TYPE_MASK) atomKK->k_type.modify<LMPHostType>();
    if (mask & MASK_MASK) atomKK->k_mask.modify<LMPHostType>();
    if (mask & IMAGE_MASK) atomKK->k_image.modify<LMPHostType>();
    if (mask & Q_MASK) atomKK->k_q.modify<LMPHostType>();
  }
}

void AtomVecChargeKokkos::sync_overlapping_device(ExecutionSpace space, unsigned int mask)
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
    if ((mask & MOLECULE_MASK) && atomKK->k_q.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_float_1d>(atomKK->k_q,space);
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
    if ((mask & MOLECULE_MASK) && atomKK->k_q.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_float_1d>(atomKK->k_q,space);
  }
}

