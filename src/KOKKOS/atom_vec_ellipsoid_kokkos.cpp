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

#include "atom_vec_ellipsoid_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm_kokkos.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory.h"
#include "memory_kokkos.h"
#include "modify.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

AtomVecEllipsoidKokkos::AtomVecEllipsoidKokkos(LAMMPS *lmp) :
    AtomVec(lmp), AtomVecKokkos(lmp), AtomVecEllipsoid(lmp), torque(nullptr)
{
  size_exchange_bonus = 8;
  datamask_bonus = ELLIPSOID_MASK|BONUS_MASK;

  k_nghost_bonus = DAT::tdual_int_scalar("atomEllipKK:k_nghost_bonus");
  k_nlocal_bonus = DAT::tdual_int_scalar("atomEllipKK:k_nlocal_bonus");

  if (((sizeof(KK_FLOAT) != sizeof(double))) && (comm->me == 0))
    error->warning(FLERR,"AtomVecEllipsoidKokkos does not (yet) fully support "
       "KK_FLOAT within bonus struct data (shape, quat). Using double for these fields.");
}

/* ---------------------------------------------------------------------- */

AtomVecEllipsoidKokkos::~AtomVecEllipsoidKokkos()
{
  if (bonus_flag) {
    memoryKK->destroy_kokkos(k_bonus,bonus);
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::init()
{
  AtomVecEllipsoid::init();

  set_atom_masks();
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::grow(int n)
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

  memoryKK->grow_kokkos(atomKK->k_rmass,atomKK->rmass,nmax,"atom:rmass");
  memoryKK->grow_kokkos(atomKK->k_angmom,atomKK->angmom,nmax,"atom:angmom");
  memoryKK->grow_kokkos(atomKK->k_torque,atomKK->torque,nmax,"atom:torque");
  memoryKK->grow_kokkos(atomKK->k_ellipsoid,atomKK->ellipsoid,nmax,"atom:ellipsoid");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);

  grow_pointers();
  atomKK->sync(Host,ALL_MASK);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::grow_pointers()
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

  rmass = atomKK->rmass;
  d_rmass = atomKK->k_rmass.view_device();
  h_rmass = atomKK->k_rmass.view_hostkk();
  angmom = atomKK->angmom;
  d_angmom = atomKK->k_angmom.view_device();
  h_angmom = atomKK->k_angmom.view_hostkk();
  torque = atomKK->torque;
  d_torque = atomKK->k_torque.view_device();
  h_torque = atomKK->k_torque.view_hostkk();
  ellipsoid = atomKK->ellipsoid;
  d_ellipsoid= atomKK->k_ellipsoid.view_device();
  h_ellipsoid = atomKK->k_ellipsoid.view_host();
}

/* ----------------------------------------------------------------------
   grow bonus data structure
------------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::grow_bonus()
{
  nmax_bonus = grow_nmax_bonus(nmax_bonus);
  if (nmax_bonus < 0) error->one(FLERR, "Per-processor system is too big");

  atomKK->sync(Device,BONUS_MASK);
  atomKK->modified(Device,BONUS_MASK);

  memoryKK->grow_kokkos(k_bonus,bonus,nmax_bonus,"atom:bonus");

  atomKK->sync(Host,BONUS_MASK);
}

/* ----------------------------------------------------------------------
   sort atom arrays on device
------------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  atomKK->sync(Device, ALL_MASK & ~F_MASK & ~TORQUE_MASK);

  Sorter.sort(LMPDeviceType(), d_tag);
  Sorter.sort(LMPDeviceType(), d_type);
  Sorter.sort(LMPDeviceType(), d_mask);
  Sorter.sort(LMPDeviceType(), d_image);
  Sorter.sort(LMPDeviceType(), d_x);
  Sorter.sort(LMPDeviceType(), d_v);
  Sorter.sort(LMPDeviceType(), d_rmass);
  Sorter.sort(LMPDeviceType(), d_angmom);
  Sorter.sort(LMPDeviceType(), d_ellipsoid);
  Sorter.sort(LMPDeviceType(), d_bonus);

  atomKK->modified(Device, ALL_MASK & ~F_MASK & ~TORQUE_MASK);
}

/* ------------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecEllipsoidKokkos_PackCommBonus {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AtomVecEllipsoidKokkosBonusArray<DeviceType>::t_bonus_1d _bonus;
  typename AT::t_int_1d _ellipsoid;
  typename AT::t_double_2d_lr_um _buf;
  typename AT::t_int_1d_const _list;
  const int _offset;

  AtomVecEllipsoidKokkos_PackCommBonus(
    const AtomKokkos* atomKK,
    const typename DAT::tdual_double_2d_lr &buf,
    const typename DEllipsoidBonusAT::tdual_bonus_1d &bonus,
    const typename DAT::tdual_int_1d &list,
    const int &offset):
      _ellipsoid(atomKK->k_ellipsoid.view<DeviceType>()),
      _bonus(bonus.view<DeviceType>()),
      _list(list.view<DeviceType>()),
      _offset(offset) {
    const int size_forward = atomKK->avecKK->size_forward;
    const size_t maxsend = (buf.view<DeviceType>().extent(0)*buf.view<DeviceType>().extent(1))/size_forward;
    const size_t elements = size_forward;
    buffer_view<DeviceType>(_buf,buf,maxsend,elements);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(i);
    int m = _offset;
    if (_ellipsoid(j) >= 0) {
      _buf(i,m++) = _bonus(_ellipsoid(j)).quat[0];
      _buf(i,m++) = _bonus(_ellipsoid(j)).quat[1];
      _buf(i,m++) = _bonus(_ellipsoid(j)).quat[2];
      _buf(i,m++) = _bonus(_ellipsoid(j)).quat[3];
    }
  }
};

/* ------------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::pack_comm_bonus_kokkos(const int &n, const DAT::tdual_int_1d &list,
                                                    const DAT::tdual_double_2d_lr &buf, int vel_flag)
{
  int offset = size_forward - size_forward_bonus;
  if (vel_flag) offset += size_velocity;

  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->sync(HostKK,datamask_bonus);
    struct AtomVecEllipsoidKokkos_PackCommBonus<LMPHostType> f(atomKK,buf,k_bonus,list,offset);
    Kokkos::parallel_for(n,f);
  } else {
    atomKK->sync(Device,datamask_bonus);
    struct AtomVecEllipsoidKokkos_PackCommBonus<LMPDeviceType> f(atomKK,buf,k_bonus,list,offset);
    Kokkos::parallel_for(n,f);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecEllipsoidKokkos_UnpackCommBonus {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AtomVecEllipsoidKokkosBonusArray<DeviceType>::t_bonus_1d _bonus;
  typename AT::t_int_1d _ellipsoid;
  typename AT::t_double_2d_lr_const _buf;
  const int _first;
  const int _offset;

  AtomVecEllipsoidKokkos_UnpackCommBonus(
    const AtomKokkos* atomKK,
    const typename DAT::tdual_double_2d_lr &buf,
    const typename DEllipsoidBonusAT::tdual_bonus_1d &bonus,
    const int& first,
    const int& offset):
      _bonus(bonus.view<DeviceType>()),
      _ellipsoid(atomKK->k_ellipsoid.view<DeviceType>()),
      _first(first),
      _offset(offset) {
    const int size_forward = atomKK->avecKK->size_forward;
    const size_t maxsend = (buf.view<DeviceType>().extent(0)*buf.view<DeviceType>().extent(1))/size_forward;
    const size_t elements = size_forward;
    buffer_view<DeviceType>(_buf,buf,maxsend,elements);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    if (_ellipsoid(i+_first) >= 0) {
      int m = _offset;
      _bonus(_ellipsoid(i+_first)).quat[0] = _buf(i,m++);
      _bonus(_ellipsoid(i+_first)).quat[1] = _buf(i,m++);
      _bonus(_ellipsoid(i+_first)).quat[2] = _buf(i,m++);
      _bonus(_ellipsoid(i+_first)).quat[3] = _buf(i,m++);
    }
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::unpack_comm_bonus_kokkos(const int &n, const int &first,
                                                      const DAT::tdual_double_2d_lr &buf, int vel_flag)
{
  int offset = size_forward - size_forward_bonus;
  if (vel_flag) offset += size_velocity;

  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->sync(HostKK,datamask_bonus);
    struct AtomVecEllipsoidKokkos_UnpackCommBonus<LMPHostType> f(
      atomKK,buf,k_bonus,first,offset);
    Kokkos::parallel_for(n,f);
    atomKK->modified(HostKK,datamask_bonus);
  } else {
    atomKK->sync(Device,datamask_bonus);
    struct AtomVecEllipsoidKokkos_UnpackCommBonus<LMPDeviceType> f(
      atomKK,buf,k_bonus,first,offset);
    Kokkos::parallel_for(n,f);
    atomKK->modified(Device,datamask_bonus);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecEllipsoidKokkos_PackCommSelfBonus {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AtomVecEllipsoidKokkosBonusArray<DeviceType>::t_bonus_1d _bonus;
  typename AT::t_int_1d _ellipsoid;

  int _nfirst;
  typename AT::t_int_1d_const _list;

  AtomVecEllipsoidKokkos_PackCommSelfBonus(
    const AtomKokkos* atomKK,
    const typename DEllipsoidBonusAT::tdual_bonus_1d &bonus,
    const int &nfirst,
    const typename DAT::tdual_int_1d &list):
      _ellipsoid(atomKK->k_ellipsoid.view<DeviceType>()),
      _bonus(bonus.view<DeviceType>()),
    _nfirst(nfirst),_list(list.view<DeviceType>()) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(i);
    if (_ellipsoid(i+_nfirst) >= 0 && _ellipsoid(j) >= 0) {
      _bonus(_ellipsoid(i+_nfirst)).quat[0] = _bonus(_ellipsoid(j)).quat[0];
      _bonus(_ellipsoid(i+_nfirst)).quat[1] = _bonus(_ellipsoid(j)).quat[1];
      _bonus(_ellipsoid(i+_nfirst)).quat[2] = _bonus(_ellipsoid(j)).quat[2];
      _bonus(_ellipsoid(i+_nfirst)).quat[3] = _bonus(_ellipsoid(j)).quat[3];
    }
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::pack_comm_self_bonus_kokkos(const int &n,
                                                         const DAT::tdual_int_1d &list,
                                                         const int nfirst) {
  // Check whether to always run forward communication on the host
  // Choose correct forward PackComm kernel

  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->sync(HostKK,datamask_bonus);
    struct AtomVecEllipsoidKokkos_PackCommSelfBonus<LMPHostType> f(
      atomKK,k_bonus,nfirst,list);
    Kokkos::parallel_for(n,f);
    atomKK->modified(HostKK,datamask_bonus);
  } else {
    atomKK->sync(Device,datamask_bonus);
    struct AtomVecEllipsoidKokkos_PackCommSelfBonus<LMPDeviceType> f(
      atomKK,k_bonus,nfirst,list);
    Kokkos::parallel_for(n,f);
    atomKK->modified(Device,datamask_bonus);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecEllipsoidKokkos_PackCommSelfFusedBonus {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AtomVecEllipsoidKokkosBonusArray<DeviceType>::t_bonus_1d _bonus;
  typename AT::t_int_1d _ellipsoid;

  typename AT::t_int_2d_lr_const _list;
  typename AT::t_int_1d_const _firstrecv;
  typename AT::t_int_1d_const _sendnum_scan;
  typename AT::t_int_1d_const _g2l;

  AtomVecEllipsoidKokkos_PackCommSelfFusedBonus(
    const AtomKokkos* atomKK,
    const typename DEllipsoidBonusAT::tdual_bonus_1d &bonus,
    const typename DAT::tdual_int_2d_lr &list,
    const typename DAT::tdual_int_1d &firstrecv,
    const typename DAT::tdual_int_1d &sendnum_scan,
    const typename DAT::tdual_int_1d &g2l):
      _ellipsoid(atomKK->k_ellipsoid.view<DeviceType>()),
      _bonus(bonus.view<DeviceType>()),
      _list(list.view<DeviceType>()),
      _firstrecv(firstrecv.view<DeviceType>()),
      _sendnum_scan(sendnum_scan.view<DeviceType>()),
      _g2l(g2l.view<DeviceType>()) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& ii) const {

    int iswap = 0;
    while (ii >= _sendnum_scan[iswap]) iswap++;
    int i = ii;
    if (iswap > 0)
      i = ii - _sendnum_scan[iswap-1];

    const int _nfirst = _firstrecv[iswap];
    const int nlocal = _firstrecv[0];

    int j = _list(iswap,i);
    if (j >= nlocal)
      j = _g2l(j-nlocal);

    if (_ellipsoid(i+_nfirst) >= 0 && _ellipsoid(j) >= 0) {
      _bonus(_ellipsoid(i+_nfirst)).quat[0] = _bonus(_ellipsoid(j)).quat[0];
      _bonus(_ellipsoid(i+_nfirst)).quat[1] = _bonus(_ellipsoid(j)).quat[1];
      _bonus(_ellipsoid(i+_nfirst)).quat[2] = _bonus(_ellipsoid(j)).quat[2];
      _bonus(_ellipsoid(i+_nfirst)).quat[3] = _bonus(_ellipsoid(j)).quat[3];
    }
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::pack_comm_self_fused_bonus_kokkos(const int &n,
                               const DAT::tdual_int_2d_lr &list,
                               const DAT::tdual_int_1d &sendnum_scan,
                               const DAT::tdual_int_1d &firstrecv,
                               const DAT::tdual_int_1d &g2l) {
  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->sync(HostKK,datamask_bonus);
    struct AtomVecEllipsoidKokkos_PackCommSelfFusedBonus<LMPHostType> f(
      atomKK,k_bonus,list,firstrecv,sendnum_scan,g2l);
    Kokkos::parallel_for(n,f);
    atomKK->modified(HostKK,datamask_bonus);
  } else {
    atomKK->sync(Device,datamask_bonus);
    struct AtomVecEllipsoidKokkos_PackCommSelfFusedBonus<LMPDeviceType> f(
      atomKK,k_bonus,list,firstrecv,sendnum_scan,g2l);
    Kokkos::parallel_for(n,f);
    atomKK->modified(Device,datamask_bonus);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecEllipsoidKokkos_PackBorderBonus {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_double_2d_lr _buf;
  const typename AT::t_int_1d_const _list;
  const typename AtomVecEllipsoidKokkosBonusArray<DeviceType>::t_bonus_1d_randomread _bonus;
  const typename AT::t_int_1d_randomread _ellipsoid;
  const int _offset;

  AtomVecEllipsoidKokkos_PackBorderBonus(
    const AtomKokkos* atomKK,
    const typename AT::t_double_2d_lr &buf,
    const typename AtomVecEllipsoidKokkosBonusArray<DeviceType>::t_bonus_1d &bonus,
    const typename AT::t_int_1d_const &list,
    const int &offset):
    _buf(buf),_list(list),_offset(offset),
    _bonus(bonus),
    _ellipsoid(atomKK->k_ellipsoid.view<DeviceType>()) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(i);
    const int j_bonus = _ellipsoid(j);
    int m = _offset;
    if (j_bonus < 0) {
      _buf(i,m) = d_ubuf(0).d;
    } else {
      _buf(i,m++) = d_ubuf(1).d;
      _buf(i,m++) = _bonus(j_bonus).shape[0];
      _buf(i,m++) = _bonus(j_bonus).shape[1];
      _buf(i,m++) = _bonus(j_bonus).shape[2];
      _buf(i,m++) = _bonus(j_bonus).quat[0];
      _buf(i,m++) = _bonus(j_bonus).quat[1];
      _buf(i,m++) = _bonus(j_bonus).quat[2];
      _buf(i,m++) = _bonus(j_bonus).quat[3];
    }
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::pack_border_bonus_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                                                      DAT::tdual_double_2d_lr &buf,
                                                      ExecutionSpace space, int vel_flag)
{
  int offset = size_border - size_border_bonus;
  if (vel_flag) offset += size_velocity;

  atomKK->sync(space,datamask_bonus);

  if (space == HostKK) {
    AtomVecEllipsoidKokkos_PackBorderBonus<LMPHostType> f(
      atomKK,buf.view_host(),k_bonus.view_host(),k_sendlist.view_host(),offset);
    Kokkos::parallel_for(n,f);
  } else {
    AtomVecEllipsoidKokkos_PackBorderBonus<LMPDeviceType> f(
      atomKK,buf.view_device(),k_bonus.view_device(),k_sendlist.view_device(),offset);
    Kokkos::parallel_for(n,f);
  }
}

/* ------------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecEllipsoidKokkos_UnpackBorderBonus {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_double_2d_lr_const _buf;
  typename AtomVecEllipsoidKokkosBonusArray<DeviceType>::t_bonus_1d _bonus;
  typename AT::t_int_1d _ellipsoid;
  const int _first;
  const int _offset;
  const int _nlocal_bonus;
  typename AT::t_int_scalar _nghost_bonus;

  AtomVecEllipsoidKokkos_UnpackBorderBonus(
    const AtomKokkos* atomKK,
    const typename AT::t_double_2d_lr_const &buf,
    const typename AtomVecEllipsoidKokkosBonusArray<DeviceType>::t_bonus_1d &bonus,
    const int& first,
    const int& offset,
    const int &nlocal_bonus,
    const typename AT::t_int_scalar &nghost_bonus):
      _buf(buf),
      _first(first),
      _offset(offset),
      _bonus(bonus),
      _ellipsoid(atomKK->k_ellipsoid.view<DeviceType>()),
      _nlocal_bonus(nlocal_bonus),
      _nghost_bonus(nghost_bonus) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {

    int ellipID = static_cast<int> (d_ubuf(_buf(i,_offset)).i);

    if (ellipID == 0 ) {
      _ellipsoid(i+_first) = -1;
    } else {
      int j = _nlocal_bonus + Kokkos::atomic_fetch_add(&_nghost_bonus(),1);
      int m = _offset + 1;
      _bonus(j).shape[0] = _buf(i,m++);
      _bonus(j).shape[1] = _buf(i,m++);
      _bonus(j).shape[2] = _buf(i,m++);
      _bonus(j).quat[0] = _buf(i,m++);
      _bonus(j).quat[1] = _buf(i,m++);
      _bonus(j).quat[2] = _buf(i,m++);
      _bonus(j).quat[3] = _buf(i,m++);
      _bonus(j).ilocal = i+_first;
      _ellipsoid(i+_first) = j;
    }
  }
};

/* ------------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::unpack_border_bonus_kokkos(const int &n, const int &first,
                                                        const DAT::tdual_double_2d_lr &buf,
                                                        ExecutionSpace space, int vel_flag) {
  while (first+n >= nmax) grow(0);
  while (n+nlocal_bonus+nghost_bonus >= nmax_bonus) grow_bonus();

  atomKK->sync(space,datamask_bonus);

  int offset = size_border - size_border_bonus;
  if (vel_flag) offset += size_velocity;

  if (space == HostKK) {
    k_nghost_bonus.view_host()() = nghost_bonus;
    struct AtomVecEllipsoidKokkos_UnpackBorderBonus<LMPHostType> f(
      atomKK,buf.view_host(),k_bonus.view_host(),first,offset,
      this->nlocal_bonus,k_nghost_bonus.view_host());
    Kokkos::parallel_for(n,f);
  } else {
    k_nghost_bonus.view_host()() = nghost_bonus;
    k_nghost_bonus.modify_host();
    k_nghost_bonus.sync_device();
    struct AtomVecEllipsoidKokkos_UnpackBorderBonus<LMPDeviceType> f(
      atomKK,buf.view_device(),k_bonus.view_device(),first,offset,
      this->nlocal_bonus, k_nghost_bonus.view_device());
    Kokkos::parallel_for(n,f);
    k_nghost_bonus.modify_device();
    k_nghost_bonus.sync_host();
  }

  atomKK->modified(space,datamask_bonus);

  nghost_bonus = k_nghost_bonus.view_host()();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecEllipsoidKokkos_PackExchangeBonus {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AtomVecEllipsoidKokkosBonusArray<DeviceType>::t_bonus_1d_randomread _bonus;
  typename AT::t_int_1d_randomread _ellipsoid;
  typename AtomVecEllipsoidKokkosBonusArray<DeviceType>::t_bonus_1d _bonusw;
  typename AT::t_int_1d _ellipsoidw;

  typename AT::t_double_2d_lr_um _buf;
  typename AT::t_int_1d_const _sendlist;
  typename AT::t_int_1d_const _copylist_bonus;
  const int _size_exchange;
  const int _offset;

  AtomVecEllipsoidKokkos_PackExchangeBonus(
    const AtomKokkos* atomKK,
    const DAT::tdual_double_2d_lr &buf,
    const typename DEllipsoidBonusAT::tdual_bonus_1d &bonus,
    DAT::tdual_int_1d &sendlist,
    DAT::tdual_int_1d &copylist_bonus,
    const int &offset):
      _bonus(bonus.template view<DeviceType>()),
      _ellipsoid(atomKK->k_ellipsoid.template view<DeviceType>()),
      _bonusw(bonus.template view<DeviceType>()),
      _ellipsoidw(atomKK->k_ellipsoid.template view<DeviceType>()),

      _size_exchange(atomKK->avecKK->size_exchange),
      _sendlist(sendlist.template view<DeviceType>()),
      _copylist_bonus(copylist_bonus.template view<DeviceType>()),
      _offset(offset) {
    const int maxsendlist = (buf.template view<DeviceType>().extent(0)*
                             buf.template view<DeviceType>().extent(1))/_size_exchange;
    buffer_view<DeviceType>(_buf,buf,maxsendlist,_size_exchange);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &mysend) const {
    const int i = _sendlist(mysend);

    int m = _offset;
    if (_ellipsoid[i] < 0)
      _buf(mysend,m++) = d_ubuf(0).d;
    else {
      _buf(mysend,m++) = d_ubuf(1).d;
      int j = _ellipsoid[i];
      _buf(mysend,m++) = _bonus(j).shape[0];
      _buf(mysend,m++) = _bonus(j).shape[1];
      _buf(mysend,m++) = _bonus(j).shape[2];
      _buf(mysend,m++) = _bonus(j).quat[0];
      _buf(mysend,m++) = _bonus(j).quat[1];
      _buf(mysend,m++) = _bonus(j).quat[2];
      _buf(mysend,m++) = _bonus(j).quat[3];
    }

    int i_bonus = _ellipsoid[i];
    int j_bonus = _copylist_bonus(mysend); // may be different than ellipsoid[j]

    if (j_bonus > -1) {

      // delete bonus data from i_bonus

      if (i_bonus > -1) {

        // copy bonus data from J to I, effectively deleting the I entry
        // also reset ellipsoid that points to J to now point to I

        _ellipsoidw[_bonus[j_bonus].ilocal] = i_bonus;
        _bonusw[i_bonus] = _bonus[j_bonus];
      }
    }
  }
};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecEllipsoidKokkos_BackfillEllipsoid {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AtomVecEllipsoidKokkosBonusArray<DeviceType>::t_bonus_1d_randomread _bonus;
  typename AT::t_int_1d_randomread _ellipsoid;
  typename AtomVecEllipsoidKokkosBonusArray<DeviceType>::t_bonus_1d _bonusw;
  typename AT::t_int_1d _ellipsoidw;

  typename AT::t_double_2d_lr_um _buf;
  typename AT::t_int_1d_const _sendlist;
  typename AT::t_int_1d_const _copylist;
  const int _size_exchange;
  const int _offset;

  AtomVecEllipsoidKokkos_BackfillEllipsoid(
    const AtomKokkos* atomKK,
    const DAT::tdual_double_2d_lr &buf,
    const typename DEllipsoidBonusAT::tdual_bonus_1d &bonus,
    DAT::tdual_int_1d &sendlist,
    DAT::tdual_int_1d &copylist,
    const int &offset):
      _bonus(bonus.template view<DeviceType>()),
      _ellipsoid(atomKK->k_ellipsoid.template view<DeviceType>()),
      _bonusw(bonus.template view<DeviceType>()),
      _ellipsoidw(atomKK->k_ellipsoid.template view<DeviceType>()),

      _size_exchange(atomKK->avecKK->size_exchange),
      _sendlist(sendlist.template view<DeviceType>()),
      _copylist(copylist.template view<DeviceType>()),
      _offset(offset) {
    const int maxsendlist = (buf.template view<DeviceType>().extent(0)*
                             buf.template view<DeviceType>().extent(1))/_size_exchange;
    buffer_view<DeviceType>(_buf,buf,maxsendlist,_size_exchange);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &mysend) const {
    const int i = _sendlist(mysend);

    // if atom J has bonus data, reset J’s bonus.ilocal to loc I

    int j = _copylist(mysend);
    if (j > -1) {
      if (_ellipsoid[j] >= 0) _bonusw[_ellipsoid[j]].ilocal = i;
      _ellipsoidw[i] = _ellipsoid[j];
    }
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::pack_exchange_bonus_kokkos(const int &nsend,
                                                        DAT::tdual_double_2d_lr &k_buf,
                                                        DAT::tdual_int_1d k_sendlist,
                                                        DAT::tdual_int_1d k_copylist,
                                                        DAT::tdual_int_1d k_copylist_bonus,
                                                        ExecutionSpace space)
{
  int offset = size_exchange - size_exchange_bonus;

  atomKK->sync(space,datamask_bonus);

  if (space == HostKK) {
    AtomVecEllipsoidKokkos_PackExchangeBonus<LMPHostType> f(atomKK,
      k_buf,k_bonus,
      k_sendlist,k_copylist_bonus,
      offset);
    Kokkos::parallel_for(nsend,f);

    // must backfill ellipsoid after pack exchange bonus in a separate
    //  functor to prevent race conditions

    AtomVecEllipsoidKokkos_BackfillEllipsoid<LMPHostType> f2(atomKK,
      k_buf,k_bonus,
      k_sendlist,k_copylist,
      offset);
    Kokkos::parallel_for(nsend,f2);
  } else {
    AtomVecEllipsoidKokkos_PackExchangeBonus<LMPDeviceType> f(atomKK,
      k_buf,k_bonus,
      k_sendlist,k_copylist_bonus,
      offset);
    Kokkos::parallel_for(nsend,f);

    // must backfill ellipsoid after pack exchange bonus in a separate
    //  functor to prevent race conditions

    AtomVecEllipsoidKokkos_BackfillEllipsoid<LMPDeviceType> f2(atomKK,
      k_buf,k_bonus,
      k_sendlist,k_copylist,
      offset);
    Kokkos::parallel_for(nsend,f2);
  }

  atomKK->modified(space,datamask_bonus);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecEllipsoidKokkos_UnpackExchangeBonus {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AtomVecEllipsoidKokkosBonusArray<DeviceType>::t_bonus_1d _bonus;
  typename AT::t_int_1d _ellipsoid;

  typename AT::t_double_2d_lr _buf;
  typename AT::t_int_1d _indices;
  typename AT::t_int_scalar _nlocal_bonus;
  int _size_exchange;
  const int _offset;

  AtomVecEllipsoidKokkos_UnpackExchangeBonus(
    const AtomKokkos* atomKK,
    const typename DAT::tdual_double_2d_lr &buf,
    const typename DEllipsoidBonusAT::tdual_bonus_1d &bonus,
    typename AT::tdual_int_scalar &nlocal_bonus,
    typename AT::tdual_int_1d &indices,
    const int &offset):
      _bonus(bonus.view<DeviceType>()),
      _ellipsoid(atomKK->k_ellipsoid.view<DeviceType>()),

      _size_exchange(atomKK->avecKK->size_exchange),
      _nlocal_bonus(nlocal_bonus.template view<DeviceType>()),
      _indices(indices.template view<DeviceType>()),
      _offset(offset) {
    const int maxsendlist = (buf.template view<DeviceType>().extent(0)*
                             buf.template view<DeviceType>().extent(1))/_size_exchange;
    buffer_view<DeviceType>(_buf,buf,maxsendlist,_size_exchange);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &myrecv) const {
    int i = _indices(myrecv);
    if (i > -1) {
      if ((int) d_ubuf(_buf(myrecv,_offset)).i == 0)
        _ellipsoid(i) = -1;
      else {
        int m = _offset + 1;
        int k = Kokkos::atomic_fetch_add(&_nlocal_bonus(),1);
        _bonus(k).shape[0] = _buf(myrecv,m++);
        _bonus(k).shape[1] = _buf(myrecv,m++);
        _bonus(k).shape[2] = _buf(myrecv,m++);
        _bonus(k).quat[0] = _buf(myrecv,m++);
        _bonus(k).quat[1] = _buf(myrecv,m++);
        _bonus(k).quat[2] = _buf(myrecv,m++);
        _bonus(k).quat[3] = _buf(myrecv,m++);
        _bonus(k).ilocal = i;
        _ellipsoid(i) = k;
      }
    }
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::unpack_exchange_bonus_kokkos(DAT::tdual_double_2d_lr &k_buf, int nrecv,
                                              ExecutionSpace space,
                                              DAT::tdual_int_1d &k_indices)
{
  while (nlocal_bonus + nrecv/size_exchange >= nmax_bonus) grow_bonus();

  int offset = size_exchange - size_exchange_bonus;

  atomKK->sync(space,datamask_bonus);

  if (space == HostKK) {
    k_nlocal_bonus.view_host()() = nlocal_bonus;

    AtomVecEllipsoidKokkos_UnpackExchangeBonus<LMPHostType> f(
      atomKK,k_buf,k_bonus,
      k_nlocal_bonus,k_indices,offset);
    Kokkos::parallel_for(nrecv/size_exchange,f);
  } else {
    k_nlocal_bonus.view_host()() = nlocal_bonus;
    k_nlocal_bonus.modify_host();
    k_nlocal_bonus.sync_device();

    struct AtomVecEllipsoidKokkos_UnpackExchangeBonus<LMPDeviceType> f(
      atomKK,k_buf,k_bonus,
      k_nlocal_bonus,k_indices,offset);
    Kokkos::parallel_for(nrecv/size_exchange,f);

    k_nlocal_bonus.modify_device();
    k_nlocal_bonus.sync_host();
  }

  atomKK->modified(space,datamask_bonus);

  nlocal_bonus = k_nlocal_bonus.view_host()();
}

/* ---------------------------------------------------------------------- */

int AtomVecEllipsoidKokkos::get_status_nlocal_bonus() {
  return nlocal_bonus;
}

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::set_status_nlocal_bonus(int nlocal_bonus) {
  this->nlocal_bonus = nlocal_bonus;
}

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::sync(ExecutionSpace space, uint64_t mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.sync_device();
    if (mask & V_MASK) atomKK->k_v.sync_device();
    if (mask & F_MASK) atomKK->k_f.sync_device();
    if (mask & TAG_MASK) atomKK->k_tag.sync_device();
    if (mask & TYPE_MASK) atomKK->k_type.sync_device();
    if (mask & MASK_MASK) atomKK->k_mask.sync_device();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_device();
    if (mask & RMASS_MASK) atomKK->k_rmass.sync_device();
    if (mask & ANGMOM_MASK) atomKK->k_angmom.sync_device();
    if (mask & TORQUE_MASK) atomKK->k_torque.sync_device();
    if (mask & ELLIPSOID_MASK) atomKK->k_ellipsoid.sync_device();
    if (mask & BONUS_MASK) k_bonus.sync_device();
  } else if (space == Host) {
    if (mask & X_MASK) atomKK->k_x.sync_host();
    if (mask & V_MASK) atomKK->k_v.sync_host();
    if (mask & F_MASK) atomKK->k_f.sync_host();
    if (mask & TAG_MASK) atomKK->k_tag.sync_host();
    if (mask & TYPE_MASK) atomKK->k_type.sync_host();
    if (mask & MASK_MASK) atomKK->k_mask.sync_host();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_host();
    if (mask & RMASS_MASK) atomKK->k_rmass.sync_host();
    if (mask & ANGMOM_MASK) atomKK->k_angmom.sync_host();
    if (mask & TORQUE_MASK) atomKK->k_torque.sync_host();
    if (mask & ELLIPSOID_MASK) atomKK->k_ellipsoid.sync_host();
    if (mask & BONUS_MASK) k_bonus.sync_host();
  } else if (space == HostKK){
    if (mask & X_MASK) atomKK->k_x.sync_hostkk();
    if (mask & V_MASK) atomKK->k_v.sync_hostkk();
    if (mask & F_MASK) atomKK->k_f.sync_hostkk();
    if (mask & TAG_MASK) atomKK->k_tag.sync_host();
    if (mask & TYPE_MASK) atomKK->k_type.sync_host();
    if (mask & MASK_MASK) atomKK->k_mask.sync_host();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_host();
    if (mask & RMASS_MASK) atomKK->k_rmass.sync_hostkk();
    if (mask & ANGMOM_MASK) atomKK->k_angmom.sync_hostkk();
    if (mask & TORQUE_MASK) atomKK->k_torque.sync_hostkk();
    if (mask & ELLIPSOID_MASK) atomKK->k_ellipsoid.sync_host();
    if (mask & BONUS_MASK) k_bonus.sync_host();
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::sync_pinned(ExecutionSpace space, uint64_t mask, int async_flag)
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
    if ((mask & RMASS_MASK) && atomKK->k_rmass.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_rmass,space,async_flag);
    if ((mask & ANGMOM_MASK) && atomKK->k_angmom.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3>(atomKK->k_angmom,space,async_flag);
    if ((mask & TORQUE_MASK) && atomKK->k_torque.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkacc_1d_3>(atomKK->k_torque,space,async_flag);
    if ((mask & ELLIPSOID_MASK) && atomKK->k_ellipsoid.need_sync_device())
      perform_pinned_copy<DAT::tdual_int_1d>(atomKK->k_ellipsoid,space,async_flag);
    if ((mask & BONUS_MASK) && k_bonus.need_sync_device())
      perform_pinned_copy<DEllipsoidBonusAT::tdual_bonus_1d>(k_bonus,space,async_flag);
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
    if ((mask & RMASS_MASK) && atomKK->k_rmass.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_rmass,space,async_flag);
    if ((mask & ANGMOM_MASK) && atomKK->k_angmom.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3>(atomKK->k_angmom,space,async_flag);
    if ((mask & TORQUE_MASK) && atomKK->k_torque.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkacc_1d_3>(atomKK->k_torque,space,async_flag);
    if ((mask & ELLIPSOID_MASK) && atomKK->k_ellipsoid.need_sync_host())
      perform_pinned_copy<DAT::tdual_int_1d>(atomKK->k_ellipsoid,space,async_flag);
    if ((mask & BONUS_MASK) && k_bonus.need_sync_host())
      perform_pinned_copy<DEllipsoidBonusAT::tdual_bonus_1d>(k_bonus,space,async_flag);
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::modified(ExecutionSpace space, uint64_t mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.modify_device();
    if (mask & V_MASK) atomKK->k_v.modify_device();
    if (mask & F_MASK) atomKK->k_f.modify_device();
    if (mask & TAG_MASK) atomKK->k_tag.modify_device();
    if (mask & TYPE_MASK) atomKK->k_type.modify_device();
    if (mask & MASK_MASK) atomKK->k_mask.modify_device();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_device();
    if (mask & RMASS_MASK) atomKK->k_rmass.modify_device();
    if (mask & ANGMOM_MASK) atomKK->k_angmom.modify_device();
    if (mask & TORQUE_MASK) atomKK->k_torque.modify_device();
    if (mask & ELLIPSOID_MASK) atomKK->k_ellipsoid.modify_device();
    if (mask & BONUS_MASK) k_bonus.modify_device();
  } else if (space == Host) {
    if (mask & X_MASK) atomKK->k_x.modify_host();
    if (mask & V_MASK) atomKK->k_v.modify_host();
    if (mask & F_MASK) atomKK->k_f.modify_host();
    if (mask & TAG_MASK) atomKK->k_tag.modify_host();
    if (mask & TYPE_MASK) atomKK->k_type.modify_host();
    if (mask & MASK_MASK) atomKK->k_mask.modify_host();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_host();
    if (mask & RMASS_MASK) atomKK->k_rmass.modify_host();
    if (mask & ANGMOM_MASK) atomKK->k_angmom.modify_host();
    if (mask & TORQUE_MASK) atomKK->k_torque.modify_host();
    if (mask & ELLIPSOID_MASK) atomKK->k_ellipsoid.modify_host();
    if (mask & BONUS_MASK) k_bonus.modify_host();
  } else if (space == HostKK) {
    if (mask & X_MASK) atomKK->k_x.modify_hostkk();
    if (mask & V_MASK) atomKK->k_v.modify_hostkk();
    if (mask & F_MASK) atomKK->k_f.modify_hostkk();
    if (mask & TAG_MASK) atomKK->k_tag.modify_host();
    if (mask & TYPE_MASK) atomKK->k_type.modify_host();
    if (mask & MASK_MASK) atomKK->k_mask.modify_host();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_host();
    if (mask & RMASS_MASK) atomKK->k_rmass.modify_hostkk();
    if (mask & ANGMOM_MASK) atomKK->k_angmom.modify_hostkk();
    if (mask & TORQUE_MASK) atomKK->k_torque.modify_hostkk();
    if (mask & ELLIPSOID_MASK) atomKK->k_ellipsoid.modify_host();
    if (mask & BONUS_MASK) k_bonus.modify_host();
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::set_size_exchange()
{
  AtomVecKokkos::set_size_exchange();
  size_exchange += size_exchange_bonus;
}
