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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)

   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018).
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics
   and molecular dynamics. Journal of Computational Physics.
------------------------------------------------------------------------- */

#include "atom_vec_spin_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "memory_kokkos.h"
#include "modify.h"

using namespace LAMMPS_NS;

static constexpr int DELTA = 10;

/* ---------------------------------------------------------------------- */

AtomVecSpinKokkos::AtomVecSpinKokkos(LAMMPS *lmp) : AtomVec(lmp),
AtomVecKokkos(lmp), AtomVecSpin(lmp)
{

}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecSpinKokkos::grow(int n)
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

  // allocating mech. quantities

  memoryKK->grow_kokkos(atomKK->k_x,atomKK->x,nmax,"atom:x");
  memoryKK->grow_kokkos(atomKK->k_v,atomKK->v,nmax,"atom:v");
  memoryKK->grow_kokkos(atomKK->k_f,atomKK->f,nmax,"atom:f");

  // allocating mag. quantities

  memoryKK->grow_kokkos(atomKK->k_sp,atomKK->sp,nmax,"atom:sp");
  memoryKK->grow_kokkos(atomKK->k_fm,atomKK->fm,nmax,"atom:fm");
  memoryKK->grow_kokkos(atomKK->k_fm_long,atomKK->fm_long,nmax,"atom:fm_long");

  grow_pointers();
  atomKK->sync(Host,ALL_MASK);

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecSpinKokkos::grow_pointers()
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

  sp = atomKK->sp;
  d_sp = atomKK->k_sp.view_device();
  h_sp = atomKK->k_sp.view_hostkk();
  fm = atomKK->fm;
  d_fm = atomKK->k_fm.view_device();
  h_fm = atomKK->k_fm.view_hostkk();
  fm_long = atomKK->fm_long;
  d_fm_long = atomKK->k_fm_long.view_device();
  h_fm_long = atomKK->k_fm_long.view_hostkk();
}

/* ----------------------------------------------------------------------
   sort atom arrays on device
------------------------------------------------------------------------- */

void AtomVecSpinKokkos::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  atomKK->sync(Device, TAG_MASK|TYPE_MASK|MASK_MASK|IMAGE_MASK|X_MASK|V_MASK|SP_MASK);

  Sorter.sort(LMPDeviceType(), d_tag);
  Sorter.sort(LMPDeviceType(), d_type);
  Sorter.sort(LMPDeviceType(), d_mask);
  Sorter.sort(LMPDeviceType(), d_image);
  Sorter.sort(LMPDeviceType(), d_x);
  Sorter.sort(LMPDeviceType(), d_v);
  Sorter.sort(LMPDeviceType(), d_sp);

  atomKK->modified(Device, TAG_MASK|TYPE_MASK|MASK_MASK|IMAGE_MASK|X_MASK|V_MASK|SP_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG,int TRICLINIC>
struct AtomVecSpinKokkos_PackComm {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr_randomread _x;
  typename AT::t_kkfloat_1d_4_randomread _sp;
  typename AT::t_double_2d_lr_um _buf;
  typename AT::t_int_1d_const _list;
  double _xprd,_yprd,_zprd,_xy,_xz,_yz;
  double _pbc[6];

  AtomVecSpinKokkos_PackComm(
      const typename DAT::ttransform_kkfloat_1d_3_lr &x,
      const typename DAT::ttransform_kkfloat_1d_4 &sp,
      const typename DAT::tdual_double_2d_lr &buf,
      const typename DAT::tdual_int_1d &list,
      const double &xprd, const double &yprd, const double &zprd,
      const double &xy, const double &xz, const double &yz, const int* const pbc):
      _x(x.view<DeviceType>()),_sp(sp.view<DeviceType>()),
      _list(list.view<DeviceType>()),
      _xprd(xprd),_yprd(yprd),_zprd(zprd),
      _xy(xy),_xz(xz),_yz(yz) {
        const size_t maxsend = (buf.view<DeviceType>().extent(0)*buf.view<DeviceType>().extent(1))/3;
        const size_t elements = 7;
        buffer_view<DeviceType>(_buf,buf,maxsend,elements);
        _pbc[0] = pbc[0]; _pbc[1] = pbc[1]; _pbc[2] = pbc[2];
        _pbc[3] = pbc[3]; _pbc[4] = pbc[4]; _pbc[5] = pbc[5];
  };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
        const int j = _list(i);
      if (PBC_FLAG == 0) {
          _buf(i,0) = _x(j,0);
          _buf(i,1) = _x(j,1);
          _buf(i,2) = _x(j,2);
          _buf(i,3) = _sp(j,0);
          _buf(i,4) = _sp(j,1);
          _buf(i,5) = _sp(j,2);
          _buf(i,6) = _sp(j,3);
      } else {
        if (TRICLINIC == 0) {
          _buf(i,0) = _x(j,0) + _pbc[0]*_xprd;
          _buf(i,1) = _x(j,1) + _pbc[1]*_yprd;
          _buf(i,2) = _x(j,2) + _pbc[2]*_zprd;
          _buf(i,3) = _sp(j,0);
          _buf(i,4) = _sp(j,1);
          _buf(i,5) = _sp(j,2);
          _buf(i,6) = _sp(j,3);
        } else {
          _buf(i,0) = _x(j,0) + _pbc[0]*_xprd + _pbc[5]*_xy + _pbc[4]*_xz;
          _buf(i,1) = _x(j,1) + _pbc[1]*_yprd + _pbc[3]*_yz;
          _buf(i,2) = _x(j,2) + _pbc[2]*_zprd;
          _buf(i,3) = _sp(j,0);
          _buf(i,4) = _sp(j,1);
          _buf(i,5) = _sp(j,2);
          _buf(i,6) = _sp(j,3);
        }
      }
  }
};

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG>
struct AtomVecSpinKokkos_PackBorder {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_double_2d_lr _buf;
  const typename AT::t_int_1d_const _list;
  const typename AT::t_kkfloat_1d_3_lr_randomread _x;
  const typename AT::t_tagint_1d _tag;
  const typename AT::t_int_1d _type;
  const typename AT::t_int_1d _mask;
  const typename AT::t_kkfloat_1d_4_randomread _sp;
  double _dx,_dy,_dz;

  AtomVecSpinKokkos_PackBorder(
      const typename AT::t_double_2d_lr &buf,
      const typename AT::t_int_1d_const &list,
      const typename AT::t_kkfloat_1d_3_lr &x,
      const typename AT::t_tagint_1d &tag,
      const typename AT::t_int_1d &type,
      const typename AT::t_int_1d &mask,
      const typename AT::t_kkfloat_1d_4 &sp,
      const double &dx, const double &dy, const double &dz):
  _buf(buf),_list(list),
    _x(x),_tag(tag),_type(type),_mask(mask),_sp(sp),
    _dx(dx),_dy(dy),_dz(dz) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
      const int j = _list(i);
      if (PBC_FLAG == 0) {
          _buf(i,0) = _x(j,0);
          _buf(i,1) = _x(j,1);
          _buf(i,2) = _x(j,2);
          _buf(i,3) = d_ubuf(_tag(j)).d;
          _buf(i,4) = d_ubuf(_type(j)).d;
          _buf(i,5) = d_ubuf(_mask(j)).d;
          _buf(i,6) = _sp(j,0);
          _buf(i,7) = _sp(j,1);
          _buf(i,8) = _sp(j,2);
          _buf(i,9) = _sp(j,3);
      } else {
          _buf(i,0) = _x(j,0) + _dx;
          _buf(i,1) = _x(j,1) + _dy;
          _buf(i,2) = _x(j,2) + _dz;
          _buf(i,3) = d_ubuf(_tag(j)).d;
          _buf(i,4) = d_ubuf(_type(j)).d;
          _buf(i,5) = d_ubuf(_mask(j)).d;
          _buf(i,6) = _sp(j,0);
          _buf(i,7) = _sp(j,1);
          _buf(i,8) = _sp(j,2);
          _buf(i,9) = _sp(j,3);
      }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecSpinKokkos::pack_border_kokkos(int n, DAT::tdual_int_1d k_sendlist, DAT::tdual_double_2d_lr buf,
                               int pbc_flag, int *pbc, ExecutionSpace space)
{
  double dx,dy,dz;

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
      AtomVecSpinKokkos_PackBorder<LMPHostType,1> f(
        buf.view_host(), k_sendlist.view_host(),
        h_x,h_tag,h_type,h_mask,h_sp,dx,dy,dz);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecSpinKokkos_PackBorder<LMPDeviceType,1> f(
        buf.view_device(), k_sendlist.view_device(),
        d_x,d_tag,d_type,d_mask,d_sp,dx,dy,dz);
      Kokkos::parallel_for(n,f);
    }

  } else {
    dx = dy = dz = 0;
    if(space==Host) {
      AtomVecSpinKokkos_PackBorder<LMPHostType,0> f(
        buf.view_host(), k_sendlist.view_host(),
        h_x,h_tag,h_type,h_mask,h_sp,dx,dy,dz);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecSpinKokkos_PackBorder<LMPDeviceType,0> f(
        buf.view_device(), k_sendlist.view_device(),
        d_x,d_tag,d_type,d_mask,d_sp,dx,dy,dz);
      Kokkos::parallel_for(n,f);
    }
  }
  return n*size_border;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecSpinKokkos_UnpackBorder {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  const typename AT::t_double_2d_lr_const _buf;
  typename AT::t_kkfloat_1d_3_lr _x;
  typename AT::t_tagint_1d _tag;
  typename AT::t_int_1d _type;
  typename AT::t_int_1d _mask;
  typename AT::t_kkfloat_1d_4 _sp;
  int _first;


  AtomVecSpinKokkos_UnpackBorder(
      const typename AT::t_double_2d_lr_const &buf,
      typename AT::t_kkfloat_1d_3_lr &x,
      typename AT::t_tagint_1d &tag,
      typename AT::t_int_1d &type,
      typename AT::t_int_1d &mask,
      typename AT::t_kkfloat_1d_4 &sp,
      const int& first):
    _buf(buf),_x(x),_tag(tag),_type(type),_mask(mask),_sp(sp),_first(first){
  };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
      _x(i+_first,0) = _buf(i,0);
      _x(i+_first,1) = _buf(i,1);
      _x(i+_first,2) = _buf(i,2);
      _tag(i+_first) = (tagint) d_ubuf(_buf(i,3)).i;
      _type(i+_first) = (int) d_ubuf(_buf(i,4)).i;
      _mask(i+_first) = (int) d_ubuf(_buf(i,5)).i;
      _sp(i+_first,0) = _buf(i,6);
      _sp(i+_first,1) = _buf(i,7);
      _sp(i+_first,2) = _buf(i,8);
      _sp(i+_first,3) = _buf(i,9);
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecSpinKokkos::unpack_border_kokkos(const int &n, const int &first,
                     const DAT::tdual_double_2d_lr &buf,ExecutionSpace space) {
  if (first+n >= nmax) {
    grow(first+n+100);
  }
  if(space==Host) {
    struct AtomVecSpinKokkos_UnpackBorder<LMPHostType>
      f(buf.view_host(),h_x,h_tag,h_type,h_mask,h_sp,first);
    Kokkos::parallel_for(n,f);
  } else {
    struct AtomVecSpinKokkos_UnpackBorder<LMPDeviceType>
      f(buf.view_device(),d_x,d_tag,d_type,d_mask,d_sp,first);
    Kokkos::parallel_for(n,f);
  }
  atomKK->modified(space,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|SP_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecSpinKokkos_PackExchangeFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_kkfloat_1d_3_lr_randomread _x;
  typename AT::t_kkfloat_1d_3_randomread _v;
  typename AT::t_tagint_1d_randomread _tag;
  typename AT::t_int_1d_randomread _type;
  typename AT::t_int_1d_randomread _mask;
  typename AT::t_imageint_1d_randomread _image;
  typename AT::t_kkfloat_1d_4_randomread _sp;
  typename AT::t_kkfloat_1d_3_lr _xw;
  typename AT::t_kkfloat_1d_3 _vw;
  typename AT::t_tagint_1d _tagw;
  typename AT::t_int_1d _typew;
  typename AT::t_int_1d _maskw;
  typename AT::t_imageint_1d _imagew;
  typename AT::t_kkfloat_1d_4 _spw;

  typename AT::t_double_2d_lr_um _buf;
  typename AT::t_int_1d_const _sendlist;
  typename AT::t_int_1d_const _copylist;
  int _size_exchange;

  AtomVecSpinKokkos_PackExchangeFunctor(
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
    _sp(atom->k_sp.view<DeviceType>()),
    _xw(atom->k_x.view<DeviceType>()),
    _vw(atom->k_v.view<DeviceType>()),
    _tagw(atom->k_tag.view<DeviceType>()),
    _typew(atom->k_type.view<DeviceType>()),
    _maskw(atom->k_mask.view<DeviceType>()),
    _imagew(atom->k_image.view<DeviceType>()),
    _spw(atom->k_sp.view<DeviceType>()),
    _sendlist(sendlist.template view<DeviceType>()),
    _copylist(copylist.template view<DeviceType>()),
    _size_exchange(atom->avecKK->size_exchange) {
    const int maxsendlist = (buf.template view<DeviceType>().extent(0)*
                             buf.template view<DeviceType>().extent(1))/_size_exchange;
    buffer_view<DeviceType>(_buf,buf,maxsendlist,_size_exchange);
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
    _buf(mysend,11) = _sp(i,0);
    _buf(mysend,12) = _sp(i,1);
    _buf(mysend,13) = _sp(i,2);
    _buf(mysend,14) = _sp(i,3);
    const int j = _copylist(mysend);

    if (j>-1) {
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
    _spw(i,0) = _sp(j,0);
    _spw(i,1) = _sp(j,1);
    _spw(i,2) = _sp(j,2);
    _spw(i,3) = _sp(j,3);
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecSpinKokkos::pack_exchange_kokkos(const int &nsend,DAT::tdual_double_2d_lr &k_buf,
                                              DAT::tdual_int_1d k_sendlist,
                                              DAT::tdual_int_1d k_copylist,
                                              ExecutionSpace space)
{
  size_exchange = 15;

  if (nsend > (int) (k_buf.view_host().extent(0)*k_buf.view_host().extent(1))/size_exchange) {
    int newsize = nsend*size_exchange/k_buf.view_host().extent(1)+1;
    k_buf.resize(newsize,k_buf.view_host().extent(1));
  }
  if (space == HostKK) {
    AtomVecSpinKokkos_PackExchangeFunctor<LMPHostType>
      f(atomKK,k_buf,k_sendlist,k_copylist);
    Kokkos::parallel_for(nsend,f);
    return nsend*size_exchange;
  } else {
    AtomVecSpinKokkos_PackExchangeFunctor<LMPDeviceType>
      f(atomKK,k_buf,k_sendlist,k_copylist);
    Kokkos::parallel_for(nsend,f);
    return nsend*size_exchange;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecSpinKokkos_UnpackExchangeFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_kkfloat_1d_3_lr _x;
  typename AT::t_kkfloat_1d_3 _v;
  typename AT::t_tagint_1d _tag;
  typename AT::t_int_1d _type;
  typename AT::t_int_1d _mask;
  typename AT::t_imageint_1d _image;
  typename AT::t_kkfloat_1d_4 _sp;
  typename AT::t_double_2d_lr_um _buf;
  typename AT::t_int_1d _nlocal;
  int _dim;
  double _lo,_hi;
  int _size_exchange;

  AtomVecSpinKokkos_UnpackExchangeFunctor(
      const AtomKokkos* atom,
      const DAT::tdual_double_2d_lr buf,
      DAT::tdual_int_1d nlocal,
      int dim, double lo, double hi):
      _x(atom->k_x.view<DeviceType>()),
      _v(atom->k_v.view<DeviceType>()),
      _tag(atom->k_tag.view<DeviceType>()),
      _type(atom->k_type.view<DeviceType>()),
      _mask(atom->k_mask.view<DeviceType>()),
      _image(atom->k_image.view<DeviceType>()),
      _sp(atom->k_sp.view<DeviceType>()),
      _nlocal(nlocal.template view<DeviceType>()),
      _dim(dim),_lo(lo),_hi(hi),_size_exchange(atom->avecKK->size_exchange) {
    const int maxsendlist = (buf.template view<DeviceType>().extent(0)*buf.template view<DeviceType>().extent(1))/_size_exchange;

    buffer_view<DeviceType>(_buf,buf,maxsendlist,_size_exchange);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &myrecv) const {
    double x = _buf(myrecv,_dim+1);
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
      _sp(i,0) = _buf(myrecv,11);
      _sp(i,1) = _buf(myrecv,12);
      _sp(i,2) = _buf(myrecv,13);
      _sp(i,3) = _buf(myrecv,14);
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecSpinKokkos::unpack_exchange_kokkos(DAT::tdual_double_2d_lr &k_buf, int nrecv, int nlocal,
                                              int dim, double lo, double hi, ExecutionSpace space,
                                              DAT::tdual_int_1d &/*k_indices*/)
{
  while (nlocal + nrecv/size_exchange >= nmax) grow(0);

  if(space == HostKK) {
    k_count.view_host()(0) = nlocal;
    AtomVecSpinKokkos_UnpackExchangeFunctor<LMPHostType> f(atomKK,k_buf,k_count,dim,lo,hi);
    Kokkos::parallel_for(nrecv/size_exchange,f);
    return k_count.view_host()(0);
  } else {
    k_count.view_host()(0) = nlocal;
    k_count.modify_host();
    k_count.sync_device();
    AtomVecSpinKokkos_UnpackExchangeFunctor<LMPDeviceType>
      f(atomKK,k_buf,k_count,dim,lo,hi);
    Kokkos::parallel_for(nrecv/size_exchange,f);
    k_count.modify_device();
    k_count.sync_host();

    return k_count.view_host()(0);
  }
}

/* ----------------------------------------------------------------------
   clear extra forces starting at atom N
   nbytes = # of bytes to clear for a per-atom vector
   include f b/c this is invoked from within SPIN pair styles
------------------------------------------------------------------------- */

void AtomVecSpinKokkos::force_clear(int /*n*/, size_t nbytes)
{
  int nzero = (double)nbytes/sizeof(double);

  if (nzero) {
    atomKK->k_fm.clear_sync_state(); // will be cleared below
    atomKK->k_fm_long.clear_sync_state(); // will be cleared below

    // local variables for lambda capture

    auto l_fm = atomKK->k_fm.view_device();
    auto l_fm_long = atomKK->k_fm_long.view_device();

    Kokkos::parallel_for(nzero, LAMMPS_LAMBDA(int i) {
      l_fm(i,0) = 0.0;
      l_fm(i,1) = 0.0;
      l_fm(i,2) = 0.0;
      l_fm_long(i,0) = 0.0;
      l_fm_long(i,1) = 0.0;
      l_fm_long(i,2) = 0.0;
    });

    atomKK->modified(Device,FM_MASK|FML_MASK);
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecSpinKokkos::sync(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.sync_device();
    if (mask & V_MASK) atomKK->k_v.sync_device();
    if (mask & F_MASK) atomKK->k_f.sync_device();
    if (mask & TAG_MASK) atomKK->k_tag.sync_device();
    if (mask & TYPE_MASK) atomKK->k_type.sync_device();
    if (mask & MASK_MASK) atomKK->k_mask.sync_device();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_device();
    if (mask & SP_MASK) atomKK->k_sp.sync_device();
    if (mask & FM_MASK) atomKK->k_fm.sync_device();
    if (mask & FML_MASK) atomKK->k_fm_long.sync_device();
  } else if (space == Host) {
    if (mask & X_MASK) atomKK->k_x.sync_host();
    if (mask & V_MASK) atomKK->k_v.sync_host();
    if (mask & F_MASK) atomKK->k_f.sync_host();
    if (mask & TAG_MASK) atomKK->k_tag.sync_host();
    if (mask & TYPE_MASK) atomKK->k_type.sync_host();
    if (mask & MASK_MASK) atomKK->k_mask.sync_host();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_host();
    if (mask & SP_MASK) atomKK->k_sp.sync_host();
    if (mask & FM_MASK) atomKK->k_fm.sync_host();
    if (mask & FML_MASK) atomKK->k_fm_long.sync_host();
  } else if (space == HostKK) {
    if (mask & X_MASK) atomKK->k_x.sync_hostkk();
    if (mask & V_MASK) atomKK->k_v.sync_hostkk();
    if (mask & F_MASK) atomKK->k_f.sync_hostkk();
    if (mask & TAG_MASK) atomKK->k_tag.sync_host();
    if (mask & TYPE_MASK) atomKK->k_type.sync_host();
    if (mask & MASK_MASK) atomKK->k_mask.sync_host();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_host();
    if (mask & SP_MASK) atomKK->k_sp.sync_hostkk();
    if (mask & FM_MASK) atomKK->k_fm.sync_hostkk();
    if (mask & FML_MASK) atomKK->k_fm_long.sync_hostkk();
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecSpinKokkos::modified(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.modify_device();
    if (mask & V_MASK) atomKK->k_v.modify_device();
    if (mask & F_MASK) atomKK->k_f.modify_device();
    if (mask & TAG_MASK) atomKK->k_tag.modify_device();
    if (mask & TYPE_MASK) atomKK->k_type.modify_device();
    if (mask & MASK_MASK) atomKK->k_mask.modify_device();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_device();
    if (mask & SP_MASK) atomKK->k_sp.modify_device();
    if (mask & FM_MASK) atomKK->k_fm.modify_device();
    if (mask & FML_MASK) atomKK->k_fm_long.modify_device();
  } else if (space == Host) {
    if (mask & X_MASK) atomKK->k_x.modify_host();
    if (mask & V_MASK) atomKK->k_v.modify_host();
    if (mask & F_MASK) atomKK->k_f.modify_host();
    if (mask & TAG_MASK) atomKK->k_tag.modify_host();
    if (mask & TYPE_MASK) atomKK->k_type.modify_host();
    if (mask & MASK_MASK) atomKK->k_mask.modify_host();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_host();
    if (mask & SP_MASK) atomKK->k_sp.modify_host();
    if (mask & FM_MASK) atomKK->k_fm.modify_host();
    if (mask & FML_MASK) atomKK->k_fm_long.modify_host();
  } else if (space == HostKK) {
    if (mask & X_MASK) atomKK->k_x.modify_hostkk();
    if (mask & V_MASK) atomKK->k_v.modify_hostkk();
    if (mask & F_MASK) atomKK->k_f.modify_hostkk();
    if (mask & TAG_MASK) atomKK->k_tag.modify_host();
    if (mask & TYPE_MASK) atomKK->k_type.modify_host();
    if (mask & MASK_MASK) atomKK->k_mask.modify_host();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_host();
    if (mask & SP_MASK) atomKK->k_sp.modify_hostkk();
    if (mask & FM_MASK) atomKK->k_fm.modify_hostkk();
    if (mask & FML_MASK) atomKK->k_fm_long.modify_hostkk();
  }
}

void AtomVecSpinKokkos::sync_pinned(ExecutionSpace space, unsigned int mask, int async_flag)
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
    if ((mask & SP_MASK) && atomKK->k_sp.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_4>(atomKK->k_sp,space,async_flag);
    if ((mask & FM_MASK) && atomKK->k_sp.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkacc_1d_3>(atomKK->k_fm,space,async_flag);
    if ((mask & FML_MASK) && atomKK->k_fm_long.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkacc_1d_3>(atomKK->k_fm_long,space,async_flag);
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
    if ((mask & SP_MASK) && atomKK->k_sp.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_4>(atomKK->k_sp,space,async_flag);
    if ((mask & FM_MASK) && atomKK->k_fm.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkacc_1d_3>(atomKK->k_fm,space,async_flag);
    if ((mask & FML_MASK) && atomKK->k_fm_long.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkacc_1d_3>(atomKK->k_fm_long,space,async_flag);
  }
}
