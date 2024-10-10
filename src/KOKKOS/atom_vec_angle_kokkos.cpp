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

#include "atom_vec_angle_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "memory_kokkos.h"
#include "modify.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecAngleKokkos::AtomVecAngleKokkos(LAMMPS *lmp) : AtomVec(lmp),
AtomVecKokkos(lmp), AtomVecAngle(lmp)
{

}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecAngleKokkos::grow(int n)
{
  auto DELTA = LMP_KOKKOS_AV_DELTA;
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

  memoryKK->grow_kokkos(atomKK->k_x,atomKK->x,nmax,"atom:x");
  memoryKK->grow_kokkos(atomKK->k_v,atomKK->v,nmax,"atom:v");
  memoryKK->grow_kokkos(atomKK->k_f,atomKK->f,nmax,"atom:f");

  memoryKK->grow_kokkos(atomKK->k_molecule,atomKK->molecule,nmax,"atom:molecule");
  memoryKK->grow_kokkos(atomKK->k_nspecial,atomKK->nspecial,nmax,3,"atom:nspecial");
  memoryKK->grow_kokkos(atomKK->k_special,atomKK->special,nmax,atomKK->maxspecial,
                      "atom:special");
  memoryKK->grow_kokkos(atomKK->k_num_bond,atomKK->num_bond,nmax,"atom:num_bond");
  memoryKK->grow_kokkos(atomKK->k_bond_type,atomKK->bond_type,nmax,atomKK->bond_per_atom,
                      "atom:bond_type");
  memoryKK->grow_kokkos(atomKK->k_bond_atom,atomKK->bond_atom,nmax,atomKK->bond_per_atom,
                      "atom:bond_atom");

  memoryKK->grow_kokkos(atomKK->k_num_angle,atomKK->num_angle,nmax,"atom:num_angle");
  memoryKK->grow_kokkos(atomKK->k_angle_type,atomKK->angle_type,nmax,atomKK->angle_per_atom,
                      "atom:angle_type");
  memoryKK->grow_kokkos(atomKK->k_angle_atom1,atomKK->angle_atom1,nmax,atomKK->angle_per_atom,
                      "atom:angle_atom1");
  memoryKK->grow_kokkos(atomKK->k_angle_atom2,atomKK->angle_atom2,nmax,atomKK->angle_per_atom,
                      "atom:angle_atom2");
  memoryKK->grow_kokkos(atomKK->k_angle_atom3,atomKK->angle_atom3,nmax,atomKK->angle_per_atom,
                      "atom:angle_atom3");

  grow_pointers();
  atomKK->sync(Host,ALL_MASK);

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecAngleKokkos::grow_pointers()
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

  molecule = atomKK->molecule;
  d_molecule = atomKK->k_molecule.d_view;
  h_molecule = atomKK->k_molecule.h_view;
  nspecial = atomKK->nspecial;
  d_nspecial = atomKK->k_nspecial.d_view;
  h_nspecial = atomKK->k_nspecial.h_view;
  special = atomKK->special;
  d_special = atomKK->k_special.d_view;
  h_special = atomKK->k_special.h_view;
  num_bond = atomKK->num_bond;
  d_num_bond = atomKK->k_num_bond.d_view;
  h_num_bond = atomKK->k_num_bond.h_view;
  bond_type = atomKK->bond_type;
  d_bond_type = atomKK->k_bond_type.d_view;
  h_bond_type = atomKK->k_bond_type.h_view;
  bond_atom = atomKK->bond_atom;
  d_bond_atom = atomKK->k_bond_atom.d_view;
  h_bond_atom = atomKK->k_bond_atom.h_view;

  num_angle = atomKK->num_angle;
  d_num_angle = atomKK->k_num_angle.d_view;
  h_num_angle = atomKK->k_num_angle.h_view;
  angle_type = atomKK->angle_type;
  d_angle_type = atomKK->k_angle_type.d_view;
  h_angle_type = atomKK->k_angle_type.h_view;
  angle_atom1 = atomKK->angle_atom1;
  d_angle_atom1 = atomKK->k_angle_atom1.d_view;
  h_angle_atom1 = atomKK->k_angle_atom1.h_view;
  angle_atom2 = atomKK->angle_atom2;
  d_angle_atom2 = atomKK->k_angle_atom2.d_view;
  h_angle_atom2 = atomKK->k_angle_atom2.h_view;
  angle_atom3 = atomKK->angle_atom3;
  d_angle_atom3 = atomKK->k_angle_atom3.d_view;
  h_angle_atom3 = atomKK->k_angle_atom3.h_view;
}

/* ----------------------------------------------------------------------
   sort atom arrays on device
------------------------------------------------------------------------- */

void AtomVecAngleKokkos::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  atomKK->sync(Device, ALL_MASK & ~F_MASK);

  Sorter.sort(LMPDeviceType(), d_tag);
  Sorter.sort(LMPDeviceType(), d_type);
  Sorter.sort(LMPDeviceType(), d_mask);
  Sorter.sort(LMPDeviceType(), d_image);
  Sorter.sort(LMPDeviceType(), d_x);
  Sorter.sort(LMPDeviceType(), d_v);
  Sorter.sort(LMPDeviceType(), d_molecule);
  Sorter.sort(LMPDeviceType(), d_num_bond);
  Sorter.sort(LMPDeviceType(), d_bond_type);
  Sorter.sort(LMPDeviceType(), d_bond_atom);
  Sorter.sort(LMPDeviceType(), d_nspecial);
  Sorter.sort(LMPDeviceType(), d_special);
  Sorter.sort(LMPDeviceType(), d_num_angle);
  Sorter.sort(LMPDeviceType(), d_angle_type);
  Sorter.sort(LMPDeviceType(), d_angle_atom1);
  Sorter.sort(LMPDeviceType(), d_angle_atom2);
  Sorter.sort(LMPDeviceType(), d_angle_atom3);

  atomKK->modified(Device, ALL_MASK & ~F_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG>
struct AtomVecAngleKokkos_PackBorder {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_xfloat_2d _buf;
  const typename AT::t_int_1d_const _list;
  const typename AT::t_x_array_randomread _x;
  const typename AT::t_tagint_1d _tag;
  const typename AT::t_int_1d _type;
  const typename AT::t_int_1d _mask;
  const typename AT::t_tagint_1d _molecule;
  X_FLOAT _dx,_dy,_dz;

  AtomVecAngleKokkos_PackBorder(
      const typename AT::t_xfloat_2d &buf,
      const typename AT::t_int_1d_const &list,
      const typename AT::t_x_array &x,
      const typename AT::t_tagint_1d &tag,
      const typename AT::t_int_1d &type,
      const typename AT::t_int_1d &mask,
      const typename AT::t_tagint_1d &molecule,
      const X_FLOAT &dx, const X_FLOAT &dy, const X_FLOAT &dz):
      _buf(buf),_list(list),
      _x(x),_tag(tag),_type(type),_mask(mask),_molecule(molecule),
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
          _buf(i,6) = d_ubuf(_molecule(j)).d;
      } else {
          _buf(i,0) = _x(j,0) + _dx;
          _buf(i,1) = _x(j,1) + _dy;
          _buf(i,2) = _x(j,2) + _dz;
          _buf(i,3) = d_ubuf(_tag(j)).d;
          _buf(i,4) = d_ubuf(_type(j)).d;
          _buf(i,5) = d_ubuf(_mask(j)).d;
          _buf(i,6) = d_ubuf(_molecule(j)).d;
      }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecAngleKokkos::pack_border_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                                          DAT::tdual_xfloat_2d buf,
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
    if (space==Host) {
      AtomVecAngleKokkos_PackBorder<LMPHostType,1> f(
        buf.view<LMPHostType>(), k_sendlist.view<LMPHostType>(),
        h_x,h_tag,h_type,h_mask,h_molecule,dx,dy,dz);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecAngleKokkos_PackBorder<LMPDeviceType,1> f(
        buf.view<LMPDeviceType>(), k_sendlist.view<LMPDeviceType>(),
        d_x,d_tag,d_type,d_mask,d_molecule,dx,dy,dz);
      Kokkos::parallel_for(n,f);
    }

  } else {
    dx = dy = dz = 0;
    if (space==Host) {
      AtomVecAngleKokkos_PackBorder<LMPHostType,0> f(
        buf.view<LMPHostType>(), k_sendlist.view<LMPHostType>(),
        h_x,h_tag,h_type,h_mask,h_molecule,dx,dy,dz);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecAngleKokkos_PackBorder<LMPDeviceType,0> f(
        buf.view<LMPDeviceType>(), k_sendlist.view<LMPDeviceType>(),
        d_x,d_tag,d_type,d_mask,d_molecule,dx,dy,dz);
      Kokkos::parallel_for(n,f);
    }
  }
  return n*size_border;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecAngleKokkos_UnpackBorder {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  const typename AT::t_xfloat_2d_const _buf;
  typename AT::t_x_array _x;
  typename AT::t_tagint_1d _tag;
  typename AT::t_int_1d _type;
  typename AT::t_int_1d _mask;
  typename AT::t_tagint_1d _molecule;
  int _first;


  AtomVecAngleKokkos_UnpackBorder(
      const typename AT::t_xfloat_2d_const &buf,
      typename AT::t_x_array &x,
      typename AT::t_tagint_1d &tag,
      typename AT::t_int_1d &type,
      typename AT::t_int_1d &mask,
      typename AT::t_tagint_1d &molecule,
      const int& first):
    _buf(buf),_x(x),_tag(tag),_type(type),_mask(mask),_molecule(molecule),
    _first(first) {
  };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
      _x(i+_first,0) = _buf(i,0);
      _x(i+_first,1) = _buf(i,1);
      _x(i+_first,2) = _buf(i,2);
      _tag(i+_first) = (tagint) d_ubuf(_buf(i,3)).i;
      _type(i+_first) = (int) d_ubuf(_buf(i,4)).i;
      _mask(i+_first) = (int) d_ubuf(_buf(i,5)).i;
      _molecule(i+_first) = (tagint) d_ubuf(_buf(i,6)).i;

  }
};

/* ---------------------------------------------------------------------- */

void AtomVecAngleKokkos::unpack_border_kokkos(const int &n, const int &first,
                                             const DAT::tdual_xfloat_2d &buf,
                                             ExecutionSpace space) {
  atomKK->modified(space,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|MOLECULE_MASK);
  while (first+n >= nmax) grow(0);
  atomKK->modified(space,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|MOLECULE_MASK);
  if (space==Host) {
    struct AtomVecAngleKokkos_UnpackBorder<LMPHostType>
      f(buf.view<LMPHostType>(),h_x,h_tag,h_type,h_mask,h_molecule,first);
    Kokkos::parallel_for(n,f);
  } else {
    struct AtomVecAngleKokkos_UnpackBorder<LMPDeviceType>
      f(buf.view<LMPDeviceType>(),d_x,d_tag,d_type,d_mask,d_molecule,first);
    Kokkos::parallel_for(n,f);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecAngleKokkos_PackExchangeFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_x_array_randomread _x;
  typename AT::t_v_array_randomread _v;
  typename AT::t_tagint_1d_randomread _tag;
  typename AT::t_int_1d_randomread _type;
  typename AT::t_int_1d_randomread _mask;
  typename AT::t_imageint_1d_randomread _image;
  typename AT::t_tagint_1d_randomread _molecule;
  typename AT::t_int_2d_randomread _nspecial;
  typename AT::t_tagint_2d_randomread _special;
  typename AT::t_int_1d_randomread _num_bond;
  typename AT::t_int_2d_randomread _bond_type;
  typename AT::t_tagint_2d_randomread _bond_atom;
  typename AT::t_int_1d_randomread _num_angle;
  typename AT::t_int_2d_randomread _angle_type;
  typename AT::t_tagint_2d_randomread _angle_atom1,_angle_atom2,_angle_atom3;
  typename AT::t_x_array _xw;
  typename AT::t_v_array _vw;
  typename AT::t_tagint_1d _tagw;
  typename AT::t_int_1d _typew;
  typename AT::t_int_1d _maskw;
  typename AT::t_imageint_1d _imagew;
  typename AT::t_tagint_1d _moleculew;
  typename AT::t_int_2d _nspecialw;
  typename AT::t_tagint_2d _specialw;
  typename AT::t_int_1d _num_bondw;
  typename AT::t_int_2d _bond_typew;
  typename AT::t_tagint_2d _bond_atomw;
  typename AT::t_int_1d _num_anglew;
  typename AT::t_int_2d _angle_typew;
  typename AT::t_tagint_2d _angle_atom1w,_angle_atom2w,_angle_atom3w;

  typename AT::t_xfloat_2d_um _buf;
  typename AT::t_int_1d_const _sendlist;
  typename AT::t_int_1d_const _copylist;
  int _size_exchange;

  AtomVecAngleKokkos_PackExchangeFunctor(
      const AtomKokkos* atom,
      const typename AT::tdual_xfloat_2d buf,
      typename AT::tdual_int_1d sendlist,
      typename AT::tdual_int_1d copylist):
    _x(atom->k_x.view<DeviceType>()),
    _v(atom->k_v.view<DeviceType>()),
    _tag(atom->k_tag.view<DeviceType>()),
    _type(atom->k_type.view<DeviceType>()),
    _mask(atom->k_mask.view<DeviceType>()),
    _image(atom->k_image.view<DeviceType>()),
    _molecule(atom->k_molecule.view<DeviceType>()),
    _nspecial(atom->k_nspecial.view<DeviceType>()),
    _special(atom->k_special.view<DeviceType>()),
    _num_bond(atom->k_num_bond.view<DeviceType>()),
    _bond_type(atom->k_bond_type.view<DeviceType>()),
    _bond_atom(atom->k_bond_atom.view<DeviceType>()),
    _num_angle(atom->k_num_angle.view<DeviceType>()),
    _angle_type(atom->k_angle_type.view<DeviceType>()),
    _angle_atom1(atom->k_angle_atom1.view<DeviceType>()),
    _angle_atom2(atom->k_angle_atom2.view<DeviceType>()),
    _angle_atom3(atom->k_angle_atom3.view<DeviceType>()),
    _xw(atom->k_x.view<DeviceType>()),
    _vw(atom->k_v.view<DeviceType>()),
    _tagw(atom->k_tag.view<DeviceType>()),
    _typew(atom->k_type.view<DeviceType>()),
    _maskw(atom->k_mask.view<DeviceType>()),
    _imagew(atom->k_image.view<DeviceType>()),
    _moleculew(atom->k_molecule.view<DeviceType>()),
    _nspecialw(atom->k_nspecial.view<DeviceType>()),
    _specialw(atom->k_special.view<DeviceType>()),
    _num_bondw(atom->k_num_bond.view<DeviceType>()),
    _bond_typew(atom->k_bond_type.view<DeviceType>()),
    _bond_atomw(atom->k_bond_atom.view<DeviceType>()),
    _num_anglew(atom->k_num_angle.view<DeviceType>()),
    _angle_typew(atom->k_angle_type.view<DeviceType>()),
    _angle_atom1w(atom->k_angle_atom1.view<DeviceType>()),
    _angle_atom2w(atom->k_angle_atom2.view<DeviceType>()),
    _angle_atom3w(atom->k_angle_atom3.view<DeviceType>()),
    _sendlist(sendlist.template view<DeviceType>()),
    _copylist(copylist.template view<DeviceType>()),
    _size_exchange(atom->avecKK->size_exchange) {
    const int maxsendlist = (buf.template view<DeviceType>().extent(0)*
                             buf.template view<DeviceType>().extent(1))/_size_exchange;
    buffer_view<DeviceType>(_buf,buf,maxsendlist,_size_exchange);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &mysend) const {
    int k;
    const int i = _sendlist(mysend);
    _buf(mysend,0) = _size_exchange;
    int m = 1;
    _buf(mysend,m++) = _x(i,0);
    _buf(mysend,m++) = _x(i,1);
    _buf(mysend,m++) = _x(i,2);
    _buf(mysend,m++) = _v(i,0);
    _buf(mysend,m++) = _v(i,1);
    _buf(mysend,m++) = _v(i,2);
    _buf(mysend,m++) = d_ubuf(_tag(i)).d;
    _buf(mysend,m++) = d_ubuf(_type(i)).d;
    _buf(mysend,m++) = d_ubuf(_mask(i)).d;
    _buf(mysend,m++) = d_ubuf(_image(i)).d;
    _buf(mysend,m++) = d_ubuf(_molecule(i)).d;
    _buf(mysend,m++) = d_ubuf(_num_bond(i)).d;
    for (k = 0; k < _num_bond(i); k++) {
      _buf(mysend,m++) = d_ubuf(_bond_type(i,k)).d;
      _buf(mysend,m++) = d_ubuf(_bond_atom(i,k)).d;
    }
    _buf(mysend,m++) = d_ubuf(_num_angle(i)).d;
    for (k = 0; k < _num_angle(i); k++) {
      _buf(mysend,m++) = d_ubuf(_angle_type(i,k)).d;
      _buf(mysend,m++) = d_ubuf(_angle_atom1(i,k)).d;
      _buf(mysend,m++) = d_ubuf(_angle_atom2(i,k)).d;
      _buf(mysend,m++) = d_ubuf(_angle_atom3(i,k)).d;
    }
    _buf(mysend,m++) = d_ubuf(_nspecial(i,0)).d;
    _buf(mysend,m++) = d_ubuf(_nspecial(i,1)).d;
    _buf(mysend,m++) = d_ubuf(_nspecial(i,2)).d;
    for (k = 0; k < _nspecial(i,2); k++)
      _buf(mysend,m++) = d_ubuf(_special(i,k)).d;

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
    _moleculew(i) = _molecule(j);
    _num_bondw(i) = _num_bond(j);
    for (k = 0; k < _num_bond(j); k++) {
      _bond_typew(i,k) = _bond_type(j,k);
      _bond_atomw(i,k) = _bond_atom(j,k);
    }
    _num_anglew(i) = _num_angle(j);
    for (k = 0; k < _num_angle(j); k++) {
      _angle_typew(i,k) = _angle_type(j,k);
      _angle_atom1w(i,k) = _angle_atom1(j,k);
      _angle_atom2w(i,k) = _angle_atom2(j,k);
      _angle_atom3w(i,k) = _angle_atom3(j,k);
    }
    _nspecialw(i,0) = _nspecial(j,0);
    _nspecialw(i,1) = _nspecial(j,1);
    _nspecialw(i,2) = _nspecial(j,2);
    for (k = 0; k < _nspecial(j,2); k++)
      _specialw(i,k) = _special(j,k);
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecAngleKokkos::pack_exchange_kokkos(const int &nsend,DAT::tdual_xfloat_2d &k_buf,
                                            DAT::tdual_int_1d k_sendlist,
                                            DAT::tdual_int_1d k_copylist,
                                            ExecutionSpace space)
{

  // maxspecial special, 1 num_bond, bond_per_atom bond_type, bond_per_atom bond_atom,
  // 1 num_angle, angle_per_atom angle_type, angle_per_atom angle_atom1, angle_atom2,
  // and angle_atom3
  // 1 to store buffer length

  size_exchange = 17+atom->maxspecial+2*atom->bond_per_atom+4*atom->angle_per_atom;

  if (nsend > (int) (k_buf.view<LMPHostType>().extent(0)*
              k_buf.view<LMPHostType>().extent(1))/size_exchange) {
    int newsize = nsend*size_exchange/k_buf.view<LMPHostType>().extent(1)+1;
    k_buf.resize(newsize,k_buf.view<LMPHostType>().extent(1));
  }
  if (space == Host) {
    AtomVecAngleKokkos_PackExchangeFunctor<LMPHostType>
      f(atomKK,k_buf,k_sendlist,k_copylist);
    Kokkos::parallel_for(nsend,f);
    return nsend*size_exchange;
  } else {
    AtomVecAngleKokkos_PackExchangeFunctor<LMPDeviceType>
      f(atomKK,k_buf,k_sendlist,k_copylist);
    Kokkos::parallel_for(nsend,f);
    return nsend*size_exchange;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecAngleKokkos_UnpackExchangeFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_x_array _x;
  typename AT::t_v_array _v;
  typename AT::t_tagint_1d _tag;
  typename AT::t_int_1d _type;
  typename AT::t_int_1d _mask;
  typename AT::t_imageint_1d _image;
  typename AT::t_tagint_1d _molecule;
  typename AT::t_int_2d _nspecial;
  typename AT::t_tagint_2d _special;
  typename AT::t_int_1d _num_bond;
  typename AT::t_int_2d _bond_type;
  typename AT::t_tagint_2d _bond_atom;
  typename AT::t_int_1d _num_angle;
  typename AT::t_int_2d _angle_type;
  typename AT::t_tagint_2d _angle_atom1,_angle_atom2,_angle_atom3;

  typename AT::t_xfloat_2d_um _buf;
  typename AT::t_int_1d _nlocal;
  int _dim;
  X_FLOAT _lo,_hi;
  int _size_exchange;

  AtomVecAngleKokkos_UnpackExchangeFunctor(
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
    _molecule(atom->k_molecule.view<DeviceType>()),
    _nspecial(atom->k_nspecial.view<DeviceType>()),
    _special(atom->k_special.view<DeviceType>()),
    _num_bond(atom->k_num_bond.view<DeviceType>()),
    _bond_type(atom->k_bond_type.view<DeviceType>()),
    _bond_atom(atom->k_bond_atom.view<DeviceType>()),
    _num_angle(atom->k_num_angle.view<DeviceType>()),
    _angle_type(atom->k_angle_type.view<DeviceType>()),
    _angle_atom1(atom->k_angle_atom1.view<DeviceType>()),
    _angle_atom2(atom->k_angle_atom2.view<DeviceType>()),
    _angle_atom3(atom->k_angle_atom3.view<DeviceType>()),
    _nlocal(nlocal.template view<DeviceType>()),
    _dim(dim),_lo(lo),_hi(hi),_size_exchange(atom->avecKK->size_exchange) {
    const int maxsendlist = (buf.template view<DeviceType>().extent(0)*
                             buf.template view<DeviceType>().extent(1))/_size_exchange;
    buffer_view<DeviceType>(_buf,buf,maxsendlist,_size_exchange);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &myrecv) const {
    X_FLOAT x = _buf(myrecv,_dim+1);
    if (x >= _lo && x < _hi) {
      int i = Kokkos::atomic_fetch_add(&_nlocal(0),1);
      int m = 1;
      _x(i,0) = _buf(myrecv,m++);
      _x(i,1) = _buf(myrecv,m++);
      _x(i,2) = _buf(myrecv,m++);
      _v(i,0) = _buf(myrecv,m++);
      _v(i,1) = _buf(myrecv,m++);
      _v(i,2) = _buf(myrecv,m++);
      _tag(i) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
      _type(i) = (int) d_ubuf(_buf(myrecv,m++)).i;
      _mask(i) = (int) d_ubuf(_buf(myrecv,m++)).i;
      _image(i) = (imageint) d_ubuf(_buf(myrecv,m++)).i;

      _molecule(i) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
      _num_bond(i) = (int) d_ubuf(_buf(myrecv,m++)).i;
      int k;
      for (k = 0; k < _num_bond(i); k++) {
        _bond_type(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
        _bond_atom(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
      }
      _num_angle(i) = (int) d_ubuf(_buf(myrecv,m++)).i;
      for (k = 0; k < _num_angle(i); k++) {
        _angle_type(i,k) = (int) d_ubuf(_buf(myrecv,m++)).i;
        _angle_atom1(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
        _angle_atom2(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
        _angle_atom3(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
      }
      _nspecial(i,0) = (int) d_ubuf(_buf(myrecv,m++)).i;
      _nspecial(i,1) = (int) d_ubuf(_buf(myrecv,m++)).i;
      _nspecial(i,2) = (int) d_ubuf(_buf(myrecv,m++)).i;
      for (k = 0; k < _nspecial(i,2); k++)
        _special(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecAngleKokkos::unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf, int nrecv, int nlocal,
                                               int dim, X_FLOAT lo, X_FLOAT hi, ExecutionSpace space,
                                               DAT::tdual_int_1d &/*k_indices*/)
{
  while (nlocal + nrecv/size_exchange >= nmax) grow(0);

  if (space == Host) {
    k_count.h_view(0) = nlocal;
    AtomVecAngleKokkos_UnpackExchangeFunctor<LMPHostType>
      f(atomKK,k_buf,k_count,dim,lo,hi);
    Kokkos::parallel_for(nrecv/size_exchange,f);
    return k_count.h_view(0);
  } else {
    k_count.h_view(0) = nlocal;
    k_count.modify<LMPHostType>();
    k_count.sync<LMPDeviceType>();
    AtomVecAngleKokkos_UnpackExchangeFunctor<LMPDeviceType>
      f(atomKK,k_buf,k_count,dim,lo,hi);
    Kokkos::parallel_for(nrecv/size_exchange,f);
    k_count.modify<LMPDeviceType>();
    k_count.sync<LMPHostType>();

    return k_count.h_view(0);
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecAngleKokkos::sync(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.sync<LMPDeviceType>();
    if (mask & V_MASK) atomKK->k_v.sync<LMPDeviceType>();
    if (mask & F_MASK) atomKK->k_f.sync<LMPDeviceType>();
    if (mask & TAG_MASK) atomKK->k_tag.sync<LMPDeviceType>();
    if (mask & TYPE_MASK) atomKK->k_type.sync<LMPDeviceType>();
    if (mask & MASK_MASK) atomKK->k_mask.sync<LMPDeviceType>();
    if (mask & IMAGE_MASK) atomKK->k_image.sync<LMPDeviceType>();
    if (mask & MOLECULE_MASK) atomKK->k_molecule.sync<LMPDeviceType>();
    if (mask & SPECIAL_MASK) {
      atomKK->k_nspecial.sync<LMPDeviceType>();
      atomKK->k_special.sync<LMPDeviceType>();
    }
    if (mask & BOND_MASK) {
      atomKK->k_num_bond.sync<LMPDeviceType>();
      atomKK->k_bond_type.sync<LMPDeviceType>();
      atomKK->k_bond_atom.sync<LMPDeviceType>();
    }
    if (mask & ANGLE_MASK) {
      atomKK->k_num_angle.sync<LMPDeviceType>();
      atomKK->k_angle_type.sync<LMPDeviceType>();
      atomKK->k_angle_atom1.sync<LMPDeviceType>();
      atomKK->k_angle_atom2.sync<LMPDeviceType>();
      atomKK->k_angle_atom3.sync<LMPDeviceType>();
    }
  } else {
    if (mask & X_MASK) atomKK->k_x.sync<LMPHostType>();
    if (mask & V_MASK) atomKK->k_v.sync<LMPHostType>();
    if (mask & F_MASK) atomKK->k_f.sync<LMPHostType>();
    if (mask & TAG_MASK) atomKK->k_tag.sync<LMPHostType>();
    if (mask & TYPE_MASK) atomKK->k_type.sync<LMPHostType>();
    if (mask & MASK_MASK) atomKK->k_mask.sync<LMPHostType>();
    if (mask & IMAGE_MASK) atomKK->k_image.sync<LMPHostType>();
    if (mask & MOLECULE_MASK) atomKK->k_molecule.sync<LMPHostType>();
    if (mask & SPECIAL_MASK) {
      atomKK->k_nspecial.sync<LMPHostType>();
      atomKK->k_special.sync<LMPHostType>();
    }
    if (mask & BOND_MASK) {
      atomKK->k_num_bond.sync<LMPHostType>();
      atomKK->k_bond_type.sync<LMPHostType>();
      atomKK->k_bond_atom.sync<LMPHostType>();
    }
    if (mask & ANGLE_MASK) {
      atomKK->k_num_angle.sync<LMPHostType>();
      atomKK->k_angle_type.sync<LMPHostType>();
      atomKK->k_angle_atom1.sync<LMPHostType>();
      atomKK->k_angle_atom2.sync<LMPHostType>();
      atomKK->k_angle_atom3.sync<LMPHostType>();
    }
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecAngleKokkos::sync_overlapping_device(ExecutionSpace space, unsigned int mask)
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
    if ((mask & MOLECULE_MASK) && atomKK->k_molecule.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_tagint_1d>(atomKK->k_molecule,space);
    if (mask & SPECIAL_MASK) {
      if (atomKK->k_nspecial.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_int_2d>(atomKK->k_nspecial,space);
      if (atomKK->k_special.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_special,space);
    }
    if (mask & BOND_MASK) {
      if (atomKK->k_num_bond.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_int_1d>(atomKK->k_num_bond,space);
      if (atomKK->k_bond_type.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_int_2d>(atomKK->k_bond_type,space);
      if (atomKK->k_bond_atom.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_bond_atom,space);
    }
    if (mask & ANGLE_MASK) {
      if (atomKK->k_num_angle.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_int_1d>(atomKK->k_num_angle,space);
      if (atomKK->k_angle_type.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_int_2d>(atomKK->k_angle_type,space);
      if (atomKK->k_angle_atom1.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_angle_atom1,space);
      if (atomKK->k_angle_atom2.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_angle_atom2,space);
      if (atomKK->k_angle_atom3.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_angle_atom3,space);
    }
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
    if ((mask & MOLECULE_MASK) && atomKK->k_molecule.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_tagint_1d>(atomKK->k_molecule,space);
    if (mask & SPECIAL_MASK) {
      if (atomKK->k_nspecial.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_int_2d>(atomKK->k_nspecial,space);
      if (atomKK->k_special.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_special,space);
    }
    if (mask & BOND_MASK) {
      if (atomKK->k_num_bond.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_int_1d>(atomKK->k_num_bond,space);
      if (atomKK->k_bond_type.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_int_2d>(atomKK->k_bond_type,space);
      if (atomKK->k_bond_atom.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_bond_atom,space);
    }
    if (mask & ANGLE_MASK) {
      if (atomKK->k_num_angle.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_int_1d>(atomKK->k_num_angle,space);
      if (atomKK->k_angle_type.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_int_2d>(atomKK->k_angle_type,space);
      if (atomKK->k_angle_atom1.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_angle_atom1,space);
      if (atomKK->k_angle_atom2.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_angle_atom2,space);
      if (atomKK->k_angle_atom3.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_angle_atom3,space);
    }
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecAngleKokkos::modified(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.modify<LMPDeviceType>();
    if (mask & V_MASK) atomKK->k_v.modify<LMPDeviceType>();
    if (mask & F_MASK) atomKK->k_f.modify<LMPDeviceType>();
    if (mask & TAG_MASK) atomKK->k_tag.modify<LMPDeviceType>();
    if (mask & TYPE_MASK) atomKK->k_type.modify<LMPDeviceType>();
    if (mask & MASK_MASK) atomKK->k_mask.modify<LMPDeviceType>();
    if (mask & IMAGE_MASK) atomKK->k_image.modify<LMPDeviceType>();
    if (mask & MOLECULE_MASK) atomKK->k_molecule.modify<LMPDeviceType>();
    if (mask & SPECIAL_MASK) {
      atomKK->k_nspecial.modify<LMPDeviceType>();
      atomKK->k_special.modify<LMPDeviceType>();
    }
    if (mask & BOND_MASK) {
      atomKK->k_num_bond.modify<LMPDeviceType>();
      atomKK->k_bond_type.modify<LMPDeviceType>();
      atomKK->k_bond_atom.modify<LMPDeviceType>();
    }
    if (mask & ANGLE_MASK) {
      atomKK->k_num_angle.modify<LMPDeviceType>();
      atomKK->k_angle_type.modify<LMPDeviceType>();
      atomKK->k_angle_atom1.modify<LMPDeviceType>();
      atomKK->k_angle_atom2.modify<LMPDeviceType>();
      atomKK->k_angle_atom3.modify<LMPDeviceType>();
    }
  } else {
    if (mask & X_MASK) atomKK->k_x.modify<LMPHostType>();
    if (mask & V_MASK) atomKK->k_v.modify<LMPHostType>();
    if (mask & F_MASK) atomKK->k_f.modify<LMPHostType>();
    if (mask & TAG_MASK) atomKK->k_tag.modify<LMPHostType>();
    if (mask & TYPE_MASK) atomKK->k_type.modify<LMPHostType>();
    if (mask & MASK_MASK) atomKK->k_mask.modify<LMPHostType>();
    if (mask & IMAGE_MASK) atomKK->k_image.modify<LMPHostType>();
    if (mask & MOLECULE_MASK) atomKK->k_molecule.modify<LMPHostType>();
    if (mask & SPECIAL_MASK) {
      atomKK->k_nspecial.modify<LMPHostType>();
      atomKK->k_special.modify<LMPHostType>();
    }
    if (mask & BOND_MASK) {
      atomKK->k_num_bond.modify<LMPHostType>();
      atomKK->k_bond_type.modify<LMPHostType>();
      atomKK->k_bond_atom.modify<LMPHostType>();
    }
    if (mask & ANGLE_MASK) {
      atomKK->k_num_angle.modify<LMPHostType>();
      atomKK->k_angle_type.modify<LMPHostType>();
      atomKK->k_angle_atom1.modify<LMPHostType>();
      atomKK->k_angle_atom2.modify<LMPHostType>();
      atomKK->k_angle_atom3.modify<LMPHostType>();
    }
  }
}
