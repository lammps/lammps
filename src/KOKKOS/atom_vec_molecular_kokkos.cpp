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

#include "atom_vec_molecular_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm_kokkos.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "memory_kokkos.h"
#include "modify.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecMolecularKokkos::AtomVecMolecularKokkos(LAMMPS *lmp) : AtomVec(lmp),
AtomVecKokkos(lmp), AtomVecMolecular(lmp)
{

}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecMolecularKokkos::grow(int n)
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

  memoryKK->grow_kokkos(atomKK->k_num_dihedral,atomKK->num_dihedral,nmax,"atom:num_dihedral");
  memoryKK->grow_kokkos(atomKK->k_dihedral_type,atomKK->dihedral_type,nmax,
                      atomKK->dihedral_per_atom,"atom:dihedral_type");
  memoryKK->grow_kokkos(atomKK->k_dihedral_atom1,atomKK->dihedral_atom1,nmax,
                      atomKK->dihedral_per_atom,"atom:dihedral_atom1");
  memoryKK->grow_kokkos(atomKK->k_dihedral_atom2,atomKK->dihedral_atom2,nmax,
                      atomKK->dihedral_per_atom,"atom:dihedral_atom2");
  memoryKK->grow_kokkos(atomKK->k_dihedral_atom3,atomKK->dihedral_atom3,nmax,
                      atomKK->dihedral_per_atom,"atom:dihedral_atom3");
  memoryKK->grow_kokkos(atomKK->k_dihedral_atom4,atomKK->dihedral_atom4,nmax,
                      atomKK->dihedral_per_atom,"atom:dihedral_atom4");

  memoryKK->grow_kokkos(atomKK->k_num_improper,atomKK->num_improper,nmax,"atom:num_improper");
  memoryKK->grow_kokkos(atomKK->k_improper_type,atomKK->improper_type,nmax,
                      atomKK->improper_per_atom,"atom:improper_type");
  memoryKK->grow_kokkos(atomKK->k_improper_atom1,atomKK->improper_atom1,nmax,
                      atomKK->improper_per_atom,"atom:improper_atom1");
  memoryKK->grow_kokkos(atomKK->k_improper_atom2,atomKK->improper_atom2,nmax,
                      atomKK->improper_per_atom,"atom:improper_atom2");
  memoryKK->grow_kokkos(atomKK->k_improper_atom3,atomKK->improper_atom3,nmax,
                      atomKK->improper_per_atom,"atom:improper_atom3");
  memoryKK->grow_kokkos(atomKK->k_improper_atom4,atomKK->improper_atom4,nmax,
                      atomKK->improper_per_atom,"atom:improper_atom4");

  grow_pointers();
  atomKK->sync(Host,ALL_MASK);

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecMolecularKokkos::grow_pointers()
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

  num_dihedral = atomKK->num_dihedral;
  d_num_dihedral = atomKK->k_num_dihedral.d_view;
  h_num_dihedral = atomKK->k_num_dihedral.h_view;
  dihedral_type = atomKK->dihedral_type;
  d_dihedral_type = atomKK->k_dihedral_type.d_view;
  h_dihedral_type = atomKK->k_dihedral_type.h_view;
  dihedral_atom1 = atomKK->dihedral_atom1;
  d_dihedral_atom1 = atomKK->k_dihedral_atom1.d_view;
  h_dihedral_atom1 = atomKK->k_dihedral_atom1.h_view;
  dihedral_atom2 = atomKK->dihedral_atom2;
  d_dihedral_atom2 = atomKK->k_dihedral_atom2.d_view;
  h_dihedral_atom2 = atomKK->k_dihedral_atom2.h_view;
  dihedral_atom3 = atomKK->dihedral_atom3;
  d_dihedral_atom3 = atomKK->k_dihedral_atom3.d_view;
  h_dihedral_atom3 = atomKK->k_dihedral_atom3.h_view;
  dihedral_atom4 = atomKK->dihedral_atom4;
  d_dihedral_atom4 = atomKK->k_dihedral_atom4.d_view;
  h_dihedral_atom4 = atomKK->k_dihedral_atom4.h_view;

  num_improper = atomKK->num_improper;
  d_num_improper = atomKK->k_num_improper.d_view;
  h_num_improper = atomKK->k_num_improper.h_view;
  improper_type = atomKK->improper_type;
  d_improper_type = atomKK->k_improper_type.d_view;
  h_improper_type = atomKK->k_improper_type.h_view;
  improper_atom1 = atomKK->improper_atom1;
  d_improper_atom1 = atomKK->k_improper_atom1.d_view;
  h_improper_atom1 = atomKK->k_improper_atom1.h_view;
  improper_atom2 = atomKK->improper_atom2;
  d_improper_atom2 = atomKK->k_improper_atom2.d_view;
  h_improper_atom2 = atomKK->k_improper_atom2.h_view;
  improper_atom3 = atomKK->improper_atom3;
  d_improper_atom3 = atomKK->k_improper_atom3.d_view;
  h_improper_atom3 = atomKK->k_improper_atom3.h_view;
  improper_atom4 = atomKK->improper_atom4;
  d_improper_atom4 = atomKK->k_improper_atom4.d_view;
  h_improper_atom4 = atomKK->k_improper_atom4.h_view;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG,int TRICLINIC>
struct AtomVecMolecularKokkos_PackComm {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_x_array_randomread _x;
  typename ArrayTypes<DeviceType>::t_xfloat_2d_um _buf;
  typename ArrayTypes<DeviceType>::t_int_2d_const _list;
  const int _iswap;
  X_FLOAT _xprd,_yprd,_zprd,_xy,_xz,_yz;
  X_FLOAT _pbc[6];

  AtomVecMolecularKokkos_PackComm(
      const typename DAT::tdual_x_array &x,
      const typename DAT::tdual_xfloat_2d &buf,
      const typename DAT::tdual_int_2d &list,
      const int & iswap,
      const X_FLOAT &xprd, const X_FLOAT &yprd, const X_FLOAT &zprd,
      const X_FLOAT &xy, const X_FLOAT &xz, const X_FLOAT &yz, const int* const pbc):
      _x(x.view<DeviceType>()),_list(list.view<DeviceType>()),_iswap(iswap),
      _xprd(xprd),_yprd(yprd),_zprd(zprd),
      _xy(xy),_xz(xz),_yz(yz) {
        const size_t maxsend = (buf.view<DeviceType>().extent(0)
                                *buf.view<DeviceType>().extent(1))/3;
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

int AtomVecMolecularKokkos::pack_comm_kokkos(const int &n,
                                             const DAT::tdual_int_2d &list,
                                             const int & iswap,
                                             const DAT::tdual_xfloat_2d &buf,
                                             const int &pbc_flag,
                                             const int* const pbc)
{
  // Check whether to always run forward communication on the host
  // Choose correct forward PackComm kernel

  if (commKK->forward_comm_on_host) {
    atomKK->sync(Host,X_MASK);
    if (pbc_flag) {
      if (domain->triclinic) {
        struct AtomVecMolecularKokkos_PackComm<LMPHostType,1,1>
          f(atomKK->k_x,buf,list,iswap,domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecMolecularKokkos_PackComm<LMPHostType,1,0>
          f(atomKK->k_x,buf,list,iswap,domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      }
    } else {
      if (domain->triclinic) {
        struct AtomVecMolecularKokkos_PackComm<LMPHostType,0,1>
          f(atomKK->k_x,buf,list,iswap,domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecMolecularKokkos_PackComm<LMPHostType,0,0>
          f(atomKK->k_x,buf,list,iswap,domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      }
    }
  } else {
    atomKK->sync(Device,X_MASK);
    if (pbc_flag) {
      if (domain->triclinic) {
        struct AtomVecMolecularKokkos_PackComm<LMPDeviceType,1,1>
          f(atomKK->k_x,buf,list,iswap,domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecMolecularKokkos_PackComm<LMPDeviceType,1,0>
          f(atomKK->k_x,buf,list,iswap,domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      }
    } else {
      if (domain->triclinic) {
        struct AtomVecMolecularKokkos_PackComm<LMPDeviceType,0,1>
          f(atomKK->k_x,buf,list,iswap,domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecMolecularKokkos_PackComm<LMPDeviceType,0,0>
          f(atomKK->k_x,buf,list,iswap,domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      }
    }
  }

        return n*size_forward;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG,int TRICLINIC>
struct AtomVecMolecularKokkos_PackCommSelf {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_x_array_randomread _x;
  typename ArrayTypes<DeviceType>::t_x_array _xw;
  int _nfirst;
  typename ArrayTypes<DeviceType>::t_int_2d_const _list;
  const int _iswap;
  X_FLOAT _xprd,_yprd,_zprd,_xy,_xz,_yz;
  X_FLOAT _pbc[6];

  AtomVecMolecularKokkos_PackCommSelf(
      const typename DAT::tdual_x_array &x,
      const int &nfirst,
      const typename DAT::tdual_int_2d &list,
      const int & iswap,
      const X_FLOAT &xprd, const X_FLOAT &yprd, const X_FLOAT &zprd,
      const X_FLOAT &xy, const X_FLOAT &xz, const X_FLOAT &yz, const int* const pbc):
    _x(x.view<DeviceType>()),_xw(x.view<DeviceType>()),_nfirst(nfirst),
    _list(list.view<DeviceType>()),_iswap(iswap),
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

  }
};

/* ---------------------------------------------------------------------- */

int AtomVecMolecularKokkos::pack_comm_self(const int &n, const DAT::tdual_int_2d &list,
                                           const int & iswap,
                                           const int nfirst, const int &pbc_flag,
                                           const int* const pbc) {
  if (commKK->forward_comm_on_host) {
    atomKK->sync(Host,X_MASK);
    atomKK->modified(Host,X_MASK);
    if (pbc_flag) {
      if (domain->triclinic) {
      struct AtomVecMolecularKokkos_PackCommSelf<LMPHostType,1,1>
        f(atomKK->k_x,nfirst,list,iswap,domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      } else {
      struct AtomVecMolecularKokkos_PackCommSelf<LMPHostType,1,0>
        f(atomKK->k_x,nfirst,list,iswap,domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      }
    } else {
      if (domain->triclinic) {
      struct AtomVecMolecularKokkos_PackCommSelf<LMPHostType,0,1>
        f(atomKK->k_x,nfirst,list,iswap,domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      } else {
      struct AtomVecMolecularKokkos_PackCommSelf<LMPHostType,0,0>
        f(atomKK->k_x,nfirst,list,iswap,domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      }
    }
  } else {
    atomKK->sync(Device,X_MASK);
    atomKK->modified(Device,X_MASK);
    if (pbc_flag) {
      if (domain->triclinic) {
      struct AtomVecMolecularKokkos_PackCommSelf<LMPDeviceType,1,1>
        f(atomKK->k_x,nfirst,list,iswap,domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      } else {
      struct AtomVecMolecularKokkos_PackCommSelf<LMPDeviceType,1,0>
        f(atomKK->k_x,nfirst,list,iswap,domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      }
    } else {
      if (domain->triclinic) {
      struct AtomVecMolecularKokkos_PackCommSelf<LMPDeviceType,0,1>
        f(atomKK->k_x,nfirst,list,iswap,domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      } else {
      struct AtomVecMolecularKokkos_PackCommSelf<LMPDeviceType,0,0>
        f(atomKK->k_x,nfirst,list,iswap,domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      }
    }
  }
        return n*3;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecMolecularKokkos_UnpackComm {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_x_array _x;
  typename ArrayTypes<DeviceType>::t_xfloat_2d_const _buf;
  int _first;

  AtomVecMolecularKokkos_UnpackComm(
      const typename DAT::tdual_x_array &x,
      const typename DAT::tdual_xfloat_2d &buf,
      const int& first):_x(x.view<DeviceType>()),_buf(buf.view<DeviceType>()),
                        _first(first) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
      _x(i+_first,0) = _buf(i,0);
      _x(i+_first,1) = _buf(i,1);
      _x(i+_first,2) = _buf(i,2);
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecMolecularKokkos::unpack_comm_kokkos(const int &n, const int &first,
    const DAT::tdual_xfloat_2d &buf) {
  if (commKK->forward_comm_on_host) {
    atomKK->sync(Host,X_MASK);
    atomKK->modified(Host,X_MASK);
    struct AtomVecMolecularKokkos_UnpackComm<LMPHostType> f(atomKK->k_x,buf,first);
    Kokkos::parallel_for(n,f);
  } else {
    atomKK->sync(Device,X_MASK);
    atomKK->modified(Device,X_MASK);
    struct AtomVecMolecularKokkos_UnpackComm<LMPDeviceType> f(atomKK->k_x,buf,first);
    Kokkos::parallel_for(n,f);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG>
struct AtomVecMolecularKokkos_PackBorder {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_xfloat_2d _buf;
  const typename AT::t_int_2d_const _list;
  const int _iswap;
  const typename AT::t_x_array_randomread _x;
  const typename AT::t_tagint_1d _tag;
  const typename AT::t_int_1d _type;
  const typename AT::t_int_1d _mask;
  const typename AT::t_tagint_1d _molecule;
  X_FLOAT _dx,_dy,_dz;

  AtomVecMolecularKokkos_PackBorder(
      const typename AT::t_xfloat_2d &buf,
      const typename AT::t_int_2d_const &list,
      const int & iswap,
      const typename AT::t_x_array &x,
      const typename AT::t_tagint_1d &tag,
      const typename AT::t_int_1d &type,
      const typename AT::t_int_1d &mask,
      const typename AT::t_tagint_1d &molecule,
      const X_FLOAT &dx, const X_FLOAT &dy, const X_FLOAT &dz):
      _buf(buf),_list(list),_iswap(iswap),
      _x(x),_tag(tag),_type(type),_mask(mask),_molecule(molecule),
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

int AtomVecMolecularKokkos::pack_border_kokkos(int n, DAT::tdual_int_2d k_sendlist,
                                               DAT::tdual_xfloat_2d buf,int iswap,
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
      AtomVecMolecularKokkos_PackBorder<LMPHostType,1> f(
        buf.view<LMPHostType>(), k_sendlist.view<LMPHostType>(),
        iswap,h_x,h_tag,h_type,h_mask,h_molecule,dx,dy,dz);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecMolecularKokkos_PackBorder<LMPDeviceType,1> f(
        buf.view<LMPDeviceType>(), k_sendlist.view<LMPDeviceType>(),
        iswap,d_x,d_tag,d_type,d_mask,d_molecule,dx,dy,dz);
      Kokkos::parallel_for(n,f);
    }

  } else {
    dx = dy = dz = 0;
    if (space==Host) {
      AtomVecMolecularKokkos_PackBorder<LMPHostType,0> f(
        buf.view<LMPHostType>(), k_sendlist.view<LMPHostType>(),
        iswap,h_x,h_tag,h_type,h_mask,h_molecule,dx,dy,dz);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecMolecularKokkos_PackBorder<LMPDeviceType,0> f(
        buf.view<LMPDeviceType>(), k_sendlist.view<LMPDeviceType>(),
        iswap,d_x,d_tag,d_type,d_mask,d_molecule,dx,dy,dz);
      Kokkos::parallel_for(n,f);
    }
  }
  return n*size_border;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecMolecularKokkos_UnpackBorder {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  const typename AT::t_xfloat_2d_const _buf;
  typename AT::t_x_array _x;
  typename AT::t_tagint_1d _tag;
  typename AT::t_int_1d _type;
  typename AT::t_int_1d _mask;
  typename AT::t_tagint_1d _molecule;
  int _first;


  AtomVecMolecularKokkos_UnpackBorder(
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

void AtomVecMolecularKokkos::unpack_border_kokkos(const int &n, const int &first,
                                                  const DAT::tdual_xfloat_2d &buf,
                                                  ExecutionSpace space) {
  atomKK->modified(space,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|MOLECULE_MASK);
  while (first+n >= nmax) grow(0);
  atomKK->modified(space,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|MOLECULE_MASK);
  if (space==Host) {
    struct AtomVecMolecularKokkos_UnpackBorder<LMPHostType>
      f(buf.view<LMPHostType>(),h_x,h_tag,h_type,h_mask,h_molecule,first);
    Kokkos::parallel_for(n,f);
  } else {
    struct AtomVecMolecularKokkos_UnpackBorder<LMPDeviceType>
      f(buf.view<LMPDeviceType>(),d_x,d_tag,d_type,d_mask,d_molecule,first);
    Kokkos::parallel_for(n,f);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecMolecularKokkos_PackExchangeFunctor {
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
  typename AT::t_int_1d_randomread _num_dihedral;
  typename AT::t_int_2d_randomread _dihedral_type;
  typename AT::t_tagint_2d_randomread _dihedral_atom1,_dihedral_atom2,
    _dihedral_atom3,_dihedral_atom4;
  typename AT::t_int_1d_randomread _num_improper;
  typename AT::t_int_2d_randomread _improper_type;
  typename AT::t_tagint_2d_randomread _improper_atom1,_improper_atom2,
    _improper_atom3,_improper_atom4;
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
  typename AT::t_int_1d _num_dihedralw;
  typename AT::t_int_2d _dihedral_typew;
  typename AT::t_tagint_2d _dihedral_atom1w,_dihedral_atom2w,
    _dihedral_atom3w,_dihedral_atom4w;
  typename AT::t_int_1d _num_improperw;
  typename AT::t_int_2d _improper_typew;
  typename AT::t_tagint_2d _improper_atom1w,_improper_atom2w,
    _improper_atom3w,_improper_atom4w;
  typename AT::t_xfloat_2d_um _buf;
  typename AT::t_int_1d_const _sendlist;
  typename AT::t_int_1d_const _copylist;
  int _nlocal,_dim;
  X_FLOAT _lo,_hi;
  size_t elements;

  AtomVecMolecularKokkos_PackExchangeFunctor(
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
    _num_dihedral(atom->k_num_dihedral.view<DeviceType>()),
    _dihedral_type(atom->k_dihedral_type.view<DeviceType>()),
    _dihedral_atom1(atom->k_dihedral_atom1.view<DeviceType>()),
    _dihedral_atom2(atom->k_dihedral_atom2.view<DeviceType>()),
    _dihedral_atom3(atom->k_dihedral_atom3.view<DeviceType>()),
    _dihedral_atom4(atom->k_dihedral_atom4.view<DeviceType>()),
    _num_improper(atom->k_num_improper.view<DeviceType>()),
    _improper_type(atom->k_improper_type.view<DeviceType>()),
    _improper_atom1(atom->k_improper_atom1.view<DeviceType>()),
    _improper_atom2(atom->k_improper_atom2.view<DeviceType>()),
    _improper_atom3(atom->k_improper_atom3.view<DeviceType>()),
    _improper_atom4(atom->k_improper_atom4.view<DeviceType>()),
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
    _num_dihedralw(atom->k_num_dihedral.view<DeviceType>()),
    _dihedral_typew(atom->k_dihedral_type.view<DeviceType>()),
    _dihedral_atom1w(atom->k_dihedral_atom1.view<DeviceType>()),
    _dihedral_atom2w(atom->k_dihedral_atom2.view<DeviceType>()),
    _dihedral_atom3w(atom->k_dihedral_atom3.view<DeviceType>()),
    _dihedral_atom4w(atom->k_dihedral_atom4.view<DeviceType>()),
    _num_improperw(atom->k_num_improper.view<DeviceType>()),
    _improper_typew(atom->k_improper_type.view<DeviceType>()),
    _improper_atom1w(atom->k_improper_atom1.view<DeviceType>()),
    _improper_atom2w(atom->k_improper_atom2.view<DeviceType>()),
    _improper_atom3w(atom->k_improper_atom3.view<DeviceType>()),
    _improper_atom4w(atom->k_improper_atom4.view<DeviceType>()),
    _sendlist(sendlist.template view<DeviceType>()),
    _copylist(copylist.template view<DeviceType>()),
    _nlocal(nlocal),_dim(dim),
    _lo(lo),_hi(hi) {
    // 3 comp of x, 3 comp of v, 1 tag, 1 type, 1 mask, 1 image, 1 molecule, 3 nspecial,
    // maxspecial special, 1 num_bond, bond_per_atom bond_type, bond_per_atom bond_atom,
    // 1 num_angle, angle_per_atom angle_type, angle_per_atom angle_atom1, angle_atom2,
    // and angle_atom3
    // 1 num_dihedral, dihedral_per_atom dihedral_type, 4*dihedral_per_atom
    // 1 num_improper, 5*improper_per_atom
    // 1 to store buffer length
    elements = 19+atom->maxspecial+2*atom->bond_per_atom+4*atom->angle_per_atom+
      5*atom->dihedral_per_atom + 5*atom->improper_per_atom;
    const int maxsendlist = (buf.template view<DeviceType>().extent(0)*
                             buf.template view<DeviceType>().extent(1))/elements;
    buffer_view<DeviceType>(_buf,buf,maxsendlist,elements);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &mysend) const {
    int k;
    const int i = _sendlist(mysend);
    _buf(mysend,0) = elements;
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
    _buf(mysend,m++) = d_ubuf(_num_dihedral(i)).d;
    for (k = 0; k < _num_dihedral(i); k++) {
      _buf(mysend,m++) = d_ubuf(_dihedral_type(i,k)).d;
      _buf(mysend,m++) = d_ubuf(_dihedral_atom1(i,k)).d;
      _buf(mysend,m++) = d_ubuf(_dihedral_atom2(i,k)).d;
      _buf(mysend,m++) = d_ubuf(_dihedral_atom3(i,k)).d;
      _buf(mysend,m++) = d_ubuf(_dihedral_atom4(i,k)).d;
    }
    _buf(mysend,m++) = d_ubuf(_num_improper(i)).d;
    for (k = 0; k < _num_improper(i); k++) {
      _buf(mysend,m++) = d_ubuf(_improper_type(i,k)).d;
      _buf(mysend,m++) = d_ubuf(_improper_atom1(i,k)).d;
      _buf(mysend,m++) = d_ubuf(_improper_atom2(i,k)).d;
      _buf(mysend,m++) = d_ubuf(_improper_atom3(i,k)).d;
      _buf(mysend,m++) = d_ubuf(_improper_atom4(i,k)).d;
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
    _num_dihedralw(i) = _num_dihedral(j);
    for (k = 0; k < _num_dihedral(j); k++) {
      _dihedral_typew(i,k) = _dihedral_type(j,k);
      _dihedral_atom1w(i,k) = _dihedral_atom1(j,k);
      _dihedral_atom2w(i,k) = _dihedral_atom2(j,k);
      _dihedral_atom3w(i,k) = _dihedral_atom3(j,k);
      _dihedral_atom4w(i,k) = _dihedral_atom4(j,k);
    }
    _num_improperw(i) = _num_improper(j);
    for (k = 0; k < _num_improper(j); k++) {
      _improper_typew(i,k) = _improper_type(j,k);
      _improper_atom1w(i,k) = _improper_atom1(j,k);
      _improper_atom2w(i,k) = _improper_atom2(j,k);
      _improper_atom3w(i,k) = _improper_atom3(j,k);
      _improper_atom4w(i,k) = _improper_atom4(j,k);
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

int AtomVecMolecularKokkos::pack_exchange_kokkos(const int &nsend,DAT::tdual_xfloat_2d &k_buf,
                                                 DAT::tdual_int_1d k_sendlist,
                                                 DAT::tdual_int_1d k_copylist,
                                                 ExecutionSpace space,int dim,X_FLOAT lo,
                                                 X_FLOAT hi )
{
  const int elements = 19+atom->maxspecial+2*atom->bond_per_atom+4*atom->angle_per_atom+
      5*atom->dihedral_per_atom + 5*atom->improper_per_atom;
  if (nsend > (int) (k_buf.view<LMPHostType>().extent(0)*
              k_buf.view<LMPHostType>().extent(1))/elements) {
    int newsize = nsend*elements/k_buf.view<LMPHostType>().extent(1)+1;
    k_buf.resize(newsize,k_buf.view<LMPHostType>().extent(1));
  }
  if (space == Host) {
    AtomVecMolecularKokkos_PackExchangeFunctor<LMPHostType>
      f(atomKK,k_buf,k_sendlist,k_copylist,atom->nlocal,dim,lo,hi);
    Kokkos::parallel_for(nsend,f);
    return nsend*elements;
  } else {
    AtomVecMolecularKokkos_PackExchangeFunctor<LMPDeviceType>
      f(atomKK,k_buf,k_sendlist,k_copylist,atom->nlocal,dim,lo,hi);
    Kokkos::parallel_for(nsend,f);
    return nsend*elements;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecMolecularKokkos_UnpackExchangeFunctor {
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
  typename AT::t_int_1d _num_dihedral;
  typename AT::t_int_2d _dihedral_type;
  typename AT::t_tagint_2d _dihedral_atom1,_dihedral_atom2,
    _dihedral_atom3,_dihedral_atom4;
  typename AT::t_int_1d _num_improper;
  typename AT::t_int_2d _improper_type;
  typename AT::t_tagint_2d _improper_atom1,_improper_atom2,
    _improper_atom3,_improper_atom4;

  typename AT::t_xfloat_2d_um _buf;
  typename AT::t_int_1d _nlocal;
  int _dim;
  X_FLOAT _lo,_hi;
  size_t elements;

  AtomVecMolecularKokkos_UnpackExchangeFunctor(
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
    _num_dihedral(atom->k_num_dihedral.view<DeviceType>()),
    _dihedral_type(atom->k_dihedral_type.view<DeviceType>()),
    _dihedral_atom1(atom->k_dihedral_atom1.view<DeviceType>()),
    _dihedral_atom2(atom->k_dihedral_atom2.view<DeviceType>()),
    _dihedral_atom3(atom->k_dihedral_atom3.view<DeviceType>()),
    _dihedral_atom4(atom->k_dihedral_atom4.view<DeviceType>()),
    _num_improper(atom->k_num_improper.view<DeviceType>()),
    _improper_type(atom->k_improper_type.view<DeviceType>()),
    _improper_atom1(atom->k_improper_atom1.view<DeviceType>()),
    _improper_atom2(atom->k_improper_atom2.view<DeviceType>()),
    _improper_atom3(atom->k_improper_atom3.view<DeviceType>()),
    _improper_atom4(atom->k_improper_atom4.view<DeviceType>()),
    _nlocal(nlocal.template view<DeviceType>()),_dim(dim),
    _lo(lo),_hi(hi) {

    elements = 19+atom->maxspecial+2*atom->bond_per_atom+4*atom->angle_per_atom+
      5*atom->dihedral_per_atom + 5*atom->improper_per_atom;
    const int maxsendlist = (buf.template view<DeviceType>().extent(0)*
                             buf.template view<DeviceType>().extent(1))/elements;
    buffer_view<DeviceType>(_buf,buf,maxsendlist,elements);
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
        _bond_type(i,k) = (int) d_ubuf(_buf(myrecv,m++)).i;
        _bond_atom(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
      }
      _num_angle(i) =  (int) d_ubuf(_buf(myrecv,m++)).i;
      for (k = 0; k < _num_angle(i); k++) {
        _angle_type(i,k) = (int) d_ubuf(_buf(myrecv,m++)).i;
        _angle_atom1(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
        _angle_atom2(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
        _angle_atom3(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
      }
      _num_dihedral(i) = d_ubuf(_buf(myrecv,m++)).i;
      for (k = 0; k < _num_dihedral(i); k++) {
        _dihedral_type(i,k) = (int) d_ubuf(_buf(myrecv,m++)).i;
        _dihedral_atom1(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
        _dihedral_atom2(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
        _dihedral_atom3(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
        _dihedral_atom4(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
      }
      _num_improper(i) =  (int) d_ubuf(_buf(myrecv,m++)).i;
      for (k = 0; k < (int) _num_improper(i); k++) {
        _improper_type(i,k) = (int) d_ubuf(_buf(myrecv,m++)).i;
        _improper_atom1(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
        _improper_atom2(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
        _improper_atom3(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
        _improper_atom4(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
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

int AtomVecMolecularKokkos::unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf,int nrecv,
                                                   int nlocal,int dim,X_FLOAT lo,X_FLOAT hi,
                                                   ExecutionSpace space) {
  const size_t elements = 19+atom->maxspecial+2*atom->bond_per_atom+4*atom->angle_per_atom+
    5*atom->dihedral_per_atom + 5*atom->improper_per_atom;

  while (nlocal + nrecv/elements >= nmax) grow(0);

  if (space == Host) {
    k_count.h_view(0) = nlocal;
    AtomVecMolecularKokkos_UnpackExchangeFunctor<LMPHostType>
      f(atomKK,k_buf,k_count,dim,lo,hi);
    Kokkos::parallel_for(nrecv/elements,f);
    return k_count.h_view(0);
  } else {
    k_count.h_view(0) = nlocal;
    k_count.modify<LMPHostType>();
    k_count.sync<LMPDeviceType>();
    AtomVecMolecularKokkos_UnpackExchangeFunctor<LMPDeviceType>
      f(atomKK,k_buf,k_count,dim,lo,hi);
    Kokkos::parallel_for(nrecv/elements,f);
    k_count.modify<LMPDeviceType>();
    k_count.sync<LMPHostType>();

    return k_count.h_view(0);
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecMolecularKokkos::sync(ExecutionSpace space, unsigned int mask)
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
    if (mask & DIHEDRAL_MASK) {
      atomKK->k_num_dihedral.sync<LMPDeviceType>();
      atomKK->k_dihedral_type.sync<LMPDeviceType>();
      atomKK->k_dihedral_atom1.sync<LMPDeviceType>();
      atomKK->k_dihedral_atom2.sync<LMPDeviceType>();
      atomKK->k_dihedral_atom3.sync<LMPDeviceType>();
      atomKK->k_dihedral_atom4.sync<LMPDeviceType>();
    }
    if (mask & IMPROPER_MASK) {
      atomKK->k_num_improper.sync<LMPDeviceType>();
      atomKK->k_improper_type.sync<LMPDeviceType>();
      atomKK->k_improper_atom1.sync<LMPDeviceType>();
      atomKK->k_improper_atom2.sync<LMPDeviceType>();
      atomKK->k_improper_atom3.sync<LMPDeviceType>();
      atomKK->k_improper_atom4.sync<LMPDeviceType>();
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
    if (mask & DIHEDRAL_MASK) {
      atomKK->k_num_dihedral.sync<LMPHostType>();
      atomKK->k_dihedral_type.sync<LMPHostType>();
      atomKK->k_dihedral_atom1.sync<LMPHostType>();
      atomKK->k_dihedral_atom2.sync<LMPHostType>();
      atomKK->k_dihedral_atom3.sync<LMPHostType>();
      atomKK->k_dihedral_atom4.sync<LMPHostType>();
    }
    if (mask & IMPROPER_MASK) {
      atomKK->k_num_improper.sync<LMPHostType>();
      atomKK->k_improper_type.sync<LMPHostType>();
      atomKK->k_improper_atom1.sync<LMPHostType>();
      atomKK->k_improper_atom2.sync<LMPHostType>();
      atomKK->k_improper_atom3.sync<LMPHostType>();
      atomKK->k_improper_atom4.sync<LMPHostType>();
    }
  }
}

void AtomVecMolecularKokkos::sync_overlapping_device(ExecutionSpace space, unsigned int mask)
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
    if (mask & DIHEDRAL_MASK) {
      if (atomKK->k_num_dihedral.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_int_1d>(atomKK->k_num_dihedral,space);
      if (atomKK->k_dihedral_type.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_int_2d>(atomKK->k_dihedral_type,space);
      if (atomKK->k_dihedral_atom1.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_dihedral_atom1,space);
      if (atomKK->k_dihedral_atom2.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_dihedral_atom2,space);
      if (atomKK->k_dihedral_atom3.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_dihedral_atom3,space);
    }
    if (mask & IMPROPER_MASK) {
      if (atomKK->k_num_improper.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_int_1d>(atomKK->k_num_improper,space);
      if (atomKK->k_improper_type.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_int_2d>(atomKK->k_improper_type,space);
      if (atomKK->k_improper_atom1.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_improper_atom1,space);
      if (atomKK->k_improper_atom2.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_improper_atom2,space);
      if (atomKK->k_improper_atom3.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_improper_atom3,space);
      if (atomKK->k_improper_atom4.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_improper_atom4,space);
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
    if (mask & DIHEDRAL_MASK) {
      if (atomKK->k_num_dihedral.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_int_1d>(atomKK->k_num_dihedral,space);
      if (atomKK->k_dihedral_type.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_int_2d>(atomKK->k_dihedral_type,space);
      if (atomKK->k_dihedral_atom1.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_dihedral_atom1,space);
      if (atomKK->k_dihedral_atom2.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_dihedral_atom2,space);
      if (atomKK->k_dihedral_atom3.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_dihedral_atom3,space);
      if (atomKK->k_dihedral_atom4.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_dihedral_atom4,space);
    }
    if (mask & IMPROPER_MASK) {
      if (atomKK->k_num_improper.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_int_1d>(atomKK->k_num_improper,space);
      if (atomKK->k_improper_type.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_int_2d>(atomKK->k_improper_type,space);
      if (atomKK->k_improper_atom1.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_improper_atom1,space);
      if (atomKK->k_improper_atom2.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_improper_atom2,space);
      if (atomKK->k_improper_atom3.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_improper_atom3,space);
      if (atomKK->k_improper_atom4.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_improper_atom4,space);
    }
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecMolecularKokkos::modified(ExecutionSpace space, unsigned int mask)
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
    if (mask & DIHEDRAL_MASK) {
      atomKK->k_num_dihedral.modify<LMPDeviceType>();
      atomKK->k_dihedral_type.modify<LMPDeviceType>();
      atomKK->k_dihedral_atom1.modify<LMPDeviceType>();
      atomKK->k_dihedral_atom2.modify<LMPDeviceType>();
      atomKK->k_dihedral_atom3.modify<LMPDeviceType>();
      atomKK->k_dihedral_atom4.modify<LMPDeviceType>();
    }
    if (mask & IMPROPER_MASK) {
      atomKK->k_num_improper.modify<LMPDeviceType>();
      atomKK->k_improper_type.modify<LMPDeviceType>();
      atomKK->k_improper_atom1.modify<LMPDeviceType>();
      atomKK->k_improper_atom2.modify<LMPDeviceType>();
      atomKK->k_improper_atom3.modify<LMPDeviceType>();
      atomKK->k_improper_atom4.modify<LMPDeviceType>();
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
    if (mask & DIHEDRAL_MASK) {
      atomKK->k_num_dihedral.modify<LMPHostType>();
      atomKK->k_dihedral_type.modify<LMPHostType>();
      atomKK->k_dihedral_atom1.modify<LMPHostType>();
      atomKK->k_dihedral_atom2.modify<LMPHostType>();
      atomKK->k_dihedral_atom3.modify<LMPHostType>();
      atomKK->k_dihedral_atom4.modify<LMPHostType>();
    }
    if (mask & IMPROPER_MASK) {
      atomKK->k_num_improper.modify<LMPHostType>();
      atomKK->k_improper_type.modify<LMPHostType>();
      atomKK->k_improper_atom1.modify<LMPHostType>();
      atomKK->k_improper_atom2.modify<LMPHostType>();
      atomKK->k_improper_atom3.modify<LMPHostType>();
      atomKK->k_improper_atom4.modify<LMPHostType>();
    }
  }
}
