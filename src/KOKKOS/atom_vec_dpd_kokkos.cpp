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

#include "atom_vec_dpd_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm_kokkos.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "modify.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecDPDKokkos::AtomVecDPDKokkos(LAMMPS *lmp) : AtomVec(lmp),
AtomVecKokkos(lmp), AtomVecDPD(lmp)
{
  no_comm_vel_flag = 1;
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecDPDKokkos::grow(int n)
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


  memoryKK->grow_kokkos(atomKK->k_rho,atomKK->rho,nmax,"atom:rho");
  memoryKK->grow_kokkos(atomKK->k_dpdTheta,atomKK->dpdTheta,nmax,"atom:dpdTheta");
  memoryKK->grow_kokkos(atomKK->k_uCond,atomKK->uCond,nmax,"atom:uCond");
  memoryKK->grow_kokkos(atomKK->k_uMech,atomKK->uMech,nmax,"atom:uMech");
  memoryKK->grow_kokkos(atomKK->k_uChem,atomKK->uChem,nmax,"atom:uChem");
  memoryKK->grow_kokkos(atomKK->k_uCG,atomKK->uCG,nmax,"atom:uCG");
  memoryKK->grow_kokkos(atomKK->k_uCGnew,atomKK->uCGnew,nmax,"atom:uCGnew");
  memoryKK->grow_kokkos(atomKK->k_duChem,atomKK->duChem,nmax,"atom:duChem");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);

  grow_pointers();
  atomKK->sync(Host,ALL_MASK);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecDPDKokkos::grow_pointers()
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

  rho = atomKK->rho;
  d_rho = atomKK->k_rho.d_view;
  h_rho = atomKK->k_rho.h_view;
  dpdTheta = atomKK->dpdTheta;
  d_dpdTheta = atomKK->k_dpdTheta.d_view;
  h_dpdTheta = atomKK->k_dpdTheta.h_view;
  uCond = atomKK->uCond;
  d_uCond = atomKK->k_uCond.d_view;
  h_uCond = atomKK->k_uCond.h_view;
  uMech = atomKK->uMech;
  d_uMech = atomKK->k_uMech.d_view;
  h_uMech = atomKK->k_uMech.h_view;
  uChem = atomKK->uChem;
  d_uChem = atomKK->k_uChem.d_view;
  h_uChem = atomKK->k_uChem.h_view;
  uCG = atomKK->uCG;
  d_uCG = atomKK->k_uCG.d_view;
  h_uCG = atomKK->k_uCG.h_view;
  uCGnew = atomKK->uCGnew;
  d_uCGnew = atomKK->k_uCGnew.d_view;
  h_uCGnew = atomKK->k_uCGnew.h_view;
  duChem = atomKK->duChem;
  d_duChem = atomKK->k_duChem.d_view;
  h_duChem = atomKK->k_duChem.h_view;
}

/* ----------------------------------------------------------------------
   sort atom arrays on device
------------------------------------------------------------------------- */

void AtomVecDPDKokkos::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  atomKK->sync(Device, ALL_MASK & ~F_MASK & ~DPDRHO_MASK & ~DUCHEM_MASK & ~DVECTOR_MASK);

  Sorter.sort(LMPDeviceType(), d_tag);
  Sorter.sort(LMPDeviceType(), d_type);
  Sorter.sort(LMPDeviceType(), d_mask);
  Sorter.sort(LMPDeviceType(), d_image);
  Sorter.sort(LMPDeviceType(), d_x);
  Sorter.sort(LMPDeviceType(), d_v);
  Sorter.sort(LMPDeviceType(), d_dpdTheta);
  Sorter.sort(LMPDeviceType(), d_uCond);
  Sorter.sort(LMPDeviceType(), d_uMech);
  Sorter.sort(LMPDeviceType(), d_uChem);
  Sorter.sort(LMPDeviceType(), d_uCG);
  Sorter.sort(LMPDeviceType(), d_uCGnew);

  atomKK->modified(Device, ALL_MASK & ~F_MASK & ~DPDRHO_MASK & ~DUCHEM_MASK & ~DVECTOR_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG,int TRICLINIC>
struct AtomVecDPDKokkos_PackComm {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_x_array_randomread _x;
  typename ArrayTypes<DeviceType>::t_efloat_1d _dpdTheta,_uCond,_uMech,_uChem;
  typename ArrayTypes<DeviceType>::t_xfloat_2d_um _buf;
  typename ArrayTypes<DeviceType>::t_int_1d_const _list;
  X_FLOAT _xprd,_yprd,_zprd,_xy,_xz,_yz;
  X_FLOAT _pbc[6];

  AtomVecDPDKokkos_PackComm(
      const typename DAT::tdual_x_array &x,
      const typename DAT::tdual_efloat_1d &dpdTheta,
      const typename DAT::tdual_efloat_1d &uCond,
      const typename DAT::tdual_efloat_1d &uMech,
      const typename DAT::tdual_efloat_1d &uChem,
      const typename DAT::tdual_xfloat_2d &buf,
      const typename DAT::tdual_int_1d &list,
      const X_FLOAT &xprd, const X_FLOAT &yprd, const X_FLOAT &zprd,
      const X_FLOAT &xy, const X_FLOAT &xz, const X_FLOAT &yz, const int* const pbc):
      _x(x.view<DeviceType>()),
      _dpdTheta(dpdTheta.view<DeviceType>()),
      _uCond(uCond.view<DeviceType>()),
      _uMech(uMech.view<DeviceType>()),
      _uChem(uChem.view<DeviceType>()),
      _list(list.view<DeviceType>()),
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
      _buf(i,3) = _dpdTheta(j);
      _buf(i,4) = _uCond(j);
      _buf(i,5) = _uMech(j);
      _buf(i,6) = _uChem(j);
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecDPDKokkos::pack_comm_kokkos(const int &n,
                                          const DAT::tdual_int_1d &list,
                                          const DAT::tdual_xfloat_2d &buf,
                                          const int &pbc_flag,
                                          const int* const pbc)
{
  // Check whether to always run forward communication on the host
  // Choose correct forward PackComm kernel

  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->sync(Host,X_MASK|DPDTHETA_MASK|UCOND_MASK|UMECH_MASK|UCHEM_MASK);
    if (pbc_flag) {
      if (domain->triclinic) {
        struct AtomVecDPDKokkos_PackComm<LMPHostType,1,1> f(atomKK->k_x,
          atomKK->k_dpdTheta,atomKK->k_uCond,atomKK->k_uMech,atomKK->k_uChem,
          buf,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecDPDKokkos_PackComm<LMPHostType,1,0> f(atomKK->k_x,
          atomKK->k_dpdTheta,atomKK->k_uCond,atomKK->k_uMech,atomKK->k_uChem,
          buf,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      }
    } else {
      if (domain->triclinic) {
        struct AtomVecDPDKokkos_PackComm<LMPHostType,0,1> f(atomKK->k_x,
          atomKK->k_dpdTheta,atomKK->k_uCond,atomKK->k_uMech,atomKK->k_uChem,
          buf,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecDPDKokkos_PackComm<LMPHostType,0,0> f(atomKK->k_x,
          atomKK->k_dpdTheta,atomKK->k_uCond,atomKK->k_uMech,atomKK->k_uChem,
          buf,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      }
    }
  } else {
    atomKK->sync(Device,X_MASK|DPDTHETA_MASK|UCOND_MASK|UMECH_MASK|UCHEM_MASK);
    if (pbc_flag) {
      if (domain->triclinic) {
        struct AtomVecDPDKokkos_PackComm<LMPDeviceType,1,1> f(atomKK->k_x,
          atomKK->k_dpdTheta,atomKK->k_uCond,atomKK->k_uMech,atomKK->k_uChem,
          buf,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecDPDKokkos_PackComm<LMPDeviceType,1,0> f(atomKK->k_x,
          atomKK->k_dpdTheta,atomKK->k_uCond,atomKK->k_uMech,atomKK->k_uChem,
          buf,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      }
    } else {
      if (domain->triclinic) {
        struct AtomVecDPDKokkos_PackComm<LMPDeviceType,0,1> f(atomKK->k_x,
          atomKK->k_dpdTheta,atomKK->k_uCond,atomKK->k_uMech,atomKK->k_uChem,
          buf,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecDPDKokkos_PackComm<LMPDeviceType,0,0> f(atomKK->k_x,
          atomKK->k_dpdTheta,atomKK->k_uCond,atomKK->k_uMech,atomKK->k_uChem,
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

template<class DeviceType,int PBC_FLAG,int TRICLINIC>
struct AtomVecDPDKokkos_PackCommSelf {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_x_array_randomread _x;
  typename ArrayTypes<DeviceType>::t_x_array _xw;
  typename ArrayTypes<DeviceType>::t_efloat_1d _dpdTheta,_uCond,_uMech,_uChem;
  int _nfirst;
  typename ArrayTypes<DeviceType>::t_int_1d_const _list;
  X_FLOAT _xprd,_yprd,_zprd,_xy,_xz,_yz;
  X_FLOAT _pbc[6];

  AtomVecDPDKokkos_PackCommSelf(
      const typename DAT::tdual_x_array &x,
      const typename DAT::tdual_efloat_1d &dpdTheta,
      const typename DAT::tdual_efloat_1d &uCond,
      const typename DAT::tdual_efloat_1d &uMech,
      const typename DAT::tdual_efloat_1d &uChem,
      const int &nfirst,
      const typename DAT::tdual_int_1d &list,
      const X_FLOAT &xprd, const X_FLOAT &yprd, const X_FLOAT &zprd,
      const X_FLOAT &xy, const X_FLOAT &xz, const X_FLOAT &yz, const int* const pbc):
      _x(x.view<DeviceType>()),_xw(x.view<DeviceType>()),
      _dpdTheta(dpdTheta.view<DeviceType>()),
      _uCond(uCond.view<DeviceType>()),
      _uMech(uMech.view<DeviceType>()),
      _uChem(uChem.view<DeviceType>()),
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
      _dpdTheta(i+_nfirst) = _dpdTheta(j);
      _uCond(i+_nfirst) = _uCond(j);
      _uMech(i+_nfirst) = _uMech(j);
      _uChem(i+_nfirst) = _uChem(j);
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecDPDKokkos::pack_comm_self(const int &n, const DAT::tdual_int_1d &list,
                                                                                const int nfirst, const int &pbc_flag, const int* const pbc) {
  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->sync(Host,X_MASK|DPDTHETA_MASK|UCOND_MASK|UMECH_MASK|UCHEM_MASK);
    atomKK->modified(Host,X_MASK|DPDTHETA_MASK|UCOND_MASK|UMECH_MASK|UCHEM_MASK);
    if (pbc_flag) {
      if (domain->triclinic) {
      struct AtomVecDPDKokkos_PackCommSelf<LMPHostType,1,1> f(atomKK->k_x,
          atomKK->k_dpdTheta,atomKK->k_uCond,atomKK->k_uMech,atomKK->k_uChem,
          nfirst,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      } else {
      struct AtomVecDPDKokkos_PackCommSelf<LMPHostType,1,0> f(atomKK->k_x,
          atomKK->k_dpdTheta,atomKK->k_uCond,atomKK->k_uMech,atomKK->k_uChem,
          nfirst,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      }
    } else {
      if (domain->triclinic) {
      struct AtomVecDPDKokkos_PackCommSelf<LMPHostType,0,1> f(atomKK->k_x,
          atomKK->k_dpdTheta,atomKK->k_uCond,atomKK->k_uMech,atomKK->k_uChem,
          nfirst,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      } else {
      struct AtomVecDPDKokkos_PackCommSelf<LMPHostType,0,0> f(atomKK->k_x,
          atomKK->k_dpdTheta,atomKK->k_uCond,atomKK->k_uMech,atomKK->k_uChem,
          nfirst,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      }
    }
  } else {
    atomKK->sync(Device,X_MASK|DPDTHETA_MASK|UCOND_MASK|UMECH_MASK|UCHEM_MASK);
    atomKK->modified(Device,X_MASK|DPDTHETA_MASK|UCOND_MASK|UMECH_MASK|UCHEM_MASK);
    if (pbc_flag) {
      if (domain->triclinic) {
      struct AtomVecDPDKokkos_PackCommSelf<LMPDeviceType,1,1> f(atomKK->k_x,
          atomKK->k_dpdTheta,atomKK->k_uCond,atomKK->k_uMech,atomKK->k_uChem,
          nfirst,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      } else {
      struct AtomVecDPDKokkos_PackCommSelf<LMPDeviceType,1,0> f(atomKK->k_x,
          atomKK->k_dpdTheta,atomKK->k_uCond,atomKK->k_uMech,atomKK->k_uChem,
          nfirst,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      }
    } else {
      if (domain->triclinic) {
      struct AtomVecDPDKokkos_PackCommSelf<LMPDeviceType,0,1> f(atomKK->k_x,
          atomKK->k_dpdTheta,atomKK->k_uCond,atomKK->k_uMech,atomKK->k_uChem,
          nfirst,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      } else {
      struct AtomVecDPDKokkos_PackCommSelf<LMPDeviceType,0,0> f(atomKK->k_x,
          atomKK->k_dpdTheta,atomKK->k_uCond,atomKK->k_uMech,atomKK->k_uChem,
          nfirst,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      }
    }
  }
        return n*3;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecDPDKokkos_UnpackComm {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_x_array _x;
  typename ArrayTypes<DeviceType>::t_efloat_1d _dpdTheta,_uCond,_uMech,_uChem;
  typename ArrayTypes<DeviceType>::t_xfloat_2d_const _buf;
  int _first;

  AtomVecDPDKokkos_UnpackComm(
      const typename DAT::tdual_x_array &x,
      const typename DAT::tdual_efloat_1d &dpdTheta,
      const typename DAT::tdual_efloat_1d &uCond,
      const typename DAT::tdual_efloat_1d &uMech,
      const typename DAT::tdual_efloat_1d &uChem,
      const typename DAT::tdual_xfloat_2d &buf,
      const int& first):_x(x.view<DeviceType>()),
                        _dpdTheta(dpdTheta.view<DeviceType>()),
                        _uCond(uCond.view<DeviceType>()),
                        _uMech(uMech.view<DeviceType>()),
                        _uChem(uChem.view<DeviceType>()),
                        _buf(buf.view<DeviceType>()),
                        _first(first) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
      _x(i+_first,0) = _buf(i,0);
      _x(i+_first,1) = _buf(i,1);
      _x(i+_first,2) = _buf(i,2);
      _dpdTheta(i+_first) = _buf(i,3);
      _uCond(i+_first) = _buf(i,4);
      _uMech(i+_first) = _buf(i,5);
      _uChem(i+_first) = _buf(i,6);
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecDPDKokkos::unpack_comm_kokkos(const int &n, const int &first,
    const DAT::tdual_xfloat_2d &buf) {
  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->sync(Host,X_MASK|DPDTHETA_MASK|UCOND_MASK|UMECH_MASK|UCHEM_MASK);
    atomKK->modified(Host,X_MASK|DPDTHETA_MASK|UCOND_MASK|UMECH_MASK|UCHEM_MASK);
    struct AtomVecDPDKokkos_UnpackComm<LMPHostType> f(atomKK->k_x,
    atomKK->k_dpdTheta,atomKK->k_uCond,atomKK->k_uMech,atomKK->k_uChem,
    buf,first);
    Kokkos::parallel_for(n,f);
  } else {
    atomKK->sync(Device,X_MASK|DPDTHETA_MASK|UCOND_MASK|UMECH_MASK|UCHEM_MASK);
    atomKK->modified(Device,X_MASK|DPDTHETA_MASK|UCOND_MASK|UMECH_MASK|UCHEM_MASK);
    struct AtomVecDPDKokkos_UnpackComm<LMPDeviceType> f(atomKK->k_x,
    atomKK->k_dpdTheta,atomKK->k_uCond,atomKK->k_uMech,atomKK->k_uChem,
    buf,first);
    Kokkos::parallel_for(n,f);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG>
struct AtomVecDPDKokkos_PackBorder {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_xfloat_2d _buf;
  const typename ArrayTypes<DeviceType>::t_int_1d_const _list;
  const typename ArrayTypes<DeviceType>::t_x_array_randomread _x;
  const typename ArrayTypes<DeviceType>::t_tagint_1d _tag;
  const typename ArrayTypes<DeviceType>::t_int_1d _type;
  const typename ArrayTypes<DeviceType>::t_int_1d _mask;
  typename ArrayTypes<DeviceType>::t_efloat_1d _dpdTheta,_uCond,_uMech,_uChem,_uCG,_uCGnew;
  X_FLOAT _dx,_dy,_dz;

  AtomVecDPDKokkos_PackBorder(
      const typename ArrayTypes<DeviceType>::t_xfloat_2d &buf,
      const typename ArrayTypes<DeviceType>::t_int_1d_const &list,
      const typename ArrayTypes<DeviceType>::t_x_array &x,
      const typename ArrayTypes<DeviceType>::t_tagint_1d &tag,
      const typename ArrayTypes<DeviceType>::t_int_1d &type,
      const typename ArrayTypes<DeviceType>::t_int_1d &mask,
      const typename ArrayTypes<DeviceType>::t_efloat_1d &dpdTheta,
      const typename ArrayTypes<DeviceType>::t_efloat_1d &uCond,
      const typename ArrayTypes<DeviceType>::t_efloat_1d &uMech,
      const typename ArrayTypes<DeviceType>::t_efloat_1d &uChem,
      const typename ArrayTypes<DeviceType>::t_efloat_1d &uCG,
      const typename ArrayTypes<DeviceType>::t_efloat_1d &uCGnew,
      const X_FLOAT &dx, const X_FLOAT &dy, const X_FLOAT &dz):
      _buf(buf),_list(list),
      _x(x),_tag(tag),_type(type),_mask(mask),
      _dpdTheta(dpdTheta),
      _uCond(uCond),
      _uMech(uMech),
      _uChem(uChem),
      _uCG(uCG),
      _uCGnew(uCGnew),
      _dx(dx),_dy(dy),_dz(dz) {}

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
      _buf(i,3) = _tag(j);
      _buf(i,4) = _type(j);
      _buf(i,5) = _mask(j);
      _buf(i,6) = _dpdTheta(j);
      _buf(i,7) = _uCond(j);
      _buf(i,8) = _uMech(j);
      _buf(i,9) = _uChem(j);
      _buf(i,10) = _uCG(j);
      _buf(i,11) = _uCGnew(j);
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecDPDKokkos::pack_border_kokkos(int n, DAT::tdual_int_1d k_sendlist, DAT::tdual_xfloat_2d buf,
                               int pbc_flag, int *pbc, ExecutionSpace space)
{
  X_FLOAT dx,dy,dz;

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
      AtomVecDPDKokkos_PackBorder<LMPHostType,1> f(
        buf.view<LMPHostType>(), k_sendlist.view<LMPHostType>(),
        h_x,h_tag,h_type,h_mask,
        h_dpdTheta,h_uCond,h_uMech,h_uChem,h_uCG,h_uCGnew,
        dx,dy,dz);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecDPDKokkos_PackBorder<LMPDeviceType,1> f(
        buf.view<LMPDeviceType>(), k_sendlist.view<LMPDeviceType>(),
        d_x,d_tag,d_type,d_mask,
        d_dpdTheta,d_uCond,d_uMech,d_uChem,d_uCG,d_uCGnew,
        dx,dy,dz);
      Kokkos::parallel_for(n,f);
    }

  } else {
    dx = dy = dz = 0;
    if (space==Host) {
      AtomVecDPDKokkos_PackBorder<LMPHostType,0> f(
        buf.view<LMPHostType>(), k_sendlist.view<LMPHostType>(),
        h_x,h_tag,h_type,h_mask,
        h_dpdTheta,h_uCond,h_uMech,h_uChem,h_uCG,h_uCGnew,
        dx,dy,dz);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecDPDKokkos_PackBorder<LMPDeviceType,0> f(
        buf.view<LMPDeviceType>(), k_sendlist.view<LMPDeviceType>(),
        d_x,d_tag,d_type,d_mask,
        d_dpdTheta,d_uCond,d_uMech,d_uChem,d_uCG,d_uCGnew,
        dx,dy,dz);
      Kokkos::parallel_for(n,f);
    }
  }
  return n*6;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecDPDKokkos_UnpackBorder {
  typedef DeviceType device_type;

  const typename ArrayTypes<DeviceType>::t_xfloat_2d_const _buf;
  typename ArrayTypes<DeviceType>::t_x_array _x;
  typename ArrayTypes<DeviceType>::t_tagint_1d _tag;
  typename ArrayTypes<DeviceType>::t_int_1d _type;
  typename ArrayTypes<DeviceType>::t_int_1d _mask;
  typename ArrayTypes<DeviceType>::t_efloat_1d _dpdTheta,_uCond,_uMech,_uChem,_uCG,_uCGnew;
  int _first;


  AtomVecDPDKokkos_UnpackBorder(
      const typename ArrayTypes<DeviceType>::t_xfloat_2d_const &buf,
      typename ArrayTypes<DeviceType>::t_x_array &x,
      typename ArrayTypes<DeviceType>::t_tagint_1d &tag,
      typename ArrayTypes<DeviceType>::t_int_1d &type,
      typename ArrayTypes<DeviceType>::t_int_1d &mask,
      const typename ArrayTypes<DeviceType>::t_efloat_1d &dpdTheta,
      const typename ArrayTypes<DeviceType>::t_efloat_1d &uCond,
      const typename ArrayTypes<DeviceType>::t_efloat_1d &uMech,
      const typename ArrayTypes<DeviceType>::t_efloat_1d &uChem,
      const typename ArrayTypes<DeviceType>::t_efloat_1d &uCG,
      const typename ArrayTypes<DeviceType>::t_efloat_1d &uCGnew,
      const int& first):
      _buf(buf),_x(x),_tag(tag),_type(type),_mask(mask),
      _dpdTheta(dpdTheta),
      _uCond(uCond),
      _uMech(uMech),
      _uChem(uChem),
      _uCG(uCG),
      _uCGnew(uCGnew),
      _first(first) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
      _x(i+_first,0) = _buf(i,0);
      _x(i+_first,1) = _buf(i,1);
      _x(i+_first,2) = _buf(i,2);
      _tag(i+_first) = static_cast<int> (_buf(i,3));
      _type(i+_first) = static_cast<int>  (_buf(i,4));
      _mask(i+_first) = static_cast<int>  (_buf(i,5));
      _dpdTheta(i+_first) = _buf(i,6);
      _uCond(i+_first) = _buf(i,7);
      _uMech(i+_first) = _buf(i,8);
      _uChem(i+_first) = _buf(i,9);
      _uCG(i+_first) = _buf(i,10);
      _uCGnew(i+_first) = _buf(i,11);
//      printf("%i %i %lf %lf %lf %i BORDER\n",_tag(i+_first),i+_first,_x(i+_first,0),_x(i+_first,1),_x(i+_first,2),_type(i+_first));
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecDPDKokkos::unpack_border_kokkos(const int &n, const int &first,
                     const DAT::tdual_xfloat_2d &buf,ExecutionSpace space) {
  atomKK->modified(space,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|
                 DPDTHETA_MASK|UCOND_MASK|UMECH_MASK|UCHEM_MASK|
                 UCG_MASK|UCGNEW_MASK);
  while (first+n >= nmax) grow(0);
  atomKK->modified(space,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|
                 DPDTHETA_MASK|UCOND_MASK|UMECH_MASK|UCHEM_MASK|
                 UCG_MASK|UCGNEW_MASK|DVECTOR_MASK);
  if (space==Host) {
    struct AtomVecDPDKokkos_UnpackBorder<LMPHostType> f(buf.view<LMPHostType>(),
      h_x,h_tag,h_type,h_mask,
      h_dpdTheta,h_uCond,h_uMech,h_uChem,h_uCG,h_uCGnew,
      first);
    Kokkos::parallel_for(n,f);
  } else {
    struct AtomVecDPDKokkos_UnpackBorder<LMPDeviceType> f(buf.view<LMPDeviceType>(),
      d_x,d_tag,d_type,d_mask,
      d_dpdTheta,d_uCond,d_uMech,d_uChem,d_uCG,d_uCGnew,
      first);
    Kokkos::parallel_for(n,f);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecDPDKokkos_PackExchangeFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_x_array_randomread _x;
  typename AT::t_v_array_randomread _v;
  typename AT::t_tagint_1d_randomread _tag;
  typename AT::t_int_1d_randomread _type;
  typename AT::t_int_1d_randomread _mask;
  typename AT::t_imageint_1d_randomread _image;
  typename AT::t_efloat_1d_randomread _dpdTheta,_uCond,_uMech,_uChem,_uCG,_uCGnew;
  typename AT::t_x_array _xw;
  typename AT::t_v_array _vw;
  typename AT::t_tagint_1d _tagw;
  typename AT::t_int_1d _typew;
  typename AT::t_int_1d _maskw;
  typename AT::t_imageint_1d _imagew;
  typename AT::t_efloat_1d _dpdThetaw,_uCondw,_uMechw,_uChemw,_uCGw,_uCGneww;

  typename AT::t_xfloat_2d_um _buf;
  typename AT::t_int_1d_const _sendlist;
  typename AT::t_int_1d_const _copylist;
  int _size_exchange;

  AtomVecDPDKokkos_PackExchangeFunctor(
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
                _dpdTheta(atom->k_dpdTheta.view<DeviceType>()),
                _uCond(atom->k_uCond.view<DeviceType>()),
                _uMech(atom->k_uMech.view<DeviceType>()),
                _uChem(atom->k_uChem.view<DeviceType>()),
                _uCG(atom->k_uCG.view<DeviceType>()),
                _uCGnew(atom->k_uCGnew.view<DeviceType>()),
                _xw(atom->k_x.view<DeviceType>()),
                _vw(atom->k_v.view<DeviceType>()),
                _tagw(atom->k_tag.view<DeviceType>()),
                _typew(atom->k_type.view<DeviceType>()),
                _maskw(atom->k_mask.view<DeviceType>()),
                _imagew(atom->k_image.view<DeviceType>()),
                _dpdThetaw(atom->k_dpdTheta.view<DeviceType>()),
                _uCondw(atom->k_uCond.view<DeviceType>()),
                _uMechw(atom->k_uMech.view<DeviceType>()),
                _uChemw(atom->k_uChem.view<DeviceType>()),
                _uCGw(atom->k_uCG.view<DeviceType>()),
                _uCGneww(atom->k_uCGnew.view<DeviceType>()),
                _sendlist(sendlist.template view<DeviceType>()),
                _copylist(copylist.template view<DeviceType>()),
                _size_exchange(atom->avecKK->size_exchange) {
    const int maxsendlist = (buf.template view<DeviceType>().extent(0)*buf.template view<DeviceType>().extent(1))/_size_exchange;

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
    _buf(mysend,7) = _tag[i];
    _buf(mysend,8) = _type[i];
    _buf(mysend,9) = _mask[i];
    _buf(mysend,10) = _image[i];
    _buf(mysend,11) = _dpdTheta[i];
    _buf(mysend,12) = _uCond[i];
    _buf(mysend,13) = _uMech[i];
    _buf(mysend,14) = _uChem[i];
    _buf(mysend,15) = _uCG[i];
    _buf(mysend,16) = _uCGnew[i];
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
    _dpdThetaw[i] = _dpdTheta(j);
    _uCondw[i] = _uCond(j);
    _uMechw[i] = _uMech(j);
    _uChemw[i] = _uChem(j);
    _uCGw[i] = _uCG(j);
    _uCGneww[i] = _uCGnew(j);
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecDPDKokkos::pack_exchange_kokkos(const int &nsend,DAT::tdual_xfloat_2d &k_buf, DAT::tdual_int_1d k_sendlist,DAT::tdual_int_1d k_copylist,ExecutionSpace space)
{
  size_exchange = 17;

  if (nsend > (int) (k_buf.view<LMPHostType>().extent(0)*k_buf.view<LMPHostType>().extent(1))/size_exchange) {
    int newsize = nsend*size_exchange/k_buf.view<LMPHostType>().extent(1)+1;
    k_buf.resize(newsize,k_buf.view<LMPHostType>().extent(1));
  }
  atomKK->sync(space,X_MASK | V_MASK | TAG_MASK | TYPE_MASK |
             MASK_MASK | IMAGE_MASK| DPDTHETA_MASK | UCOND_MASK |
             UMECH_MASK | UCHEM_MASK | UCG_MASK | UCGNEW_MASK |
             DVECTOR_MASK);
  if (space == Host) {
    AtomVecDPDKokkos_PackExchangeFunctor<LMPHostType> f(atomKK,k_buf,k_sendlist,k_copylist);
    Kokkos::parallel_for(nsend,f);
  } else {
    AtomVecDPDKokkos_PackExchangeFunctor<LMPDeviceType> f(atomKK,k_buf,k_sendlist,k_copylist);
    Kokkos::parallel_for(nsend,f);
  }
  return nsend*size_exchange;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecDPDKokkos_UnpackExchangeFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_x_array _x;
  typename AT::t_v_array _v;
  typename AT::t_tagint_1d _tag;
  typename AT::t_int_1d _type;
  typename AT::t_int_1d _mask;
  typename AT::t_imageint_1d _image;
  typename AT::t_efloat_1d _dpdTheta;
  typename AT::t_efloat_1d _uCond;
  typename AT::t_efloat_1d _uMech;
  typename AT::t_efloat_1d _uChem;
  typename AT::t_efloat_1d _uCG;
  typename AT::t_efloat_1d _uCGnew;

  typename AT::t_xfloat_2d_um _buf;
  typename AT::t_int_1d _nlocal;
  int _dim;
  X_FLOAT _lo,_hi;
  int _size_exchange;

  AtomVecDPDKokkos_UnpackExchangeFunctor(
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
                _nlocal(nlocal.template view<DeviceType>()),
                _dim(dim),_lo(lo),_hi(hi),_size_exchange(atom->avecKK->size_exchange) {
    const int maxsendlist = (buf.template view<DeviceType>().extent(0)*buf.template view<DeviceType>().extent(1))/_size_exchange;

    buffer_view<DeviceType>(_buf,buf,maxsendlist,_size_exchange);
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
      _tag[i] = _buf(myrecv,7);
      _type[i] = _buf(myrecv,8);
      _mask[i] = _buf(myrecv,9);
      _image[i] = _buf(myrecv,10);
      _dpdTheta[i] = _buf(myrecv,11);
      _uCond[i] = _buf(myrecv,12);
      _uMech[i] = _buf(myrecv,13);
      _uChem[i] = _buf(myrecv,14);
      _uCG[i] = _buf(myrecv,15);
      _uCGnew[i] = _buf(myrecv,16);
    }
  }
};

/* ---------------------------------------------------------------------- */
int AtomVecDPDKokkos::unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf, int nrecv, int nlocal,
                                             int dim, X_FLOAT lo, X_FLOAT hi, ExecutionSpace space,
                                             DAT::tdual_int_1d &/*k_indices*/)
{
  while (nlocal + nrecv/size_exchange >= nmax) grow(0);

  if (space == Host) {
    k_count.h_view(0) = nlocal;
    AtomVecDPDKokkos_UnpackExchangeFunctor<LMPHostType> f(atomKK,k_buf,k_count,dim,lo,hi);
    Kokkos::parallel_for(nrecv/size_exchange,f);
  } else {
    k_count.h_view(0) = nlocal;
    k_count.modify<LMPHostType>();
    k_count.sync<LMPDeviceType>();
    AtomVecDPDKokkos_UnpackExchangeFunctor<LMPDeviceType> f(atomKK,k_buf,k_count,dim,lo,hi);
    Kokkos::parallel_for(nrecv/size_exchange,f);
    k_count.modify<LMPDeviceType>();
    k_count.sync<LMPHostType>();
  }

  atomKK->modified(space,X_MASK | V_MASK | TAG_MASK | TYPE_MASK |
                 MASK_MASK | IMAGE_MASK| DPDTHETA_MASK | UCOND_MASK |
                 UMECH_MASK | UCHEM_MASK | UCG_MASK | UCGNEW_MASK |
                 DVECTOR_MASK);

  return k_count.h_view(0);
}

/* ---------------------------------------------------------------------- */

void AtomVecDPDKokkos::sync(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.sync<LMPDeviceType>();
    if (mask & V_MASK) atomKK->k_v.sync<LMPDeviceType>();
    if (mask & F_MASK) atomKK->k_f.sync<LMPDeviceType>();
    if (mask & TAG_MASK) atomKK->k_tag.sync<LMPDeviceType>();
    if (mask & TYPE_MASK) atomKK->k_type.sync<LMPDeviceType>();
    if (mask & MASK_MASK) atomKK->k_mask.sync<LMPDeviceType>();
    if (mask & IMAGE_MASK) atomKK->k_image.sync<LMPDeviceType>();
    if (mask & DPDRHO_MASK) atomKK->k_rho.sync<LMPDeviceType>();
    if (mask & DPDTHETA_MASK) atomKK->k_dpdTheta.sync<LMPDeviceType>();
    if (mask & UCOND_MASK) atomKK->k_uCond.sync<LMPDeviceType>();
    if (mask & UMECH_MASK) atomKK->k_uMech.sync<LMPDeviceType>();
    if (mask & UCHEM_MASK) atomKK->k_uChem.sync<LMPDeviceType>();
    if (mask & UCG_MASK) atomKK->k_uCG.sync<LMPDeviceType>();
    if (mask & UCGNEW_MASK) atomKK->k_uCGnew.sync<LMPDeviceType>();
    if (mask & DUCHEM_MASK) atomKK->k_duChem.sync<LMPDeviceType>();
  } else {
    if (mask & X_MASK) atomKK->k_x.sync<LMPHostType>();
    if (mask & V_MASK) atomKK->k_v.sync<LMPHostType>();
    if (mask & F_MASK) atomKK->k_f.sync<LMPHostType>();
    if (mask & TAG_MASK) atomKK->k_tag.sync<LMPHostType>();
    if (mask & TYPE_MASK) atomKK->k_type.sync<LMPHostType>();
    if (mask & MASK_MASK) atomKK->k_mask.sync<LMPHostType>();
    if (mask & IMAGE_MASK) atomKK->k_image.sync<LMPHostType>();
    if (mask & DPDRHO_MASK) atomKK->k_rho.sync<LMPHostType>();
    if (mask & DPDTHETA_MASK) atomKK->k_dpdTheta.sync<LMPHostType>();
    if (mask & UCOND_MASK) atomKK->k_uCond.sync<LMPHostType>();
    if (mask & UMECH_MASK) atomKK->k_uMech.sync<LMPHostType>();
    if (mask & UCHEM_MASK) atomKK->k_uChem.sync<LMPHostType>();
    if (mask & UCG_MASK) atomKK->k_uCG.sync<LMPHostType>();
    if (mask & UCGNEW_MASK) atomKK->k_uCGnew.sync<LMPHostType>();
    if (mask & DUCHEM_MASK) atomKK->k_duChem.sync<LMPHostType>();
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecDPDKokkos::sync_overlapping_device(ExecutionSpace space, unsigned int mask)
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
    if ((mask & DPDRHO_MASK) && atomKK->k_rho.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_efloat_1d>(atomKK->k_rho,space);
    if ((mask & DPDTHETA_MASK) && atomKK->k_dpdTheta.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_efloat_1d>(atomKK->k_dpdTheta,space);
    if ((mask & UCOND_MASK) && atomKK->k_uCond.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_efloat_1d>(atomKK->k_uCond,space);
    if ((mask & UMECH_MASK) && atomKK->k_uMech.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_efloat_1d>(atomKK->k_uMech,space);
    if ((mask & UCHEM_MASK) && atomKK->k_uChem.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_efloat_1d>(atomKK->k_uChem,space);
    if ((mask & UCG_MASK) && atomKK->k_uCG.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_efloat_1d>(atomKK->k_uCG,space);
    if ((mask & UCGNEW_MASK) && atomKK->k_uCGnew.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_efloat_1d>(atomKK->k_uCGnew,space);
    if ((mask & DUCHEM_MASK) && atomKK->k_duChem.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_efloat_1d>(atomKK->k_duChem,space);
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
    if ((mask & DPDRHO_MASK) && atomKK->k_rho.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_efloat_1d>(atomKK->k_rho,space);
    if ((mask & DPDTHETA_MASK) && atomKK->k_dpdTheta.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_efloat_1d>(atomKK->k_dpdTheta,space);
    if ((mask & UCOND_MASK) && atomKK->k_uCond.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_efloat_1d>(atomKK->k_uCond,space);
    if ((mask & UMECH_MASK) && atomKK->k_uMech.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_efloat_1d>(atomKK->k_uMech,space);
    if ((mask & UCHEM_MASK) && atomKK->k_uChem.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_efloat_1d>(atomKK->k_uChem,space);
    if ((mask & UCG_MASK) && atomKK->k_uCG.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_efloat_1d>(atomKK->k_uCG,space);
    if ((mask & UCGNEW_MASK) && atomKK->k_uCGnew.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_efloat_1d>(atomKK->k_uCGnew,space);
    if ((mask & DUCHEM_MASK) && atomKK->k_duChem.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_efloat_1d>(atomKK->k_duChem,space);
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecDPDKokkos::modified(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.modify<LMPDeviceType>();
    if (mask & V_MASK) atomKK->k_v.modify<LMPDeviceType>();
    if (mask & F_MASK) atomKK->k_f.modify<LMPDeviceType>();
    if (mask & TAG_MASK) atomKK->k_tag.modify<LMPDeviceType>();
    if (mask & TYPE_MASK) atomKK->k_type.modify<LMPDeviceType>();
    if (mask & MASK_MASK) atomKK->k_mask.modify<LMPDeviceType>();
    if (mask & IMAGE_MASK) atomKK->k_image.modify<LMPDeviceType>();
    if (mask & DPDRHO_MASK) atomKK->k_rho.modify<LMPDeviceType>();
    if (mask & DPDTHETA_MASK) atomKK->k_dpdTheta.modify<LMPDeviceType>();
    if (mask & UCOND_MASK) atomKK->k_uCond.modify<LMPDeviceType>();
    if (mask & UMECH_MASK) atomKK->k_uMech.modify<LMPDeviceType>();
    if (mask & UCHEM_MASK) atomKK->k_uChem.modify<LMPDeviceType>();
    if (mask & UCG_MASK) atomKK->k_uCG.modify<LMPDeviceType>();
    if (mask & UCGNEW_MASK) atomKK->k_uCGnew.modify<LMPDeviceType>();
    if (mask & DUCHEM_MASK) atomKK->k_duChem.modify<LMPDeviceType>();
  } else {
    if (mask & X_MASK) atomKK->k_x.modify<LMPHostType>();
    if (mask & V_MASK) atomKK->k_v.modify<LMPHostType>();
    if (mask & F_MASK) atomKK->k_f.modify<LMPHostType>();
    if (mask & TAG_MASK) atomKK->k_tag.modify<LMPHostType>();
    if (mask & TYPE_MASK) atomKK->k_type.modify<LMPHostType>();
    if (mask & MASK_MASK) atomKK->k_mask.modify<LMPHostType>();
    if (mask & IMAGE_MASK) atomKK->k_image.modify<LMPHostType>();
    if (mask & DPDRHO_MASK) atomKK->k_rho.modify<LMPHostType>();
    if (mask & DPDTHETA_MASK) atomKK->k_dpdTheta.modify<LMPHostType>();
    if (mask & UCOND_MASK) atomKK->k_uCond.modify<LMPHostType>();
    if (mask & UMECH_MASK) atomKK->k_uMech.modify<LMPHostType>();
    if (mask & UCHEM_MASK) atomKK->k_uChem.modify<LMPHostType>();
    if (mask & UCG_MASK) atomKK->k_uCG.modify<LMPHostType>();
    if (mask & UCGNEW_MASK) atomKK->k_uCGnew.modify<LMPHostType>();
    if (mask & DUCHEM_MASK) atomKK->k_duChem.modify<LMPHostType>();
  }
}
