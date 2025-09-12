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
  h_x = atomKK->k_x.h_viewkk;
  v = atomKK->v;
  d_v = atomKK->k_v.d_view;
  h_v = atomKK->k_v.h_viewkk;
  f = atomKK->f;
  d_f = atomKK->k_f.d_view;
  h_f = atomKK->k_f.h_viewkk;

  rho = atomKK->rho;
  d_rho = atomKK->k_rho.d_view;
  h_rho = atomKK->k_rho.h_viewkk;
  dpdTheta = atomKK->dpdTheta;
  d_dpdTheta = atomKK->k_dpdTheta.d_view;
  h_dpdTheta = atomKK->k_dpdTheta.h_viewkk;
  uCond = atomKK->uCond;
  d_uCond = atomKK->k_uCond.d_view;
  h_uCond = atomKK->k_uCond.h_viewkk;
  uMech = atomKK->uMech;
  d_uMech = atomKK->k_uMech.d_view;
  h_uMech = atomKK->k_uMech.h_viewkk;
  uChem = atomKK->uChem;
  d_uChem = atomKK->k_uChem.d_view;
  h_uChem = atomKK->k_uChem.h_viewkk;
  uCG = atomKK->uCG;
  d_uCG = atomKK->k_uCG.d_view;
  h_uCG = atomKK->k_uCG.h_viewkk;
  uCGnew = atomKK->uCGnew;
  d_uCGnew = atomKK->k_uCGnew.d_view;
  h_uCGnew = atomKK->k_uCGnew.h_viewkk;
  duChem = atomKK->duChem;
  d_duChem = atomKK->k_duChem.d_view;
  h_duChem = atomKK->k_duChem.h_viewkk;
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
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr_randomread _x;
  typename AT::t_kkfloat_1d _dpdTheta,_uCond,_uMech,_uChem;
  typename AT::t_double_2d_lr_um _buf;
  typename AT::t_int_1d_const _list;
  double _xprd,_yprd,_zprd,_xy,_xz,_yz;
  double _pbc[6];

  AtomVecDPDKokkos_PackComm(
      const typename DAT::ttransform_kkfloat_1d_3_lr &x,
      const typename DAT::ttransform_kkfloat_1d &dpdTheta,
      const typename DAT::ttransform_kkfloat_1d &uCond,
      const typename DAT::ttransform_kkfloat_1d &uMech,
      const typename DAT::ttransform_kkfloat_1d &uChem,
      const typename DAT::tdual_double_2d_lr &buf,
      const typename DAT::tdual_int_1d &list,
      const double &xprd, const double &yprd, const double &zprd,
      const double &xy, const double &xz, const double &yz, const int* const pbc):
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
                                          const DAT::tdual_double_2d_lr &buf,
                                          const int &pbc_flag,
                                          const int* const pbc)
{
  // Check whether to always run forward communication on the host
  // Choose correct forward PackComm kernel

  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->sync(HostKK,X_MASK|DPDTHETA_MASK|UCOND_MASK|UMECH_MASK|UCHEM_MASK);
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
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr_randomread _x;
  typename AT::t_kkfloat_1d_3_lr _xw;
  typename AT::t_kkfloat_1d _dpdTheta,_uCond,_uMech,_uChem;
  int _nfirst;
  typename AT::t_int_1d_const _list;
  double _xprd,_yprd,_zprd,_xy,_xz,_yz;
  double _pbc[6];

  AtomVecDPDKokkos_PackCommSelf(
      const typename DAT::ttransform_kkfloat_1d_3_lr &x,
      const typename DAT::ttransform_kkfloat_1d &dpdTheta,
      const typename DAT::ttransform_kkfloat_1d &uCond,
      const typename DAT::ttransform_kkfloat_1d &uMech,
      const typename DAT::ttransform_kkfloat_1d &uChem,
      const int &nfirst,
      const typename DAT::tdual_int_1d &list,
      const double &xprd, const double &yprd, const double &zprd,
      const double &xy, const double &xz, const double &yz, const int* const pbc):
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
    atomKK->sync(HostKK,X_MASK|DPDTHETA_MASK|UCOND_MASK|UMECH_MASK|UCHEM_MASK);
    atomKK->modified(HostKK,X_MASK|DPDTHETA_MASK|UCOND_MASK|UMECH_MASK|UCHEM_MASK);
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
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr _x;
  typename AT::t_kkfloat_1d _dpdTheta,_uCond,_uMech,_uChem;
  typename AT::t_double_2d_lr_const _buf;
  int _first;

  AtomVecDPDKokkos_UnpackComm(
      const typename DAT::ttransform_kkfloat_1d_3_lr &x,
      const typename DAT::ttransform_kkfloat_1d &dpdTheta,
      const typename DAT::ttransform_kkfloat_1d &uCond,
      const typename DAT::ttransform_kkfloat_1d &uMech,
      const typename DAT::ttransform_kkfloat_1d &uChem,
      const typename DAT::tdual_double_2d_lr &buf,
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
    const DAT::tdual_double_2d_lr &buf) {
  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->sync(HostKK,X_MASK|DPDTHETA_MASK|UCOND_MASK|UMECH_MASK|UCHEM_MASK);
    atomKK->modified(HostKK,X_MASK|DPDTHETA_MASK|UCOND_MASK|UMECH_MASK|UCHEM_MASK);
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
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_double_2d_lr _buf;
  const typename AT::t_int_1d_const _list;
  const typename AT::t_kkfloat_1d_3_lr_randomread _x;
  const typename AT::t_tagint_1d _tag;
  const typename AT::t_int_1d _type;
  const typename AT::t_int_1d _mask;
  typename AT::t_kkfloat_1d _dpdTheta,_uCond,_uMech,_uChem,_uCG,_uCGnew;
  double _dx,_dy,_dz;

  AtomVecDPDKokkos_PackBorder(
      const typename AT::t_double_2d_lr &buf,
      const typename AT::t_int_1d_const &list,
      const typename AT::t_kkfloat_1d_3_lr &x,
      const typename AT::t_tagint_1d &tag,
      const typename AT::t_int_1d &type,
      const typename AT::t_int_1d &mask,
      const typename AT::t_kkfloat_1d &dpdTheta,
      const typename AT::t_kkfloat_1d &uCond,
      const typename AT::t_kkfloat_1d &uMech,
      const typename AT::t_kkfloat_1d &uChem,
      const typename AT::t_kkfloat_1d &uCG,
      const typename AT::t_kkfloat_1d &uCGnew,
      const double &dx, const double &dy, const double &dz):
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

int AtomVecDPDKokkos::pack_border_kokkos(int n, DAT::tdual_int_1d k_sendlist, DAT::tdual_double_2d_lr buf,
                               int pbc_flag, int *pbc, ExecutionSpace space)
{
  double dx,dy,dz;

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
        buf.h_view, k_sendlist.h_view,
        h_x,h_tag,h_type,h_mask,
        h_dpdTheta,h_uCond,h_uMech,h_uChem,h_uCG,h_uCGnew,
        dx,dy,dz);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecDPDKokkos_PackBorder<LMPDeviceType,1> f(
        buf.d_view, k_sendlist.d_view,
        d_x,d_tag,d_type,d_mask,
        d_dpdTheta,d_uCond,d_uMech,d_uChem,d_uCG,d_uCGnew,
        dx,dy,dz);
      Kokkos::parallel_for(n,f);
    }

  } else {
    dx = dy = dz = 0;
    if (space==Host) {
      AtomVecDPDKokkos_PackBorder<LMPHostType,0> f(
        buf.h_view, k_sendlist.h_view,
        h_x,h_tag,h_type,h_mask,
        h_dpdTheta,h_uCond,h_uMech,h_uChem,h_uCG,h_uCGnew,
        dx,dy,dz);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecDPDKokkos_PackBorder<LMPDeviceType,0> f(
        buf.d_view, k_sendlist.d_view,
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
  typedef ArrayTypes<DeviceType> AT;

  const typename AT::t_double_2d_lr_const _buf;
  typename AT::t_kkfloat_1d_3_lr _x;
  typename AT::t_tagint_1d _tag;
  typename AT::t_int_1d _type;
  typename AT::t_int_1d _mask;
  typename AT::t_kkfloat_1d _dpdTheta,_uCond,_uMech,_uChem,_uCG,_uCGnew;
  int _first;


  AtomVecDPDKokkos_UnpackBorder(
      const typename AT::t_double_2d_lr_const &buf,
      typename AT::t_kkfloat_1d_3_lr &x,
      typename AT::t_tagint_1d &tag,
      typename AT::t_int_1d &type,
      typename AT::t_int_1d &mask,
      const typename AT::t_kkfloat_1d &dpdTheta,
      const typename AT::t_kkfloat_1d &uCond,
      const typename AT::t_kkfloat_1d &uMech,
      const typename AT::t_kkfloat_1d &uChem,
      const typename AT::t_kkfloat_1d &uCG,
      const typename AT::t_kkfloat_1d &uCGnew,
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
                     const DAT::tdual_double_2d_lr &buf,ExecutionSpace space) {
  atomKK->modified(space,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|
                 DPDTHETA_MASK|UCOND_MASK|UMECH_MASK|UCHEM_MASK|
                 UCG_MASK|UCGNEW_MASK);
  while (first+n >= nmax) grow(0);
  atomKK->modified(space,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|
                 DPDTHETA_MASK|UCOND_MASK|UMECH_MASK|UCHEM_MASK|
                 UCG_MASK|UCGNEW_MASK|DVECTOR_MASK);
  if (space==Host) {
    struct AtomVecDPDKokkos_UnpackBorder<LMPHostType> f(buf.h_view,
      h_x,h_tag,h_type,h_mask,
      h_dpdTheta,h_uCond,h_uMech,h_uChem,h_uCG,h_uCGnew,
      first);
    Kokkos::parallel_for(n,f);
  } else {
    struct AtomVecDPDKokkos_UnpackBorder<LMPDeviceType> f(buf.d_view,
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
  typename AT::t_kkfloat_1d_3_lr_randomread _x;
  typename AT::t_kkfloat_1d_3_randomread _v;
  typename AT::t_tagint_1d_randomread _tag;
  typename AT::t_int_1d_randomread _type;
  typename AT::t_int_1d_randomread _mask;
  typename AT::t_imageint_1d_randomread _image;
  typename AT::t_kkfloat_1d_randomread _dpdTheta,_uCond,_uMech,_uChem,_uCG,_uCGnew;
  typename AT::t_kkfloat_1d_3_lr _xw;
  typename AT::t_kkfloat_1d_3 _vw;
  typename AT::t_tagint_1d _tagw;
  typename AT::t_int_1d _typew;
  typename AT::t_int_1d _maskw;
  typename AT::t_imageint_1d _imagew;
  typename AT::t_kkfloat_1d _dpdThetaw,_uCondw,_uMechw,_uChemw,_uCGw,_uCGneww;

  typename AT::t_double_2d_lr_um _buf;
  typename AT::t_int_1d_const _sendlist;
  typename AT::t_int_1d_const _copylist;
  int _size_exchange;

  AtomVecDPDKokkos_PackExchangeFunctor(
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

int AtomVecDPDKokkos::pack_exchange_kokkos(const int &nsend,DAT::tdual_double_2d_lr &k_buf, DAT::tdual_int_1d k_sendlist,DAT::tdual_int_1d k_copylist,ExecutionSpace space)
{
  size_exchange = 17;

  if (nsend > (int) (k_buf.h_view.extent(0)*k_buf.h_view.extent(1))/size_exchange) {
    int newsize = nsend*size_exchange/k_buf.h_view.extent(1)+1;
    k_buf.resize(newsize,k_buf.h_view.extent(1));
  }
  atomKK->sync(space,X_MASK | V_MASK | TAG_MASK | TYPE_MASK |
             MASK_MASK | IMAGE_MASK| DPDTHETA_MASK | UCOND_MASK |
             UMECH_MASK | UCHEM_MASK | UCG_MASK | UCGNEW_MASK |
             DVECTOR_MASK);
  if (space == HostKK) {
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
  typename AT::t_kkfloat_1d_3_lr _x;
  typename AT::t_kkfloat_1d_3 _v;
  typename AT::t_tagint_1d _tag;
  typename AT::t_int_1d _type;
  typename AT::t_int_1d _mask;
  typename AT::t_imageint_1d _image;
  typename AT::t_kkfloat_1d _dpdTheta;
  typename AT::t_kkfloat_1d _uCond;
  typename AT::t_kkfloat_1d _uMech;
  typename AT::t_kkfloat_1d _uChem;
  typename AT::t_kkfloat_1d _uCG;
  typename AT::t_kkfloat_1d _uCGnew;

  typename AT::t_double_2d_lr_um _buf;
  typename AT::t_int_1d _nlocal;
  int _dim;
  double _lo,_hi;
  int _size_exchange;

  AtomVecDPDKokkos_UnpackExchangeFunctor(
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
int AtomVecDPDKokkos::unpack_exchange_kokkos(DAT::tdual_double_2d_lr &k_buf, int nrecv, int nlocal,
                                             int dim, double lo, double hi, ExecutionSpace space,
                                             DAT::tdual_int_1d &/*k_indices*/)
{
  while (nlocal + nrecv/size_exchange >= nmax) grow(0);

  if (space == HostKK) {
    k_count.h_view(0) = nlocal;
    AtomVecDPDKokkos_UnpackExchangeFunctor<LMPHostType> f(atomKK,k_buf,k_count,dim,lo,hi);
    Kokkos::parallel_for(nrecv/size_exchange,f);
  } else {
    k_count.h_view(0) = nlocal;
    k_count.modify_host();
    k_count.sync_device();
    AtomVecDPDKokkos_UnpackExchangeFunctor<LMPDeviceType> f(atomKK,k_buf,k_count,dim,lo,hi);
    Kokkos::parallel_for(nrecv/size_exchange,f);
    k_count.modify_device();
    k_count.sync_host();
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
    if (mask & X_MASK) atomKK->k_x.sync_device();
    if (mask & V_MASK) atomKK->k_v.sync_device();
    if (mask & F_MASK) atomKK->k_f.sync_device();
    if (mask & TAG_MASK) atomKK->k_tag.sync_device();
    if (mask & TYPE_MASK) atomKK->k_type.sync_device();
    if (mask & MASK_MASK) atomKK->k_mask.sync_device();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_device();
    if (mask & DPDRHO_MASK) atomKK->k_rho.sync_device();
    if (mask & DPDTHETA_MASK) atomKK->k_dpdTheta.sync_device();
    if (mask & UCOND_MASK) atomKK->k_uCond.sync_device();
    if (mask & UMECH_MASK) atomKK->k_uMech.sync_device();
    if (mask & UCHEM_MASK) atomKK->k_uChem.sync_device();
    if (mask & UCG_MASK) atomKK->k_uCG.sync_device();
    if (mask & UCGNEW_MASK) atomKK->k_uCGnew.sync_device();
    if (mask & DUCHEM_MASK) atomKK->k_duChem.sync_device();
  } else if (space == Host) {
    if (mask & X_MASK) atomKK->k_x.sync_host();
    if (mask & V_MASK) atomKK->k_v.sync_host();
    if (mask & F_MASK) atomKK->k_f.sync_host();
    if (mask & TAG_MASK) atomKK->k_tag.sync_host();
    if (mask & TYPE_MASK) atomKK->k_type.sync_host();
    if (mask & MASK_MASK) atomKK->k_mask.sync_host();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_host();
    if (mask & DPDRHO_MASK) atomKK->k_rho.sync_host();
    if (mask & DPDTHETA_MASK) atomKK->k_dpdTheta.sync_host();
    if (mask & UCOND_MASK) atomKK->k_uCond.sync_host();
    if (mask & UMECH_MASK) atomKK->k_uMech.sync_host();
    if (mask & UCHEM_MASK) atomKK->k_uChem.sync_host();
    if (mask & UCG_MASK) atomKK->k_uCG.sync_host();
    if (mask & UCGNEW_MASK) atomKK->k_uCGnew.sync_host();
    if (mask & DUCHEM_MASK) atomKK->k_duChem.sync_host();
  } else if (space == HostKK) {
    if (mask & X_MASK) atomKK->k_x.sync_hostkk();
    if (mask & V_MASK) atomKK->k_v.sync_hostkk();
    if (mask & F_MASK) atomKK->k_f.sync_hostkk();
    if (mask & TAG_MASK) atomKK->k_tag.sync_host();
    if (mask & TYPE_MASK) atomKK->k_type.sync_host();
    if (mask & MASK_MASK) atomKK->k_mask.sync_host();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_host();
    if (mask & DPDRHO_MASK) atomKK->k_rho.sync_hostkk();
    if (mask & DPDTHETA_MASK) atomKK->k_dpdTheta.sync_hostkk();
    if (mask & UCOND_MASK) atomKK->k_uCond.sync_hostkk();
    if (mask & UMECH_MASK) atomKK->k_uMech.sync_hostkk();
    if (mask & UCHEM_MASK) atomKK->k_uChem.sync_hostkk();
    if (mask & UCG_MASK) atomKK->k_uCG.sync_hostkk();
    if (mask & UCGNEW_MASK) atomKK->k_uCGnew.sync_hostkk();
    if (mask & DUCHEM_MASK) atomKK->k_duChem.sync_hostkk();
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecDPDKokkos::sync_pinned(ExecutionSpace space, unsigned int mask, int async_flag)
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
    if ((mask & DPDRHO_MASK) && atomKK->k_rho.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_rho,space,async_flag);
    if ((mask & DPDTHETA_MASK) && atomKK->k_dpdTheta.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_dpdTheta,space,async_flag);
    if ((mask & UCOND_MASK) && atomKK->k_uCond.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_uCond,space,async_flag);
    if ((mask & UMECH_MASK) && atomKK->k_uMech.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_uMech,space,async_flag);
    if ((mask & UCHEM_MASK) && atomKK->k_uChem.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_uChem,space,async_flag);
    if ((mask & UCG_MASK) && atomKK->k_uCG.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_uCG,space,async_flag);
    if ((mask & UCGNEW_MASK) && atomKK->k_uCGnew.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_uCGnew,space,async_flag);
    if ((mask & DUCHEM_MASK) && atomKK->k_duChem.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_duChem,space,async_flag);
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
    if ((mask & DPDRHO_MASK) && atomKK->k_rho.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_rho,space,async_flag);
    if ((mask & DPDTHETA_MASK) && atomKK->k_dpdTheta.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_dpdTheta,space,async_flag);
    if ((mask & UCOND_MASK) && atomKK->k_uCond.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_uCond,space,async_flag);
    if ((mask & UMECH_MASK) && atomKK->k_uMech.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_uMech,space,async_flag);
    if ((mask & UCHEM_MASK) && atomKK->k_uChem.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_uChem,space,async_flag);
    if ((mask & UCG_MASK) && atomKK->k_uCG.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_uCG,space,async_flag);
    if ((mask & UCGNEW_MASK) && atomKK->k_uCGnew.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_uCGnew,space,async_flag);
    if ((mask & DUCHEM_MASK) && atomKK->k_duChem.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_duChem,space,async_flag);
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecDPDKokkos::modified(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.modify_device();
    if (mask & V_MASK) atomKK->k_v.modify_device();
    if (mask & F_MASK) atomKK->k_f.modify_device();
    if (mask & TAG_MASK) atomKK->k_tag.modify_device();
    if (mask & TYPE_MASK) atomKK->k_type.modify_device();
    if (mask & MASK_MASK) atomKK->k_mask.modify_device();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_device();
    if (mask & DPDRHO_MASK) atomKK->k_rho.modify_device();
    if (mask & DPDTHETA_MASK) atomKK->k_dpdTheta.modify_device();
    if (mask & UCOND_MASK) atomKK->k_uCond.modify_device();
    if (mask & UMECH_MASK) atomKK->k_uMech.modify_device();
    if (mask & UCHEM_MASK) atomKK->k_uChem.modify_device();
    if (mask & UCG_MASK) atomKK->k_uCG.modify_device();
    if (mask & UCGNEW_MASK) atomKK->k_uCGnew.modify_device();
    if (mask & DUCHEM_MASK) atomKK->k_duChem.modify_device();
  } else if (space == Host) {
    if (mask & X_MASK) atomKK->k_x.modify_host();
    if (mask & V_MASK) atomKK->k_v.modify_host();
    if (mask & F_MASK) atomKK->k_f.modify_host();
    if (mask & TAG_MASK) atomKK->k_tag.modify_host();
    if (mask & TYPE_MASK) atomKK->k_type.modify_host();
    if (mask & MASK_MASK) atomKK->k_mask.modify_host();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_host();
    if (mask & DPDRHO_MASK) atomKK->k_rho.modify_host();
    if (mask & DPDTHETA_MASK) atomKK->k_dpdTheta.modify_host();
    if (mask & UCOND_MASK) atomKK->k_uCond.modify_host();
    if (mask & UMECH_MASK) atomKK->k_uMech.modify_host();
    if (mask & UCHEM_MASK) atomKK->k_uChem.modify_host();
    if (mask & UCG_MASK) atomKK->k_uCG.modify_host();
    if (mask & UCGNEW_MASK) atomKK->k_uCGnew.modify_host();
    if (mask & DUCHEM_MASK) atomKK->k_duChem.modify_host();
  } else if (space == HostKK) {
    if (mask & X_MASK) atomKK->k_x.modify_hostkk();
    if (mask & V_MASK) atomKK->k_v.modify_hostkk();
    if (mask & F_MASK) atomKK->k_f.modify_hostkk();
    if (mask & TAG_MASK) atomKK->k_tag.modify_host();
    if (mask & TYPE_MASK) atomKK->k_type.modify_host();
    if (mask & MASK_MASK) atomKK->k_mask.modify_host();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_host();
    if (mask & DPDRHO_MASK) atomKK->k_rho.modify_hostkk();
    if (mask & DPDTHETA_MASK) atomKK->k_dpdTheta.modify_hostkk();
    if (mask & UCOND_MASK) atomKK->k_uCond.modify_hostkk();
    if (mask & UMECH_MASK) atomKK->k_uMech.modify_hostkk();
    if (mask & UCHEM_MASK) atomKK->k_uChem.modify_hostkk();
    if (mask & UCG_MASK) atomKK->k_uCG.modify_hostkk();
    if (mask & UCGNEW_MASK) atomKK->k_uCGnew.modify_hostkk();
    if (mask & DUCHEM_MASK) atomKK->k_duChem.modify_hostkk();
  }
}
