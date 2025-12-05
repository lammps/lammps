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

#include "atom_vec_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm_kokkos.h"
#include "domain.h"
#include "error.h"
#include "kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecKokkos::AtomVecKokkos(LAMMPS *lmp) : AtomVec(lmp)
{
  if (!lmp->kokkos || !lmp->kokkos->kokkos_exists)
    error->all(FLERR, Error::NOLASTLINE, "Cannot use KOKKOS styles without enabling KOKKOS");

  kokkosable = 1;
  buffer = nullptr;
  buffer_size = 0;
  size_exchange = 0;

  datamask_grow = datamask_comm = datamask_comm_vel = datamask_reverse =
    datamask_border = datamask_border_vel = datamask_exchange = EMPTY_MASK;

  k_count = DAT::tdual_int_1d("atom:k_count",1);
  atomKK = (AtomKokkos *) atom;
}

/* ---------------------------------------------------------------------- */

AtomVecKokkos::~AtomVecKokkos()
{
  // Kokkos already deallocated host memory

  ngrow = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG,int TRICLINIC,int DEFAULT>
struct AtomVecKokkos_PackComm {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr_randomread _x;
  typename AT::t_kkfloat_1d_4_randomread _mu;
  typename AT::t_kkfloat_1d_4_randomread _sp;
  typename AT::t_kkfloat_1d_randomread _dpdTheta,_uCond,_uMech,_uChem;
  typename AT::t_double_2d_lr_um _buf;
  typename AT::t_int_1d_const _list;
  double _xprd,_yprd,_zprd,_xy,_xz,_yz;
  double _pbc[6];
  uint64_t _datamask;

  AtomVecKokkos_PackComm(
    const AtomKokkos* atomKK,
    const typename DAT::tdual_double_2d_lr &buf,
    const typename DAT::tdual_int_1d &list,
    const double &xprd, const double &yprd, const double &zprd,
    const double &xy, const double &xz, const double &yz, const int* const pbc,
    const uint64_t &datamask):
    _x(atomKK->k_x.view<DeviceType>()),
    _mu(atomKK->k_mu.view<DeviceType>()),
    _sp(atomKK->k_sp.view<DeviceType>()),
    _dpdTheta(atomKK->k_dpdTheta.view<DeviceType>()),
    _uCond(atomKK->k_uCond.view<DeviceType>()),
    _uMech(atomKK->k_uMech.view<DeviceType>()),
    _uChem(atomKK->k_uChem.view<DeviceType>()),
    _list(list.view<DeviceType>()),
    _xprd(xprd),_yprd(yprd),_zprd(zprd),
    _xy(xy),_xz(xz),_yz(yz),_datamask(datamask) {
      const int size_forward = atomKK->avecKK->size_forward;
      const size_t maxsend = (buf.view<DeviceType>().extent(0)*buf.view<DeviceType>().extent(1))/size_forward;
      const size_t elements = size_forward;
      buffer_view<DeviceType>(_buf,buf,maxsend,elements);
      _pbc[0] = pbc[0]; _pbc[1] = pbc[1]; _pbc[2] = pbc[2];
      _pbc[3] = pbc[3]; _pbc[4] = pbc[4]; _pbc[5] = pbc[5];
    };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(i);
    int m = 0;
    if constexpr (PBC_FLAG == 0) {
      _buf(i,m++) = _x(j,0);
      _buf(i,m++) = _x(j,1);
      _buf(i,m++) = _x(j,2);
    } else {
      if (TRICLINIC == 0) {
        _buf(i,m++) = _x(j,0) + _pbc[0]*_xprd;
        _buf(i,m++) = _x(j,1) + _pbc[1]*_yprd;
        _buf(i,m++) = _x(j,2) + _pbc[2]*_zprd;
      } else {
        _buf(i,m++) = _x(j,0) + _pbc[0]*_xprd + _pbc[5]*_xy + _pbc[4]*_xz;
        _buf(i,m++) = _x(j,1) + _pbc[1]*_yprd + _pbc[3]*_yz;
        _buf(i,m++) = _x(j,2) + _pbc[2]*_zprd;
      }
    }

    if constexpr (!DEFAULT) {

      // DIPOLE package

      if (_datamask & MU_MASK) {
        _buf(i,m++) = _mu(j,0);
        _buf(i,m++) = _mu(j,1);
        _buf(i,m++) = _mu(j,2);
      }

      // SPIN package

      if (_datamask & SP_MASK) {
        _buf(i,m++) = _sp(j,0);
        _buf(i,m++) = _sp(j,1);
        _buf(i,m++) = _sp(j,2);
        _buf(i,m++) = _sp(j,3);
      }

      // DPD-REACT package

      if (_datamask & DPDTHETA_MASK) {
        _buf(i,m++) = _dpdTheta(j);
        _buf(i,m++) = _uCond(j);
        _buf(i,m++) = _uMech(j);
        _buf(i,m++) = _uChem(j);
      }
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecKokkos::pack_comm_kokkos(const int &n,
                                          const DAT::tdual_int_1d &list,
                                          const DAT::tdual_double_2d_lr &buf,
                                          const int &pbc_flag,
                                          const int* const pbc)
{
  // Check whether to always run forward communication on the host
  // Choose correct forward PackComm kernel

  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->sync(HostKK,datamask_comm);
    if (pbc_flag) {
      if (domain->triclinic) {
        if (comm_x_only) {
          struct AtomVecKokkos_PackComm<LMPHostType,1,1,1> f(atomKK,buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackComm<LMPHostType,1,1,0> f(atomKK,buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        }
      } else {
        if (comm_x_only) {
          struct AtomVecKokkos_PackComm<LMPHostType,1,0,1> f(atomKK,buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackComm<LMPHostType,1,0,0> f(atomKK,buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        }
      }
    } else {
      if (domain->triclinic) {
        if (comm_x_only) {
          struct AtomVecKokkos_PackComm<LMPHostType,0,1,1> f(atomKK,buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackComm<LMPHostType,0,1,0> f(atomKK,buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        }
      } else {
        if (comm_x_only) {
          struct AtomVecKokkos_PackComm<LMPHostType,0,0,1> f(atomKK,buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackComm<LMPHostType,0,0,0> f(atomKK,buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        }
      }
    }
  } else {
    atomKK->sync(Device,datamask_comm);
    if (pbc_flag) {
      if (domain->triclinic) {
        if (comm_x_only) {
          struct AtomVecKokkos_PackComm<LMPDeviceType,1,1,1> f(atomKK,buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackComm<LMPDeviceType,1,1,0> f(atomKK,buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        }
      } else {
        if (comm_x_only) {
          struct AtomVecKokkos_PackComm<LMPDeviceType,1,0,1> f(atomKK,buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackComm<LMPDeviceType,1,0,0> f(atomKK,buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        }
      }
    } else {
      if (domain->triclinic) {
        if (comm_x_only) {
          struct AtomVecKokkos_PackComm<LMPDeviceType,0,1,1> f(atomKK,buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackComm<LMPDeviceType,0,1,0> f(atomKK,buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        }
      } else {
        if (comm_x_only) {
          struct AtomVecKokkos_PackComm<LMPDeviceType,0,0,1> f(atomKK,buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackComm<LMPDeviceType,0,0,0> f(atomKK,buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        }
      }
    }
  }

  return n*size_forward;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG,int TRICLINIC,int DEFAULT>
struct AtomVecKokkos_PackCommSelf {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr _x;
  typename AT::t_kkfloat_1d_4 _mu;
  typename AT::t_kkfloat_1d_4 _sp;
  typename AT::t_kkfloat_1d _dpdTheta,_uCond,_uMech,_uChem;
  int _nfirst;
  typename AT::t_int_1d_const _list;
  double _xprd,_yprd,_zprd,_xy,_xz,_yz;
  double _pbc[6];
  uint64_t _datamask;

  AtomVecKokkos_PackCommSelf(
    const AtomKokkos* atomKK,
    const int &nfirst,
    const typename DAT::tdual_int_1d &list,
    const double &xprd, const double &yprd, const double &zprd,
    const double &xy, const double &xz, const double &yz, const int* const pbc,
    const uint64_t datamask):
    _x(atomKK->k_x.view<DeviceType>()),
    _mu(atomKK->k_mu.view<DeviceType>()),
    _sp(atomKK->k_sp.view<DeviceType>()),
    _dpdTheta(atomKK->k_dpdTheta.view<DeviceType>()),
    _uCond(atomKK->k_uCond.view<DeviceType>()),
    _uMech(atomKK->k_uMech.view<DeviceType>()),
    _uChem(atomKK->k_uChem.view<DeviceType>()),
    _nfirst(nfirst),_list(list.view<DeviceType>()),
    _xprd(xprd),_yprd(yprd),_zprd(zprd),
    _xy(xy),_xz(xz),_yz(yz),_datamask(datamask) {
      _pbc[0] = pbc[0]; _pbc[1] = pbc[1]; _pbc[2] = pbc[2];
      _pbc[3] = pbc[3]; _pbc[4] = pbc[4]; _pbc[5] = pbc[5];
  };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(i);
    if constexpr (PBC_FLAG == 0) {
      _x(i+_nfirst,0) = _x(j,0);
      _x(i+_nfirst,1) = _x(j,1);
      _x(i+_nfirst,2) = _x(j,2);
    } else {
      if (TRICLINIC == 0) {
        _x(i+_nfirst,0) = _x(j,0) + _pbc[0]*_xprd;
        _x(i+_nfirst,1) = _x(j,1) + _pbc[1]*_yprd;
        _x(i+_nfirst,2) = _x(j,2) + _pbc[2]*_zprd;
      } else {
        _x(i+_nfirst,0) = _x(j,0) + _pbc[0]*_xprd + _pbc[5]*_xy + _pbc[4]*_xz;
        _x(i+_nfirst,1) = _x(j,1) + _pbc[1]*_yprd + _pbc[3]*_yz;
        _x(i+_nfirst,2) = _x(j,2) + _pbc[2]*_zprd;
      }
    }

    if constexpr (!DEFAULT) {

      // DIPOLE package

      if (_datamask & MU_MASK) {
        _mu(i+_nfirst,0) = _mu(j,0);
        _mu(i+_nfirst,1) = _mu(j,1);
        _mu(i+_nfirst,2) = _mu(j,2);
      }

      // SPIN package

      if (_datamask & SP_MASK) {
        _sp(i+_nfirst,0) = _sp(j,0);
        _sp(i+_nfirst,1) = _sp(j,1);
        _sp(i+_nfirst,2) = _sp(j,2);
        _sp(i+_nfirst,3) = _sp(j,3);
      }

      // DPD-REACT package

      if (_datamask & DPDTHETA_MASK) {
        _dpdTheta(i+_nfirst) = _dpdTheta(j);
        _uCond(i+_nfirst) = _uCond(j);
        _uMech(i+_nfirst) = _uMech(j);
        _uChem(i+_nfirst) = _uChem(j);
      }
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecKokkos::pack_comm_self(const int &n, const DAT::tdual_int_1d &list,
                                        const int nfirst, const int &pbc_flag, const int* const pbc) {
  // Check whether to always run forward communication on the host
  // Choose correct forward PackComm kernel

  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->sync(HostKK,datamask_comm);
    if (pbc_flag) {
      if (domain->triclinic) {
        if (comm_x_only) {
          struct AtomVecKokkos_PackCommSelf<LMPHostType,1,1,1> f(atomKK,nfirst,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackCommSelf<LMPHostType,1,1,0> f(atomKK,nfirst,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        }
      } else {
        if (comm_x_only) {
          struct AtomVecKokkos_PackCommSelf<LMPHostType,1,0,1> f(atomKK,nfirst,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackCommSelf<LMPHostType,1,0,0> f(atomKK,nfirst,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        }
      }
    } else {
      if (domain->triclinic) {
        if (comm_x_only) {
          struct AtomVecKokkos_PackCommSelf<LMPHostType,0,1,1> f(atomKK,nfirst,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackCommSelf<LMPHostType,0,1,0> f(atomKK,nfirst,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        }
      } else {
        if (comm_x_only) {
          struct AtomVecKokkos_PackCommSelf<LMPHostType,0,0,1> f(atomKK,nfirst,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackCommSelf<LMPHostType,0,0,0> f(atomKK,nfirst,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        }
      }
    }
    atomKK->modified(HostKK,datamask_comm);
  } else {
    atomKK->sync(Device,datamask_comm);
    if (pbc_flag) {
      if (domain->triclinic) {
        if (comm_x_only) {
          struct AtomVecKokkos_PackCommSelf<LMPDeviceType,1,1,1> f(atomKK,nfirst,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackCommSelf<LMPDeviceType,1,1,0> f(atomKK,nfirst,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        }
      } else {
        if (comm_x_only) {
          struct AtomVecKokkos_PackCommSelf<LMPDeviceType,1,0,1> f(atomKK,nfirst,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackCommSelf<LMPDeviceType,1,0,0> f(atomKK,nfirst,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        }
      }
    } else {
      if (domain->triclinic) {
        if (comm_x_only) {
          struct AtomVecKokkos_PackCommSelf<LMPDeviceType,0,1,1> f(atomKK,nfirst,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackCommSelf<LMPDeviceType,0,1,0> f(atomKK,nfirst,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        }
      } else {
        if (comm_x_only) {
          struct AtomVecKokkos_PackCommSelf<LMPDeviceType,0,0,1> f(atomKK,nfirst,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackCommSelf<LMPDeviceType,0,0,0> f(atomKK,nfirst,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,datamask_comm);
          Kokkos::parallel_for(n,f);
        }
      }
    }
    atomKK->modified(Device,datamask_comm);
  }

  return n*size_forward;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int TRICLINIC,int DEFAULT>
struct AtomVecKokkos_PackCommSelfFused {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr _x;
  typename AT::t_kkfloat_1d_4 _mu;
  typename AT::t_kkfloat_1d_4 _sp;
  typename AT::t_kkfloat_1d _dpdTheta,_uCond,_uMech,_uChem;
  typename AT::t_int_2d_lr_const _list;
  typename AT::t_int_2d_const _pbc;
  typename AT::t_int_1d_const _pbc_flag;
  typename AT::t_int_1d_const _firstrecv;
  typename AT::t_int_1d_const _sendnum_scan;
  typename AT::t_int_1d_const _g2l;
  double _xprd,_yprd,_zprd,_xy,_xz,_yz;
  uint64_t _datamask;

  AtomVecKokkos_PackCommSelfFused(
      const AtomKokkos* atomKK,
      const typename DAT::tdual_int_2d_lr &list,
      const typename DAT::tdual_int_2d &pbc,
      const typename DAT::tdual_int_1d &pbc_flag,
      const typename DAT::tdual_int_1d &firstrecv,
      const typename DAT::tdual_int_1d &sendnum_scan,
      const typename DAT::tdual_int_1d &g2l,
      const double &xprd, const double &yprd, const double &zprd,
      const double &xy, const double &xz, const double &yz,
      const uint64_t datamask):
      _x(atomKK->k_x.view<DeviceType>()),
      _mu(atomKK->k_mu.view<DeviceType>()),
      _sp(atomKK->k_sp.view<DeviceType>()),
      _dpdTheta(atomKK->k_dpdTheta.view<DeviceType>()),
      _uCond(atomKK->k_uCond.view<DeviceType>()),
      _uMech(atomKK->k_uMech.view<DeviceType>()),
      _uChem(atomKK->k_uChem.view<DeviceType>()),
      _list(list.view<DeviceType>()),
      _pbc(pbc.view<DeviceType>()),
      _pbc_flag(pbc_flag.view<DeviceType>()),
      _firstrecv(firstrecv.view<DeviceType>()),
      _sendnum_scan(sendnum_scan.view<DeviceType>()),
      _g2l(g2l.view<DeviceType>()),
      _xprd(xprd),_yprd(yprd),_zprd(zprd),
      _xy(xy),_xz(xz),_yz(yz),_datamask(datamask) {};

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

    if (_pbc_flag(ii) == 0) {
      _x(i+_nfirst,0) = _x(j,0);
      _x(i+_nfirst,1) = _x(j,1);
      _x(i+_nfirst,2) = _x(j,2);
    } else {
      if (TRICLINIC == 0) {
        _x(i+_nfirst,0) = _x(j,0) + _pbc(ii,0)*_xprd;
        _x(i+_nfirst,1) = _x(j,1) + _pbc(ii,1)*_yprd;
        _x(i+_nfirst,2) = _x(j,2) + _pbc(ii,2)*_zprd;
      } else {
        _x(i+_nfirst,0) = _x(j,0) + _pbc(ii,0)*_xprd + _pbc(ii,5)*_xy + _pbc(ii,4)*_xz;
        _x(i+_nfirst,1) = _x(j,1) + _pbc(ii,1)*_yprd + _pbc(ii,3)*_yz;
        _x(i+_nfirst,2) = _x(j,2) + _pbc(ii,2)*_zprd;
      }
    }

    if constexpr (!DEFAULT) {

      // DIPOLE package

      if (_datamask & MU_MASK) {
        _mu(i+_nfirst,0) = _mu(j,0);
        _mu(i+_nfirst,1) = _mu(j,1);
        _mu(i+_nfirst,2) = _mu(j,2);
      }

      // SPIN package

      if (_datamask & SP_MASK) {
        _sp(i+_nfirst,0) = _sp(j,0);
        _sp(i+_nfirst,1) = _sp(j,1);
        _sp(i+_nfirst,2) = _sp(j,2);
        _sp(i+_nfirst,3) = _sp(j,3);
      }

      // DPD-REACT package

      if (_datamask & DPDTHETA_MASK) {
        _dpdTheta(i+_nfirst) = _dpdTheta(j);
        _uCond(i+_nfirst) = _uCond(j);
        _uMech(i+_nfirst) = _uMech(j);
        _uChem(i+_nfirst) = _uChem(j);
      }
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecKokkos::pack_comm_self_fused(const int &n, const DAT::tdual_int_2d_lr &list, const DAT::tdual_int_1d &sendnum_scan,
                                         const DAT::tdual_int_1d &firstrecv, const DAT::tdual_int_1d &pbc_flag, const DAT::tdual_int_2d &pbc,
                                         const DAT::tdual_int_1d &g2l) {
  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->sync(HostKK,datamask_comm);
    if (domain->triclinic) {
      if (comm_x_only) {
        struct AtomVecKokkos_PackCommSelfFused<LMPHostType,1,1> f(atomKK,list,pbc,pbc_flag,firstrecv,sendnum_scan,g2l,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,datamask_comm);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecKokkos_PackCommSelfFused<LMPHostType,1,0> f(atomKK,list,pbc,pbc_flag,firstrecv,sendnum_scan,g2l,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,datamask_comm);
        Kokkos::parallel_for(n,f);
      }
    } else {
      if (comm_x_only) {
        struct AtomVecKokkos_PackCommSelfFused<LMPHostType,0,1> f(atomKK,list,pbc,pbc_flag,firstrecv,sendnum_scan,g2l,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,datamask_comm);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecKokkos_PackCommSelfFused<LMPHostType,0,0> f(atomKK,list,pbc,pbc_flag,firstrecv,sendnum_scan,g2l,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,datamask_comm);
        Kokkos::parallel_for(n,f);
      }
    }
    atomKK->modified(HostKK,datamask_comm);
  } else {
    atomKK->sync(Device,datamask_comm);
    if (domain->triclinic) {
      if (comm_x_only) {
        struct AtomVecKokkos_PackCommSelfFused<LMPDeviceType,1,1> f(atomKK,list,pbc,pbc_flag,firstrecv,sendnum_scan,g2l,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,datamask_comm);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecKokkos_PackCommSelfFused<LMPDeviceType,1,0> f(atomKK,list,pbc,pbc_flag,firstrecv,sendnum_scan,g2l,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,datamask_comm);
        Kokkos::parallel_for(n,f);
      }
    } else {
      if (comm_x_only) {
        struct AtomVecKokkos_PackCommSelfFused<LMPDeviceType,0,1> f(atomKK,list,pbc,pbc_flag,firstrecv,sendnum_scan,g2l,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,datamask_comm);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecKokkos_PackCommSelfFused<LMPDeviceType,0,0> f(atomKK,list,pbc,pbc_flag,firstrecv,sendnum_scan,g2l,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,datamask_comm);
        Kokkos::parallel_for(n,f);
      }
    }
    atomKK->modified(Device,datamask_comm);
  }

  return n*size_forward;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int DEFAULT>
struct AtomVecKokkos_UnpackComm {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr _x;
  typename AT::t_kkfloat_1d_4 _mu;
  typename AT::t_kkfloat_1d_4 _sp;
  typename AT::t_kkfloat_1d _dpdTheta,_uCond,_uMech,_uChem;
  typename AT::t_double_2d_lr_const _buf;
  int _first;
  uint64_t _datamask;

  AtomVecKokkos_UnpackComm(
    const AtomKokkos* atomKK,
    const typename DAT::tdual_double_2d_lr &buf,
    const int &first, const uint64_t &datamask):
      _x(atomKK->k_x.view<DeviceType>()),
      _mu(atomKK->k_mu.view<DeviceType>()),
      _sp(atomKK->k_sp.view<DeviceType>()),
      _dpdTheta(atomKK->k_dpdTheta.view<DeviceType>()),
      _uCond(atomKK->k_uCond.view<DeviceType>()),
      _uMech(atomKK->k_uMech.view<DeviceType>()),
      _uChem(atomKK->k_uChem.view<DeviceType>()),
      _first(first),_datamask(datamask) {
        const int size_forward = atomKK->avecKK->size_forward;
        const size_t maxsend = (buf.view<DeviceType>().extent(0)*buf.view<DeviceType>().extent(1))/size_forward;
        const size_t elements = size_forward;
        buffer_view<DeviceType>(_buf,buf,maxsend,elements);
      };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    int m = 0;
    _x(i+_first,0) = _buf(i,m++);
    _x(i+_first,1) = _buf(i,m++);
    _x(i+_first,2) = _buf(i,m++);

    if constexpr (!DEFAULT) {

      // DIPOLE package

      if (_datamask & MU_MASK) {
        _mu(i+_first,0) = _buf(i,m++);
        _mu(i+_first,1) = _buf(i,m++);
        _mu(i+_first,2) = _buf(i,m++);
      }

      // SPIN package

      if (_datamask & SP_MASK) {
        _sp(i+_first,0) = _buf(i,m++);
        _sp(i+_first,1) = _buf(i,m++);
        _sp(i+_first,2) = _buf(i,m++);
        _sp(i+_first,3) = _buf(i,m++);
      }

      // DPD-REACT package

      if (_datamask & DPDTHETA_MASK) {
        _dpdTheta(i+_first) = _buf(i,m++);
        _uCond(i+_first) = _buf(i,m++);
        _uMech(i+_first) = _buf(i,m++);
        _uChem(i+_first) = _buf(i,m++);
      }
    }
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecKokkos::unpack_comm_kokkos(const int &n, const int &first,
    const DAT::tdual_double_2d_lr &buf) {
  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->sync(HostKK,datamask_comm);
    if (comm_x_only) {
      struct AtomVecKokkos_UnpackComm<LMPHostType,1> f(atomKK,buf,first,datamask_comm);
      Kokkos::parallel_for(n,f);
    } else {
      struct AtomVecKokkos_UnpackComm<LMPHostType,0> f(atomKK,buf,first,datamask_comm);
      Kokkos::parallel_for(n,f);
    }
    atomKK->modified(HostKK,datamask_comm);
  } else {
    atomKK->sync(Device,datamask_comm);
    if (comm_x_only) {
      struct AtomVecKokkos_UnpackComm<LMPDeviceType,1> f(atomKK,buf,first,datamask_comm);
      Kokkos::parallel_for(n,f);
    } else {
      struct AtomVecKokkos_UnpackComm<LMPDeviceType,0> f(atomKK,buf,first,datamask_comm);
      Kokkos::parallel_for(n,f);
    }
    atomKK->modified(Device,datamask_comm);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG,int TRICLINIC,int DEFORM_VREMAP>
struct AtomVecKokkos_PackCommVel {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr_randomread _x;
  typename AT::t_int_1d_randomread _mask;
  typename AT::t_kkfloat_1d_3_randomread _v;
  typename AT::t_kkfloat_1d_4_randomread _mu;
  typename AT::t_kkfloat_1d_4_randomread _sp;
  typename AT::t_kkfloat_1d_3_randomread _omega;
  typename AT::t_kkfloat_1d_randomread _dpdTheta,_uCond,_uMech,_uChem;
  typename AT::t_double_2d_lr_um _buf;
  typename AT::t_int_1d_const _list;
  double _xprd,_yprd,_zprd,_xy,_xz,_yz;
  double _pbc[6];
  double _h_rate[6];
  const int _deform_vremap;
  uint64_t _datamask;

  AtomVecKokkos_PackCommVel(
    const AtomKokkos* atomKK,
    const typename DAT::tdual_double_2d_lr &buf,
    const typename DAT::tdual_int_1d &list,
    const double &xprd, const double &yprd, const double &zprd,
    const double &xy, const double &xz, const double &yz, const int* const pbc,
    const double * const h_rate,
    const int &deform_vremap,
    const uint64_t &datamask):
    _x(atomKK->k_x.view<DeviceType>()),
    _mask(atomKK->k_mask.view<DeviceType>()),
    _v(atomKK->k_v.view<DeviceType>()),
    _omega(atomKK->k_omega.view<DeviceType>()),
    _dpdTheta(atomKK->k_dpdTheta.view<DeviceType>()),
    _uCond(atomKK->k_uCond.view<DeviceType>()),
    _uMech(atomKK->k_uMech.view<DeviceType>()),
    _uChem(atomKK->k_uChem.view<DeviceType>()),
    _list(list.view<DeviceType>()),
    _xprd(xprd),_yprd(yprd),_zprd(zprd),
    _xy(xy),_xz(xz),_yz(yz),
    _deform_vremap(deform_vremap),
    _datamask(datamask)
  {
    const size_t elements = atomKK->avecKK->size_forward + atomKK->avecKK->size_velocity;
    const int maxsend = (buf.template view<DeviceType>().extent(0)*buf.template view<DeviceType>().extent(1))/elements;
    buffer_view<DeviceType>(_buf,buf,maxsend,elements);
    _pbc[0] = pbc[0]; _pbc[1] = pbc[1]; _pbc[2] = pbc[2];
    _pbc[3] = pbc[3]; _pbc[4] = pbc[4]; _pbc[5] = pbc[5];
    _h_rate[0] = h_rate[0]; _h_rate[1] = h_rate[1]; _h_rate[2] = h_rate[2];
    _h_rate[3] = h_rate[3]; _h_rate[4] = h_rate[4]; _h_rate[5] = h_rate[5];
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    int m = 0;
    const int j = _list(i);
    if constexpr (PBC_FLAG == 0) {
      _buf(i,m++) = _x(j,0);
      _buf(i,m++) = _x(j,1);
      _buf(i,m++) = _x(j,2);
      _buf(i,m++) = _v(j,0);
      _buf(i,m++) = _v(j,1);
      _buf(i,m++) = _v(j,2);
    } else {
      if (TRICLINIC == 0) {
        _buf(i,m++) = _x(j,0) + _pbc[0]*_xprd;
        _buf(i,m++) = _x(j,1) + _pbc[1]*_yprd;
        _buf(i,m++) = _x(j,2) + _pbc[2]*_zprd;
             } else {
        _buf(i,m++) = _x(j,0) + _pbc[0]*_xprd + _pbc[5]*_xy + _pbc[4]*_xz;
        _buf(i,m++) = _x(j,1) + _pbc[1]*_yprd + _pbc[3]*_yz;
        _buf(i,m++) = _x(j,2) + _pbc[2]*_zprd;
      }

      if constexpr (DEFORM_VREMAP == 0) {
        _buf(i,m++) = _v(j,0);
        _buf(i,m++) = _v(j,1);
        _buf(i,m++) = _v(j,2);
      } else {
        if (_mask(i) & _deform_vremap) {
          _buf(i,m++) = _v(j,0) + _pbc[0]*_h_rate[0] + _pbc[5]*_h_rate[5] + _pbc[4]*_h_rate[4];
          _buf(i,m++) = _v(j,1) + _pbc[1]*_h_rate[1] + _pbc[3]*_h_rate[3];
          _buf(i,m++) = _v(j,2) + _pbc[2]*_h_rate[2];
        } else {
          _buf(i,m++) = _v(j,0);
          _buf(i,m++) = _v(j,1);
          _buf(i,m++) = _v(j,2);
        }
      }
    }

    // DIPOLE package

    if (_datamask & MU_MASK) {
      _buf(i,m++) = _mu(j,0);
      _buf(i,m++) = _mu(j,1);
      _buf(i,m++) = _mu(j,2);
    }

    // SPIN package

    if (_datamask & SP_MASK) {
      _buf(i,m++) = _sp(j,0);
      _buf(i,m++) = _sp(j,1);
      _buf(i,m++) = _sp(j,2);
      _buf(i,m++) = _sp(j,3);
    }

    // SPHERE package

    if (_datamask & OMEGA_MASK) {
      _buf(i,m++) = _omega(j,0);
      _buf(i,m++) = _omega(j,1);
      _buf(i,m++) = _omega(j,2);
    }

      // DPD-REACT package

    if (_datamask & DPDTHETA_MASK) {
      _buf(i,m++) = _dpdTheta(j);
      _buf(i,m++) = _uCond(j);
      _buf(i,m++) = _uMech(j);
      _buf(i,m++) = _uChem(j);
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecKokkos::pack_comm_vel_kokkos(
  const int &n,
  const DAT::tdual_int_1d &list,
  const DAT::tdual_double_2d_lr &buf,
  const int &pbc_flag,
  const int* const pbc)
{
  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->sync(HostKK,datamask_comm_vel);
    if (pbc_flag) {
      if (deform_vremap) {
        if (domain->triclinic) {
          struct AtomVecKokkos_PackCommVel<LMPHostType,1,1,1> f(
            atomKK,
            buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap,
            datamask_comm_vel);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackCommVel<LMPHostType,1,0,1> f(
            atomKK,
            buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap,
            datamask_comm_vel);
          Kokkos::parallel_for(n,f);
        }
      } else {
        if (domain->triclinic) {
          struct AtomVecKokkos_PackCommVel<LMPHostType,1,1,0> f(
            atomKK,
            buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap,
            datamask_comm_vel);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackCommVel<LMPHostType,1,0,0> f(
            atomKK,
            buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap,
            datamask_comm_vel);
          Kokkos::parallel_for(n,f);
        }
      }
    } else {
      if (domain->triclinic) {
        struct AtomVecKokkos_PackCommVel<LMPHostType,0,1,0> f(
          atomKK,
          buf,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap,
          datamask_comm_vel);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecKokkos_PackCommVel<LMPHostType,0,0,0> f(
          atomKK,
          buf,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap,
          datamask_comm_vel);
        Kokkos::parallel_for(n,f);
      }
    }
  } else {
    atomKK->sync(Device,datamask_comm_vel);
    if (pbc_flag) {
      if (deform_vremap) {
        if (domain->triclinic) {
          struct AtomVecKokkos_PackCommVel<LMPDeviceType,1,1,1> f(
            atomKK,
            buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap,
            datamask_comm_vel);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackCommVel<LMPDeviceType,1,0,1> f(
            atomKK,
            buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap,
            datamask_comm_vel);
          Kokkos::parallel_for(n,f);
        }
      } else {
        if (domain->triclinic) {
          struct AtomVecKokkos_PackCommVel<LMPDeviceType,1,1,0> f(
            atomKK,
            buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap,
            datamask_comm_vel);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackCommVel<LMPDeviceType,1,0,0> f(
            atomKK,
            buf,list,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap,
            datamask_comm_vel);
          Kokkos::parallel_for(n,f);
        }
      }
    } else {
      if (domain->triclinic) {
        struct AtomVecKokkos_PackCommVel<LMPDeviceType,0,1,0> f(
          atomKK,
          buf,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap,
          datamask_comm_vel);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecKokkos_PackCommVel<LMPDeviceType,0,0,0> f(
          atomKK,
          buf,list,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap,
          datamask_comm_vel);
        Kokkos::parallel_for(n,f);
      }
    }
  }

  return n*(size_forward + size_velocity);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int DEFAULT>
struct AtomVecKokkos_UnpackCommVel {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr _x;
  typename AT::t_kkfloat_1d_3 _v;
  typename AT::t_kkfloat_1d_4 _mu;
  typename AT::t_kkfloat_1d_4 _sp;
  typename AT::t_kkfloat_1d_3 _omega;
  typename AT::t_kkfloat_1d _dpdTheta,_uCond,_uMech,_uChem;
  typename AT::t_double_2d_lr_const _buf;
  int _first;
  uint64_t _datamask;

  AtomVecKokkos_UnpackCommVel(
    const AtomKokkos* atomKK,
    const typename DAT::tdual_double_2d_lr &buf,
    const int &first, const uint64_t &datamask):
    _x(atomKK->k_x.view<DeviceType>()),
    _v(atomKK->k_v.view<DeviceType>()),
    _mu(atomKK->k_mu.view<DeviceType>()),
    _sp(atomKK->k_sp.view<DeviceType>()),
    _omega(atomKK->k_omega.view<DeviceType>()),
    _dpdTheta(atomKK->k_dpdTheta.view<DeviceType>()),
    _uCond(atomKK->k_uCond.view<DeviceType>()),
    _uMech(atomKK->k_uMech.view<DeviceType>()),
    _uChem(atomKK->k_uChem.view<DeviceType>()),
    _first(first),_datamask(datamask)
  {
    const size_t elements = atomKK->avecKK->size_forward + atomKK->avecKK->size_velocity;
    const int maxsend = (buf.template view<DeviceType>().extent(0)*buf.template view<DeviceType>().extent(1))/elements;
    buffer_view<DeviceType>(_buf,buf,maxsend,elements);
  };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    int m = 0;
    _x(i+_first,0) = _buf(i,m++);
    _x(i+_first,1) = _buf(i,m++);
    _x(i+_first,2) = _buf(i,m++);
    _v(i+_first,0) = _buf(i,m++);
    _v(i+_first,1) = _buf(i,m++);
    _v(i+_first,2) = _buf(i,m++);

    if constexpr (!DEFAULT) {

      // DIPOLE package

      if (_datamask & MU_MASK) {
        _mu(i+_first,0) = _buf(i,m++);
        _mu(i+_first,1) = _buf(i,m++);
        _mu(i+_first,2) = _buf(i,m++);
      }

      // SPIN package

      if (_datamask & SP_MASK) {
        _sp(i+_first,0) = _buf(i,m++);
        _sp(i+_first,1) = _buf(i,m++);
        _sp(i+_first,2) = _buf(i,m++);
        _sp(i+_first,3) = _buf(i,m++);
      }

      // SPHERE package

      if (_datamask & OMEGA_MASK) {
        _omega(i+_first,0) = _buf(i,m++);
        _omega(i+_first,1) = _buf(i,m++);
        _omega(i+_first,2) = _buf(i,m++);
      }

      // DPD-REACT package

      if (_datamask & DPDTHETA_MASK) {
        _dpdTheta(i+_first) = _buf(i,m++);
        _uCond(i+_first) = _buf(i,m++);
        _uMech(i+_first) = _buf(i,m++);
        _uChem(i+_first) = _buf(i,m++);
      }
    }
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecKokkos::unpack_comm_vel_kokkos(const int &n, const int &first,
    const DAT::tdual_double_2d_lr &buf) {

  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->sync(HostKK,datamask_comm_vel);
    if (!ncomm_vel) {
      struct AtomVecKokkos_UnpackCommVel<LMPHostType,1> f(atomKK,buf,first,datamask_comm_vel);
      Kokkos::parallel_for(n,f);
    } else {
      struct AtomVecKokkos_UnpackCommVel<LMPHostType,0> f(atomKK,buf,first,datamask_comm_vel);
      Kokkos::parallel_for(n,f);
    }
    atomKK->modified(HostKK,datamask_comm_vel);
  } else {
    atomKK->sync(Device,datamask_comm_vel);
    if (!ncomm_vel) {
      struct AtomVecKokkos_UnpackCommVel<LMPDeviceType,1> f(atomKK,buf,first,datamask_comm_vel);
      Kokkos::parallel_for(n,f);
    } else {
      struct AtomVecKokkos_UnpackCommVel<LMPDeviceType,0> f(atomKK,buf,first,datamask_comm_vel);
      Kokkos::parallel_for(n,f);
    }
    atomKK->modified(Device,datamask_comm_vel);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int DEFAULT>
struct AtomVecKokkos_PackReverse {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkacc_1d_3_randomread _f,_fm,_fm_long;
  typename AT::t_kkacc_1d_3_randomread _torque;
  typename AT::t_double_2d_lr _buf;
  int _first;
  uint64_t _datamask;

  AtomVecKokkos_PackReverse(
    const AtomKokkos* atomKK,
    const typename DAT::tdual_double_2d_lr &buf,
    const int &first, const uint64_t &datamask):
      _f(atomKK->k_f.view<DeviceType>()),
      _torque(atomKK->k_torque.view<DeviceType>()),
      _fm(atomKK->k_fm.view<DeviceType>()),
      _fm_long(atomKK->k_fm_long.view<DeviceType>()),
      _first(first),_datamask(datamask) {
        const size_t elements = atomKK->avecKK->size_reverse;
        const size_t maxsend = (buf.view<DeviceType>().extent(0)*buf.view<DeviceType>().extent(1))/elements;
        buffer_view<DeviceType>(_buf,buf,maxsend,elements);
      };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    int m = 0;
    _buf(i,m++) = _f(i+_first,0);
    _buf(i,m++) = _f(i+_first,1);
    _buf(i,m++) = _f(i+_first,2);

    if constexpr (!DEFAULT) {

      // DIPLE package

      if (_datamask & TORQUE_MASK) {
        _buf(i,m++) = _torque(i+_first,0);
        _buf(i,m++) = _torque(i+_first,1);
        _buf(i,m++) = _torque(i+_first,2);
      }

      // SPIN package

      if (_datamask & FM_MASK) {
        _buf(i,m++) = _fm(i+_first,0);
        _buf(i,m++) = _fm(i+_first,1);
        _buf(i,m++) = _fm(i+_first,2);

        _buf(i,m++) = _fm_long(i+_first,0);
        _buf(i,m++) = _fm_long(i+_first,1);
        _buf(i,m++) = _fm_long(i+_first,2);
      }
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecKokkos::pack_reverse_kokkos(const int &n, const int &first,
    const DAT::tdual_double_2d_lr &buf) {
  if (lmp->kokkos->reverse_comm_on_host) {
    atomKK->sync(HostKK,datamask_reverse);
    if (comm_f_only) {
      struct AtomVecKokkos_PackReverse<LMPHostType,1> f(atomKK,buf,first,datamask_reverse);
      Kokkos::parallel_for(n,f);
    } else {
      struct AtomVecKokkos_PackReverse<LMPHostType,0> f(atomKK,buf,first,datamask_reverse);
      Kokkos::parallel_for(n,f);
    }
  } else {
    atomKK->sync(Device,datamask_reverse);
    if (comm_f_only) {
      struct AtomVecKokkos_PackReverse<LMPDeviceType,1> f(atomKK,buf,first,datamask_reverse);
      Kokkos::parallel_for(n,f);
    } else {
      struct AtomVecKokkos_PackReverse<LMPDeviceType,0> f(atomKK,buf,first,datamask_reverse);
      Kokkos::parallel_for(n,f);
    }
  }

  return n*size_reverse;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int DEFAULT>
struct AtomVecKokkos_UnPackReverseSelf {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkacc_1d_3 _f,_fm,_fm_long;
  typename AT::t_kkacc_1d_3 _torque;
  typename AT::t_int_1d_const _list;
  int _nfirst;
  uint64_t _datamask;

  AtomVecKokkos_UnPackReverseSelf(
    const AtomKokkos* atomKK,
    const int &nfirst,
    const typename DAT::tdual_int_1d &list,
    const uint64_t &datamask):
      _f(atomKK->k_f.view<DeviceType>()),
      _torque(atomKK->k_torque.view<DeviceType>()),
      _fm(atomKK->k_fm.view<DeviceType>()),
      _fm_long(atomKK->k_fm_long.view<DeviceType>()),
      _nfirst(nfirst),_list(list.view<DeviceType>()),
      _datamask(datamask) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(i);
    _f(j,0) += _f(i+_nfirst,0);
    _f(j,1) += _f(i+_nfirst,1);
    _f(j,2) += _f(i+_nfirst,2);

    if constexpr (!DEFAULT) {

      // DIPOLE package

      if (_datamask & TORQUE_MASK) {
        _torque(j,0) += _torque(i+_nfirst,0);
        _torque(j,1) += _torque(i+_nfirst,1);
        _torque(j,2) += _torque(i+_nfirst,2);
      }

      // SPIN package

      if (_datamask & FM_MASK) {
        _fm(j,0) += _fm(i+_nfirst,0);
        _fm(j,1) += _fm(i+_nfirst,1);
        _fm(j,2) += _fm(i+_nfirst,2);

        _fm_long(j,0) += _fm_long(i+_nfirst,0);
        _fm_long(j,1) += _fm_long(i+_nfirst,1);
        _fm_long(j,2) += _fm_long(i+_nfirst,2);
      }
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecKokkos::pack_reverse_self(const int &n, const DAT::tdual_int_1d &list,
                                     const int nfirst) {
  if (lmp->kokkos->reverse_comm_on_host) {
    atomKK->sync(HostKK,datamask_reverse);
    if (comm_f_only) {
      struct AtomVecKokkos_UnPackReverseSelf<LMPHostType,1> f(atomKK,nfirst,list,datamask_reverse);
      Kokkos::parallel_for(n,f);
    } else {
      struct AtomVecKokkos_UnPackReverseSelf<LMPHostType,0> f(atomKK,nfirst,list,datamask_reverse);
      Kokkos::parallel_for(n,f);
    }
    atomKK->modified(HostKK,datamask_reverse);
  } else {
    atomKK->sync(Device,datamask_reverse);
    if (comm_f_only) {
      struct AtomVecKokkos_UnPackReverseSelf<LMPDeviceType,1> f(atomKK,nfirst,list,datamask_reverse);
      Kokkos::parallel_for(n,f);
    } else {
      struct AtomVecKokkos_UnPackReverseSelf<LMPDeviceType,0> f(atomKK,nfirst,list,datamask_reverse);
      Kokkos::parallel_for(n,f);
    }
    atomKK->modified(Device,datamask_reverse);
  }

  return n*size_reverse;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int DEFAULT>
struct AtomVecKokkos_UnPackReverse {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkacc_1d_3 _f,_fm,_fm_long;
  typename AT::t_kkacc_1d_3 _torque;
  typename AT::t_double_2d_lr_const _buf;
  typename AT::t_int_1d_const _list;
  uint64_t _datamask;

  AtomVecKokkos_UnPackReverse(
    const AtomKokkos* atomKK,
    const typename DAT::tdual_double_2d_lr &buf,
    const typename DAT::tdual_int_1d &list,
    const uint64_t datamask):
      _f(atomKK->k_f.view<DeviceType>()),
      _torque(atomKK->k_torque.view<DeviceType>()),
      _fm(atomKK->k_fm.view<DeviceType>()),
      _fm_long(atomKK->k_fm_long.view<DeviceType>()),
      _list(list.view<DeviceType>()),
      _datamask(datamask) {
        const size_t elements = atomKK->avecKK->size_reverse;
        const size_t maxsend = (buf.view<DeviceType>().extent(0)*buf.view<DeviceType>().extent(1))/elements;
        buffer_view<DeviceType>(_buf,buf,maxsend,elements);
      };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    int m = 0;
    const int j = _list(i);
    _f(j,0) += _buf(i,m++);
    _f(j,1) += _buf(i,m++);
    _f(j,2) += _buf(i,m++);

    if constexpr (!DEFAULT) {

      // DIPOLE package

      if (_datamask & TORQUE_MASK) {
        _torque(j,0) += _buf(i,m++);
        _torque(j,1) += _buf(i,m++);
        _torque(j,2) += _buf(i,m++);
      }

      // SPIN package

      if (_datamask & FM_MASK) {
        _fm(j,0) += _buf(i,m++);
        _fm(j,1) += _buf(i,m++);
        _fm(j,2) += _buf(i,m++);

        _fm_long(j,0) += _buf(i,m++);
        _fm_long(j,1) += _buf(i,m++);
        _fm_long(j,2) += _buf(i,m++);
      }
    }
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecKokkos::unpack_reverse_kokkos(const int &n,
                                          const DAT::tdual_int_1d &list,
                                          const DAT::tdual_double_2d_lr &buf)
{
  // Check whether to always run reverse communication on the host
  // Choose correct reverse UnPackReverse kernel

  if (lmp->kokkos->reverse_comm_on_host) {
    atomKK->sync(HostKK,datamask_reverse);
    if (comm_f_only) {
      struct AtomVecKokkos_UnPackReverse<LMPHostType,1> f(atomKK,buf,list,datamask_reverse);
      Kokkos::parallel_for(n,f);
    } else {
      struct AtomVecKokkos_UnPackReverse<LMPHostType,0> f(atomKK,buf,list,datamask_reverse);
      Kokkos::parallel_for(n,f);
    }
    atomKK->modified(HostKK,datamask_reverse);
  } else {
    atomKK->sync(Device,datamask_reverse);
    if (comm_f_only) {
      struct AtomVecKokkos_UnPackReverse<LMPDeviceType,1> f(atomKK,buf,list,datamask_reverse);
      Kokkos::parallel_for(n,f);
    } else {
      struct AtomVecKokkos_UnPackReverse<LMPDeviceType,0> f(atomKK,buf,list,datamask_reverse);
      Kokkos::parallel_for(n,f);
    }
    atomKK->modified(Device,datamask_reverse);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG,int DEFAULT>
struct AtomVecKokkos_PackBorder {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_double_2d_lr _buf;
  const typename AT::t_int_1d_const _list;
  const typename AT::t_kkfloat_1d_3_lr_randomread _x;
  const typename AT::t_tagint_1d_randomread _tag;
  const typename AT::t_int_1d_randomread _type;
  const typename AT::t_int_1d_randomread _mask;
  const typename AT::t_tagint_1d_randomread _molecule;
  const typename AT::t_kkfloat_1d_randomread _q;
  const typename AT::t_kkfloat_1d_4_randomread _mu;
  const typename AT::t_kkfloat_1d_4_randomread _sp;
  typename AT::t_kkfloat_1d_randomread _radius,_rmass;
  typename AT::t_kkfloat_1d_randomread _dpdTheta,_uCond,_uMech,_uChem,_uCG,_uCGnew;
  double _dx,_dy,_dz;
  uint64_t _datamask;

  AtomVecKokkos_PackBorder(
    const AtomKokkos* atomKK,
    const typename AT::t_double_2d_lr &buf,
    const typename AT::t_int_1d_const &list,
    const double &dx, const double &dy, const double &dz,
    const uint64_t &datamask):
      _buf(buf),_list(list),
      _x(atomKK->k_x.view<DeviceType>()),
      _tag(atomKK->k_tag.view<DeviceType>()),
      _type(atomKK->k_type.view<DeviceType>()),
      _mask(atomKK->k_mask.view<DeviceType>()),
      _molecule(atomKK->k_molecule.view<DeviceType>()),
      _q(atomKK->k_q.view<DeviceType>()),
      _mu(atomKK->k_mu.view<DeviceType>()),
      _sp(atomKK->k_sp.view<DeviceType>()),
      _radius(atomKK->k_radius.view<DeviceType>()),
      _rmass(atomKK->k_rmass.view<DeviceType>()),
      _dpdTheta(atomKK->k_dpdTheta.view<DeviceType>()),
      _uCond(atomKK->k_uCond.view<DeviceType>()),
      _uMech(atomKK->k_uMech.view<DeviceType>()),
      _uChem(atomKK->k_uChem.view<DeviceType>()),
      _uCG(atomKK->k_uCG.view<DeviceType>()),
      _uCGnew(atomKK->k_uCGnew.view<DeviceType>()),
      _dx(dx),_dy(dy),_dz(dz),_datamask(datamask) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(i);
    int m = 0;
    if constexpr (PBC_FLAG == 0) {
      _buf(i,m++) = _x(j,0);
      _buf(i,m++) = _x(j,1);
      _buf(i,m++) = _x(j,2);
    } else {
      _buf(i,m++) = _x(j,0) + _dx;
      _buf(i,m++) = _x(j,1) + _dy;
      _buf(i,m++) = _x(j,2) + _dz;
    }

    _buf(i,m++) = d_ubuf(_tag(j)).d;
    _buf(i,m++) = d_ubuf(_type(j)).d;
    _buf(i,m++) = d_ubuf(_mask(j)).d;

    if constexpr (!DEFAULT) {

      if (_datamask & MOLECULE_MASK)
        _buf(i,m++) = d_ubuf(_molecule(j)).d;

      if (_datamask & Q_MASK)
        _buf(i,m++) = _q(j);

      if (_datamask & MU_MASK) {
        _buf(i,m++) = _mu(j,0);
        _buf(i,m++) = _mu(j,1);
        _buf(i,m++) = _mu(j,2);
        _buf(i,m++) = _mu(j,3);
      }

      if (_datamask & SP_MASK) {
        _buf(i,m++) = _sp(j,0);
        _buf(i,m++) = _sp(j,1);
        _buf(i,m++) = _sp(j,2);
        _buf(i,m++) = _sp(j,3);
      }

      if (_datamask & RADIUS_MASK)
        _buf(i,m++) = _radius(j);

      if (_datamask & RMASS_MASK)
        _buf(i,m++) = _rmass(j);

      // DPD-REACT package

      if (_datamask & DPDTHETA_MASK) {
        _buf(i,m++) = _dpdTheta(j);
        _buf(i,m++) = _uCond(j);
        _buf(i,m++) = _uMech(j);
        _buf(i,m++) = _uChem(j);
        _buf(i,m++) = _uCG(j);
        _buf(i,m++) = _uCGnew(j);
      }
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecKokkos::pack_border_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                                               DAT::tdual_double_2d_lr buf,
                                               int pbc_flag, int *pbc, ExecutionSpace space)
{
  atomKK->sync(space,datamask_border);

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
    if (space == HostKK) {
      if (!nborder) {
        AtomVecKokkos_PackBorder<LMPHostType,1,1> f(
          atomKK,buf.view_host(), k_sendlist.view_host(),
          dx,dy,dz,datamask_border);
        Kokkos::parallel_for(n,f);
      } else {
        AtomVecKokkos_PackBorder<LMPHostType,1,0> f(
          atomKK,buf.view_host(), k_sendlist.view_host(),
          dx,dy,dz,datamask_border);
        Kokkos::parallel_for(n,f);
      }
    } else {
      if (!nborder) {
        AtomVecKokkos_PackBorder<LMPDeviceType,1,1> f(
          atomKK,buf.view_device(), k_sendlist.view_device(),
          dx,dy,dz,datamask_border);
        Kokkos::parallel_for(n,f);
      } else {
        AtomVecKokkos_PackBorder<LMPDeviceType,1,0> f(
          atomKK,buf.view_device(), k_sendlist.view_device(),
          dx,dy,dz,datamask_border);
        Kokkos::parallel_for(n,f);
      }
    }
  } else {
    dx = dy = dz = 0;
    if (space == HostKK) {
      if (!nborder) {
        AtomVecKokkos_PackBorder<LMPHostType,0,1> f(
          atomKK,buf.view_host(), k_sendlist.view_host(),
          dx,dy,dz,datamask_border);
        Kokkos::parallel_for(n,f);
      } else {
        AtomVecKokkos_PackBorder<LMPHostType,0,0> f(
          atomKK,buf.view_host(), k_sendlist.view_host(),
          dx,dy,dz,datamask_border);
        Kokkos::parallel_for(n,f);
      }
    } else {
      if (!nborder) {
        AtomVecKokkos_PackBorder<LMPDeviceType,0,1> f(
          atomKK,buf.view_device(), k_sendlist.view_device(),
          dx,dy,dz,datamask_border);
        Kokkos::parallel_for(n,f);
      } else {
        AtomVecKokkos_PackBorder<LMPDeviceType,0,0> f(
          atomKK,buf.view_device(), k_sendlist.view_device(),
          dx,dy,dz,datamask_border);
        Kokkos::parallel_for(n,f);
      }
    }
  }
  return n*size_border;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int DEFAULT>
struct AtomVecKokkos_UnpackBorder {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  const typename AT::t_double_2d_lr_const _buf;
  typename AT::t_kkfloat_1d_3_lr _x;
  typename AT::t_tagint_1d _tag;
  typename AT::t_int_1d _type;
  typename AT::t_int_1d _mask;
  typename AT::t_tagint_1d _molecule;
  typename AT::t_kkfloat_1d _q;
  typename AT::t_kkfloat_1d_4 _mu;
  typename AT::t_kkfloat_1d_4 _sp;
  typename AT::t_kkfloat_1d _radius,_rmass;
  typename AT::t_kkfloat_1d _dpdTheta,_uCond,_uMech,_uChem,_uCG,_uCGnew;
  int _first;
  uint64_t _datamask;

  AtomVecKokkos_UnpackBorder(
    const AtomKokkos* atomKK,
    const typename AT::t_double_2d_lr_const &buf,
    const int &first, const uint64_t &datamask):
    _buf(buf),
    _x(atomKK->k_x.view<DeviceType>()),
    _tag(atomKK->k_tag.view<DeviceType>()),
    _type(atomKK->k_type.view<DeviceType>()),
    _mask(atomKK->k_mask.view<DeviceType>()),
    _molecule(atomKK->k_molecule.view<DeviceType>()),
    _q(atomKK->k_q.view<DeviceType>()),
    _mu(atomKK->k_mu.view<DeviceType>()),
    _sp(atomKK->k_sp.view<DeviceType>()),
    _radius(atomKK->k_radius.view<DeviceType>()),
    _rmass(atomKK->k_rmass.view<DeviceType>()),
    _dpdTheta(atomKK->k_dpdTheta.view<DeviceType>()),
    _uCond(atomKK->k_uCond.view<DeviceType>()),
    _uMech(atomKK->k_uMech.view<DeviceType>()),
    _uChem(atomKK->k_uChem.view<DeviceType>()),
    _uCG(atomKK->k_uCG.view<DeviceType>()),
    _uCGnew(atomKK->k_uCGnew.view<DeviceType>()),
    _first(first),_datamask(datamask) {
  };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    int m = 0;
    _x(i+_first,0) = _buf(i,m++);
    _x(i+_first,1) = _buf(i,m++);
    _x(i+_first,2) = _buf(i,m++);
    _tag(i+_first) = (tagint) d_ubuf(_buf(i,m++)).i;
    _type(i+_first) = (int) d_ubuf(_buf(i,m++)).i;
    _mask(i+_first) = (int) d_ubuf(_buf(i,m++)).i;

    if constexpr (!DEFAULT) {

      if (_datamask & MOLECULE_MASK)
        _molecule(i+_first) = (tagint) d_ubuf(_buf(i,m++)).i;

      if (_datamask & Q_MASK)
        _q(i+_first) = _buf(i,m++);

      if (_datamask & MU_MASK) {
        _mu(i+_first,0) = _buf(i,m++);
        _mu(i+_first,1) = _buf(i,m++);
        _mu(i+_first,2) = _buf(i,m++);
        _mu(i+_first,3) = _buf(i,m++);
      }

      if (_datamask & SP_MASK) {
        _sp(i+_first,0) = _buf(i,m++);
        _sp(i+_first,1) = _buf(i,m++);
        _sp(i+_first,2) = _buf(i,m++);
        _sp(i+_first,3) = _buf(i,m++);
      }

      if (_datamask & RADIUS_MASK)
        _radius(i+_first) = _buf(i,m++);

      if (_datamask & RMASS_MASK)
        _rmass(i+_first) = _buf(i,m++);

      // DPD-REACT package

      if (_datamask & DPDTHETA_MASK) {
        _dpdTheta(i+_first) = _buf(i,m++);
        _uCond(i+_first) = _buf(i,m++);
        _uMech(i+_first) = _buf(i,m++);
        _uChem(i+_first) = _buf(i,m++);
        _uCG(i+_first) = _buf(i,m++);
        _uCGnew(i+_first) = _buf(i,m++);
      }
    }
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecKokkos::unpack_border_kokkos(const int &n, const int &first,
                                               const DAT::tdual_double_2d_lr &buf,
                                               ExecutionSpace space) {
  while (first+n >= nmax) grow(0);

  atomKK->sync(space,datamask_border);

  if (space == HostKK) {
    if (!nborder) {
      struct AtomVecKokkos_UnpackBorder<LMPHostType,1>
        f(atomKK,buf.view_host(),first,datamask_border);
      Kokkos::parallel_for(n,f);
    } else {
      struct AtomVecKokkos_UnpackBorder<LMPHostType,0>
        f(atomKK,buf.view_host(),first,datamask_border);
      Kokkos::parallel_for(n,f);
    }
  } else {
    if (!nborder) {
      struct AtomVecKokkos_UnpackBorder<LMPDeviceType,1>
        f(atomKK,buf.view_device(),first,datamask_border);
      Kokkos::parallel_for(n,f);
    } else {
      struct AtomVecKokkos_UnpackBorder<LMPDeviceType,0>
        f(atomKK,buf.view_device(),first,datamask_border);
      Kokkos::parallel_for(n,f);
    }
  }

  atomKK->modified(space,datamask_border);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG,int DEFORM_VREMAP>
struct AtomVecKokkos_PackBorderVel {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_double_2d_lr_um _buf;
  const typename AT::t_int_1d_const _list;
  const typename AT::t_kkfloat_1d_3_lr_randomread _x;
  typename AT::t_kkfloat_1d_3_randomread _v;
  const typename AT::t_tagint_1d_randomread _tag;
  const typename AT::t_int_1d_randomread _type;
  const typename AT::t_int_1d_randomread _mask;
  const typename AT::t_tagint_1d_randomread _molecule;
  const typename AT::t_kkfloat_1d_randomread _q;
  const typename AT::t_kkfloat_1d_4_randomread _mu;
  const typename AT::t_kkfloat_1d_4_randomread _sp;
  typename AT::t_kkfloat_1d_randomread _radius,_rmass;
  typename AT::t_kkfloat_1d_3_randomread _omega;
  typename AT::t_kkfloat_1d_randomread _dpdTheta,_uCond,_uMech,_uChem,_uCG,_uCGnew;
  double _dx,_dy,_dz, _dvx, _dvy, _dvz;
  const int _deform_groupbit;
  const uint64_t _datamask;

  AtomVecKokkos_PackBorderVel(
    const AtomKokkos* atomKK,
    const typename AT::t_double_2d_lr &buf,
    const typename AT::t_int_1d_const &list,
    const double &dx, const double &dy, const double &dz,
    const double &dvx, const double &dvy, const double &dvz,
    const int &deform_groupbit,
    const uint64_t &datamask):
      _buf(buf),_list(list),_datamask(datamask),
      _x(atomKK->k_x.view<DeviceType>()),
      _tag(atomKK->k_tag.view<DeviceType>()),
      _type(atomKK->k_type.view<DeviceType>()),
      _mask(atomKK->k_mask.view<DeviceType>()),
      _molecule(atomKK->k_molecule.view<DeviceType>()),
      _q(atomKK->k_q.view<DeviceType>()),
      _v(atomKK->k_v.view<DeviceType>()),
      _mu(atomKK->k_mu.view<DeviceType>()),
      _sp(atomKK->k_sp.view<DeviceType>()),
      _radius(atomKK->k_radius.view<DeviceType>()),
      _rmass(atomKK->k_rmass.view<DeviceType>()),
      _omega(atomKK->k_omega.view<DeviceType>()),
      _dpdTheta(atomKK->k_dpdTheta.view<DeviceType>()),
      _uCond(atomKK->k_uCond.view<DeviceType>()),
      _uMech(atomKK->k_uMech.view<DeviceType>()),
      _uChem(atomKK->k_uChem.view<DeviceType>()),
      _uCG(atomKK->k_uCG.view<DeviceType>()),
      _uCGnew(atomKK->k_uCGnew.view<DeviceType>()),
      _dx(dx),_dy(dy),_dz(dz),
      _dvx(dvx),_dvy(dvy),_dvz(dvz),
      _deform_groupbit(deform_groupbit) {
        const size_t elements = atomKK->avecKK->size_border + atomKK->avecKK->size_velocity;
        const int maxsend = (buf.extent(0)*buf.extent(1))/elements;
        _buf = typename AT::t_double_2d_lr_um(buf.data(),maxsend,elements);
      }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    int m = 0;
    const int j = _list(i);
    if constexpr (PBC_FLAG == 0) {
      _buf(i,m++) = _x(j,0);
      _buf(i,m++) = _x(j,1);
      _buf(i,m++) = _x(j,2);
    } else {
      _buf(i,m++) = _x(j,0) + _dx;
      _buf(i,m++) = _x(j,1) + _dy;
      _buf(i,m++) = _x(j,2) + _dz;
    }
    _buf(i,m++) = d_ubuf(_tag(j)).d;
    _buf(i,m++) = d_ubuf(_type(j)).d;
    _buf(i,m++) = d_ubuf(_mask(j)).d;

    if constexpr (DEFORM_VREMAP) {
      if (_mask(i) & _deform_groupbit) {
        _buf(i,m++) = _v(j,0) + _dvx;
        _buf(i,m++) = _v(j,1) + _dvy;
        _buf(i,m++) = _v(j,2) + _dvz;
      }
    } else {
      _buf(i,m++) = _v(j,0);
      _buf(i,m++) = _v(j,1);
      _buf(i,m++) = _v(j,2);
    }

    if (_datamask & MOLECULE_MASK)
      _buf(i,m++) = d_ubuf(_molecule(j)).d;

    if (_datamask & Q_MASK)
      _buf(i,m++) = _q(j);

    if (_datamask & MU_MASK) {
      _buf(i,m++) = _mu(j,0);
      _buf(i,m++) = _mu(j,1);
      _buf(i,m++) = _mu(j,2);
      _buf(i,m++) = _mu(j,3);
    }

    if (_datamask & SP_MASK) {
      _buf(i,m++) = _sp(j,0);
      _buf(i,m++) = _sp(j,1);
      _buf(i,m++) = _sp(j,2);
      _buf(i,m++) = _sp(j,3);
    }

    if (_datamask & RADIUS_MASK)
      _buf(i,m++) = _radius(j);

    if (_datamask & RMASS_MASK)
      _buf(i,m++) = _rmass(j);

    if (_datamask & OMEGA_MASK) {
      _buf(i,m++) = _omega(j,0);
      _buf(i,m++) = _omega(j,1);
      _buf(i,m++) = _omega(j,2);
    }

    // DPD-REACT package

    if (_datamask & DPDTHETA_MASK) {
      _buf(i,m++) = _dpdTheta(j);
      _buf(i,m++) = _uCond(j);
      _buf(i,m++) = _uMech(j);
      _buf(i,m++) = _uChem(j);
      _buf(i,m++) = _uCG(j);
      _buf(i,m++) = _uCGnew(j);
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecKokkos::pack_border_vel_kokkos(
  int n, DAT::tdual_int_1d k_sendlist, DAT::tdual_double_2d_lr buf,
  int pbc_flag, int *pbc, ExecutionSpace space)
{
  double dx = 0, dy = 0, dz = 0;
  double dvx = 0, dvy = 0, dvz = 0;

  atomKK->sync(space,datamask_border_vel);

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
      if (space == HostKK) {
        AtomVecKokkos_PackBorderVel<LMPHostType,1,0> f(
          atomKK,
          buf.view_host(), k_sendlist.view_host(),
          dx,dy,dz,dvx,dvy,dvz,
          deform_groupbit,datamask_border_vel);
        Kokkos::parallel_for(n,f);
      } else {
        AtomVecKokkos_PackBorderVel<LMPDeviceType,1,0> f(
          atomKK,
          buf.view_device(), k_sendlist.view_device(),
          dx,dy,dz,dvx,dvy,dvz,
          deform_groupbit,datamask_border_vel);
        Kokkos::parallel_for(n,f);
      }
    }
    else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      if (space == HostKK) {
        AtomVecKokkos_PackBorderVel<LMPHostType,1,1> f(
          atomKK,
          buf.view_host(), k_sendlist.view_host(),
          dx,dy,dz,dvx,dvy,dvz,
          deform_groupbit,datamask_border_vel);
        Kokkos::parallel_for(n,f);
      } else {
        AtomVecKokkos_PackBorderVel<LMPDeviceType,1,1> f(
          atomKK,
          buf.view_device(), k_sendlist.view_device(),
          dx,dy,dz,dvx,dvy,dvz,
          deform_groupbit,datamask_border_vel);
        Kokkos::parallel_for(n,f);
      }
    }
  } else {
    if (space == HostKK) {
      AtomVecKokkos_PackBorderVel<LMPHostType,0,0> f(
        atomKK,
        buf.view_host(), k_sendlist.view_host(),
        dx,dy,dz,dvx,dvy,dvz,
        deform_groupbit,datamask_border_vel);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecKokkos_PackBorderVel<LMPDeviceType,0,0> f(
        atomKK,
        buf.view_device(), k_sendlist.view_device(),
        dx,dy,dz,dvx,dvy,dvz,
        deform_groupbit,datamask_border_vel);
      Kokkos::parallel_for(n,f);
    }
  }

  return n*(size_border + size_velocity);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int DEFAULT>
struct AtomVecKokkos_UnpackBorderVel {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_double_2d_lr_const_um _buf;
  typename AT::t_kkfloat_1d_3_lr _x;
  typename AT::t_tagint_1d _tag;
  typename AT::t_int_1d _type;
  typename AT::t_int_1d _mask;
  typename AT::t_tagint_1d _molecule;
  typename AT::t_kkfloat_1d _q;
  typename AT::t_kkfloat_1d_3 _v;
  typename AT::t_kkfloat_1d_4 _mu;
  typename AT::t_kkfloat_1d_4 _sp;
  typename AT::t_kkfloat_1d _radius,_rmass;
  typename AT::t_kkfloat_1d_3 _omega;
  typename AT::t_kkfloat_1d _dpdTheta,_uCond,_uMech,_uChem,_uCG,_uCGnew;
  int _first;
  uint64_t _datamask;

  AtomVecKokkos_UnpackBorderVel(
    const AtomKokkos* atomKK,
    const typename AT::t_double_2d_lr_const &buf,
    const int &first,
    const uint64_t &datamask):
    _buf(buf),
    _x(atomKK->k_x.view<DeviceType>()),
    _tag(atomKK->k_tag.view<DeviceType>()),
    _type(atomKK->k_type.view<DeviceType>()),
    _mask(atomKK->k_mask.view<DeviceType>()),
    _molecule(atomKK->k_molecule.view<DeviceType>()),
    _q(atomKK->k_q.view<DeviceType>()),
    _v(atomKK->k_v.view<DeviceType>()),
    _mu(atomKK->k_mu.view<DeviceType>()),
    _sp(atomKK->k_sp.view<DeviceType>()),
    _radius(atomKK->k_radius.view<DeviceType>()),
    _rmass(atomKK->k_rmass.view<DeviceType>()),
    _omega(atomKK->k_omega.view<DeviceType>()),
    _dpdTheta(atomKK->k_dpdTheta.view<DeviceType>()),
    _uCond(atomKK->k_uCond.view<DeviceType>()),
    _uMech(atomKK->k_uMech.view<DeviceType>()),
    _uChem(atomKK->k_uChem.view<DeviceType>()),
    _uCG(atomKK->k_uCG.view<DeviceType>()),
    _uCGnew(atomKK->k_uCGnew.view<DeviceType>()),
    _first(first),_datamask(datamask)
  {
    const size_t elements = atomKK->avecKK->size_border + atomKK->avecKK->size_velocity;
    const int maxsend = (buf.extent(0)*buf.extent(1))/elements;
    _buf = typename AT::t_double_2d_lr_const_um(buf.data(),maxsend,elements);
  };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    int m = 0;
    _x(i+_first,0) = _buf(i,m++);
    _x(i+_first,1) = _buf(i,m++);
    _x(i+_first,2) = _buf(i,m++);
    _tag(i+_first) = static_cast<tagint>(d_ubuf(_buf(i,m++)).i);
    _type(i+_first) = static_cast<int>(d_ubuf(_buf(i,m++)).i);
    _mask(i+_first) = static_cast<int>(d_ubuf(_buf(i,m++)).i);
    _v(i+_first,0) = _buf(i,m++);
    _v(i+_first,1) = _buf(i,m++);
    _v(i+_first,2) = _buf(i,m++);

    if constexpr (!DEFAULT) {

      if (_datamask & MOLECULE_MASK)
        _molecule(i+_first) = (tagint) d_ubuf(_buf(i,m++)).i;

      if (_datamask & Q_MASK)
        _q(i+_first) = _buf(i,m++);

      if (_datamask & MU_MASK) {
        _mu(i+_first,0) = _buf(i,m++);
        _mu(i+_first,1) = _buf(i,m++);
        _mu(i+_first,2) = _buf(i,m++);
        _mu(i+_first,3) = _buf(i,m++);
      }

      if (_datamask & SP_MASK) {
        _sp(i+_first,0) = _buf(i,m++);
        _sp(i+_first,1) = _buf(i,m++);
        _sp(i+_first,2) = _buf(i,m++);
        _sp(i+_first,3) = _buf(i,m++);
      }

      if (_datamask & RADIUS_MASK)
        _radius(i+_first) = _buf(i,m++);

      if (_datamask & RMASS_MASK)
        _rmass(i+_first) = _buf(i,m++);

      if (_datamask & OMEGA_MASK) {
        _omega(i+_first,0) = _buf(i,m++);
        _omega(i+_first,1) = _buf(i,m++);
        _omega(i+_first,2) = _buf(i,m++);
      }

      // DPD-REACT package

      if (_datamask & DPDTHETA_MASK) {
        _dpdTheta(i+_first) = _buf(i,m++);
        _uCond(i+_first) = _buf(i,m++);
        _uMech(i+_first) = _buf(i,m++);
        _uChem(i+_first) = _buf(i,m++);
        _uCG(i+_first) = _buf(i,m++);
        _uCGnew(i+_first) = _buf(i,m++);
      }
    }
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecKokkos::unpack_border_vel_kokkos(
    const int &n, const int &first,
    const DAT::tdual_double_2d_lr &buf,ExecutionSpace space) {

  while (first+n >= nmax) grow(0);

  atomKK->sync(space,datamask_border_vel);

  if (space == HostKK) {
    if (!ncomm_vel) {
      struct AtomVecKokkos_UnpackBorderVel<LMPHostType,1> f(
        atomKK,
        buf.view_host(),
        first,datamask_border_vel);
      Kokkos::parallel_for(n,f);
    } else {
      struct AtomVecKokkos_UnpackBorderVel<LMPHostType,0> f(
        atomKK,
        buf.view_host(),
        first,datamask_border_vel);
      Kokkos::parallel_for(n,f);
    }
  } else {
    if (!ncomm_vel) {
      struct AtomVecKokkos_UnpackBorderVel<LMPDeviceType,1> f(
        atomKK,
        buf.view_device(),
        first,datamask_border_vel);
      Kokkos::parallel_for(n,f);
    } else {
      struct AtomVecKokkos_UnpackBorderVel<LMPDeviceType,0> f(
        atomKK,
        buf.view_device(),
        first,datamask_border_vel);
      Kokkos::parallel_for(n,f);
    }
  }

  atomKK->modified(space,datamask_border_vel);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int DEFAULT>
struct AtomVecKokkos_PackExchangeFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr _x;
  typename AT::t_kkfloat_1d_3 _v;
  typename AT::t_tagint_1d _tag;
  typename AT::t_int_1d _type;
  typename AT::t_int_1d _mask;
  typename AT::t_imageint_1d _image;
  typename AT::t_kkfloat_1d _q;
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
  typename AT::t_kkfloat_1d_4 _mu;
  typename AT::t_kkfloat_1d_4 _sp;
  typename AT::t_kkfloat_1d _radius,_rmass;
  typename AT::t_kkfloat_1d_3 _omega;
  typename AT::t_kkfloat_1d _dpdTheta,_uCond,_uMech,_uChem,_uCG,_uCGnew;

  typename AT::t_double_2d_lr_um _buf;
  typename AT::t_int_1d_const _sendlist;
  typename AT::t_int_1d_const _copylist;
  int _size_exchange;
  uint64_t _datamask;

  AtomVecKokkos_PackExchangeFunctor(
    const AtomKokkos* atomKK,
    const DAT::tdual_double_2d_lr buf,
    DAT::tdual_int_1d sendlist,
    DAT::tdual_int_1d copylist,
    const uint64_t datamask):
      _x(atomKK->k_x.view<DeviceType>()),
      _v(atomKK->k_v.view<DeviceType>()),
      _tag(atomKK->k_tag.view<DeviceType>()),
      _type(atomKK->k_type.view<DeviceType>()),
      _mask(atomKK->k_mask.view<DeviceType>()),
      _image(atomKK->k_image.view<DeviceType>()),
      _q(atomKK->k_q.view<DeviceType>()),
      _molecule(atomKK->k_molecule.view<DeviceType>()),
      _nspecial(atomKK->k_nspecial.view<DeviceType>()),
      _special(atomKK->k_special.view<DeviceType>()),
      _num_bond(atomKK->k_num_bond.view<DeviceType>()),
      _bond_type(atomKK->k_bond_type.view<DeviceType>()),
      _bond_atom(atomKK->k_bond_atom.view<DeviceType>()),
      _num_angle(atomKK->k_num_angle.view<DeviceType>()),
      _angle_type(atomKK->k_angle_type.view<DeviceType>()),
      _angle_atom1(atomKK->k_angle_atom1.view<DeviceType>()),
      _angle_atom2(atomKK->k_angle_atom2.view<DeviceType>()),
      _angle_atom3(atomKK->k_angle_atom3.view<DeviceType>()),
      _num_dihedral(atomKK->k_num_dihedral.view<DeviceType>()),
      _dihedral_type(atomKK->k_dihedral_type.view<DeviceType>()),
      _dihedral_atom1(atomKK->k_dihedral_atom1.view<DeviceType>()),
      _dihedral_atom2(atomKK->k_dihedral_atom2.view<DeviceType>()),
      _dihedral_atom3(atomKK->k_dihedral_atom3.view<DeviceType>()),
      _dihedral_atom4(atomKK->k_dihedral_atom4.view<DeviceType>()),
      _num_improper(atomKK->k_num_improper.view<DeviceType>()),
      _improper_type(atomKK->k_improper_type.view<DeviceType>()),
      _improper_atom1(atomKK->k_improper_atom1.view<DeviceType>()),
      _improper_atom2(atomKK->k_improper_atom2.view<DeviceType>()),
      _improper_atom3(atomKK->k_improper_atom3.view<DeviceType>()),
      _improper_atom4(atomKK->k_improper_atom4.view<DeviceType>()),
      _mu(atomKK->k_mu.view<DeviceType>()),
      _sp(atomKK->k_sp.view<DeviceType>()),
      _radius(atomKK->k_radius.view<DeviceType>()),
      _rmass(atomKK->k_rmass.view<DeviceType>()),
      _omega(atomKK->k_omega.view<DeviceType>()),
      _dpdTheta(atomKK->k_dpdTheta.view<DeviceType>()),
      _uCond(atomKK->k_uCond.view<DeviceType>()),
      _uMech(atomKK->k_uMech.view<DeviceType>()),
      _uChem(atomKK->k_uChem.view<DeviceType>()),
      _uCG(atomKK->k_uCG.view<DeviceType>()),
      _uCGnew(atomKK->k_uCGnew.view<DeviceType>()),

      _sendlist(sendlist.template view<DeviceType>()),
      _copylist(copylist.template view<DeviceType>()),
      _size_exchange(atomKK->avecKK->size_exchange),
      _datamask(datamask) {
        const int maxsendlist = (buf.template view<DeviceType>().extent(0)*
                                 buf.template view<DeviceType>().extent(1))/_size_exchange;
        buffer_view<DeviceType>(_buf,buf,maxsendlist,_size_exchange);
      }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &mysend) const {
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

    if constexpr (!DEFAULT) {

      if (_datamask & Q_MASK)
        _buf(mysend,m++) = _q(i);

      if (_datamask & MOLECULE_MASK)
        _buf(mysend,m++) = d_ubuf(_molecule(i)).d;

      if (_datamask & BOND_MASK) {
        _buf(mysend,m++) = d_ubuf(_num_bond(i)).d;
        for (int k = 0; k < _num_bond(i); k++) {
          _buf(mysend,m++) = d_ubuf(_bond_type(i,k)).d;
          _buf(mysend,m++) = d_ubuf(_bond_atom(i,k)).d;
        }
      }

      if (_datamask & ANGLE_MASK) {
        _buf(mysend,m++) = d_ubuf(_num_angle(i)).d;
        for (int k = 0; k < _num_angle(i); k++) {
          _buf(mysend,m++) = d_ubuf(_angle_type(i,k)).d;
          _buf(mysend,m++) = d_ubuf(_angle_atom1(i,k)).d;
          _buf(mysend,m++) = d_ubuf(_angle_atom2(i,k)).d;
          _buf(mysend,m++) = d_ubuf(_angle_atom3(i,k)).d;
        }
      }

      if (_datamask & DIHEDRAL_MASK) {
        _buf(mysend,m++) = d_ubuf(_num_dihedral(i)).d;
        for (int k = 0; k < _num_dihedral(i); k++) {
          _buf(mysend,m++) = d_ubuf(_dihedral_type(i,k)).d;
          _buf(mysend,m++) = d_ubuf(_dihedral_atom1(i,k)).d;
          _buf(mysend,m++) = d_ubuf(_dihedral_atom2(i,k)).d;
          _buf(mysend,m++) = d_ubuf(_dihedral_atom3(i,k)).d;
          _buf(mysend,m++) = d_ubuf(_dihedral_atom4(i,k)).d;
        }
      }

      if (_datamask & IMPROPER_MASK) {
        _buf(mysend,m++) = d_ubuf(_num_improper(i)).d;
        for (int k = 0; k < _num_improper(i); k++) {
          _buf(mysend,m++) = d_ubuf(_improper_type(i,k)).d;
          _buf(mysend,m++) = d_ubuf(_improper_atom1(i,k)).d;
          _buf(mysend,m++) = d_ubuf(_improper_atom2(i,k)).d;
          _buf(mysend,m++) = d_ubuf(_improper_atom3(i,k)).d;
          _buf(mysend,m++) = d_ubuf(_improper_atom4(i,k)).d;
        }
      }

      if (_datamask & SPECIAL_MASK) {
        _buf(mysend,m++) = d_ubuf(_nspecial(i,0)).d;
        _buf(mysend,m++) = d_ubuf(_nspecial(i,1)).d;
        _buf(mysend,m++) = d_ubuf(_nspecial(i,2)).d;
        for (int k = 0; k < _nspecial(i,2); k++)
          _buf(mysend,m++) = d_ubuf(_special(i,k)).d;
      }

      if (_datamask & MU_MASK) {
        _buf(mysend,m++) = _mu(i,0);
        _buf(mysend,m++) = _mu(i,1);
        _buf(mysend,m++) = _mu(i,2);
        _buf(mysend,m++) = _mu(i,3);
      }

      if (_datamask & SP_MASK) {
        _buf(mysend,m++) = _sp(i,0);
        _buf(mysend,m++) = _sp(i,1);
        _buf(mysend,m++) = _sp(i,2);
        _buf(mysend,m++) = _sp(i,3);
      }

      if (_datamask & RADIUS_MASK)
        _buf(mysend,m++) = _radius(i);

      if (_datamask & RMASS_MASK)
        _buf(mysend,m++) = _rmass(i);

      if (_datamask & OMEGA_MASK) {
        _buf(mysend,m++) = _omega(i,0);
        _buf(mysend,m++) = _omega(i,1);
        _buf(mysend,m++) = _omega(i,2);
      }

      // DPD-REACT package

      if (_datamask & DPDTHETA_MASK) {
        _buf(mysend,m++) = _dpdTheta(i);
        _buf(mysend,m++) = _uCond(i);
        _buf(mysend,m++) = _uMech(i);
        _buf(mysend,m++) = _uChem(i);
        _buf(mysend,m++) = _uCG(i);
        _buf(mysend,m++) = _uCGnew(i);
      }
    }

    const int j = _copylist(mysend);

    if (j > -1) {
      _x(i,0) = _x(j,0);
      _x(i,1) = _x(j,1);
      _x(i,2) = _x(j,2);
      _v(i,0) = _v(j,0);
      _v(i,1) = _v(j,1);
      _v(i,2) = _v(j,2);
      _tag(i) = _tag(j);
      _type(i) = _type(j);
      _mask(i) = _mask(j);
      _image(i) = _image(j);

      if constexpr (!DEFAULT) {

        if (_datamask & Q_MASK)
          _q(i) = _q(j);

        if (_datamask & MOLECULE_MASK)
          _molecule(i) = _molecule(j);

        if (_datamask & BOND_MASK) {
          _num_bond(i) = _num_bond(j);
          for (int k = 0; k < _num_bond(j); k++) {
            _bond_type(i,k) = _bond_type(j,k);
            _bond_atom(i,k) = _bond_atom(j,k);
          }
        }

        if (_datamask & ANGLE_MASK) {
          _num_angle(i) = _num_angle(j);
          for (int k = 0; k < _num_angle(j); k++) {
            _angle_type(i,k) = _angle_type(j,k);
            _angle_atom1(i,k) = _angle_atom1(j,k);
            _angle_atom2(i,k) = _angle_atom2(j,k);
            _angle_atom3(i,k) = _angle_atom3(j,k);
          }
        }

        if (_datamask & DIHEDRAL_MASK) {
          _num_dihedral(i) = _num_dihedral(j);
          for (int k = 0; k < _num_dihedral(j); k++) {
            _dihedral_type(i,k) = _dihedral_type(j,k);
            _dihedral_atom1(i,k) = _dihedral_atom1(j,k);
            _dihedral_atom2(i,k) = _dihedral_atom2(j,k);
            _dihedral_atom3(i,k) = _dihedral_atom3(j,k);
            _dihedral_atom4(i,k) = _dihedral_atom4(j,k);
          }
        }

        if (_datamask & IMPROPER_MASK) {
          _num_improper(i) = _num_improper(j);
          for (int k = 0; k < _num_improper(j); k++) {
            _improper_type(i,k) = _improper_type(j,k);
            _improper_atom1(i,k) = _improper_atom1(j,k);
            _improper_atom2(i,k) = _improper_atom2(j,k);
            _improper_atom3(i,k) = _improper_atom3(j,k);
            _improper_atom4(i,k) = _improper_atom4(j,k);
          }
        }

        if (_datamask & SPECIAL_MASK) {
          _nspecial(i,0) = _nspecial(j,0);
          _nspecial(i,1) = _nspecial(j,1);
          _nspecial(i,2) = _nspecial(j,2);
          for (int k = 0; k < _nspecial(j,2); k++)
            _special(i,k) = _special(j,k);
        }

        if (_datamask & MU_MASK) {
          _mu(i,0) = _mu(j,0);
          _mu(i,1) = _mu(j,1);
          _mu(i,2) = _mu(j,2);
          _mu(i,3) = _mu(j,3);
        }

        if (_datamask & SP_MASK) {
          _sp(i,0) = _sp(j,0);
          _sp(i,1) = _sp(j,1);
          _sp(i,2) = _sp(j,2);
          _sp(i,3) = _sp(j,3);
        }

        if (_datamask & RADIUS_MASK)
          _radius(i) = _radius(j);

        if (_datamask & RMASS_MASK)
          _rmass(i) = _rmass(j);

        if (_datamask & OMEGA_MASK) {
          _omega(i,0) = _omega(j,0);
          _omega(i,1) = _omega(j,1);
          _omega(i,2) = _omega(j,2);
        }

        // DPD-REACT package

        if (_datamask & DPDTHETA_MASK) {
           _dpdTheta(i) = _dpdTheta(j);
          _uCond(i) = _uCond(j);
          _uMech(i) = _uMech(j);
          _uChem(i) = _uChem(j);
          _uCG(i) = _uCG(j);
          _uCGnew(i) = _uCGnew(j);
        }
      }
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecKokkos::pack_exchange_kokkos(const int &nsend,DAT::tdual_double_2d_lr &k_buf,
                                                 DAT::tdual_int_1d k_sendlist,
                                                 DAT::tdual_int_1d k_copylist,
                                                 ExecutionSpace space)
{
  atomKK->sync(space,datamask_exchange);
  set_size_exchange();

  if (nsend > (int) (k_buf.view_host().extent(0)*
              k_buf.view_host().extent(1))/size_exchange) {
    int newsize = nsend*size_exchange/k_buf.view_host().extent(1)+1;
    k_buf.resize(newsize,k_buf.view_host().extent(1));
  }

  if (space == HostKK) {
    if (size_exchange == size_exchange_default) {
      AtomVecKokkos_PackExchangeFunctor<LMPHostType,1>
        f(atomKK,k_buf,k_sendlist,k_copylist,datamask_exchange);
      Kokkos::parallel_for(nsend,f);
    } else {
      AtomVecKokkos_PackExchangeFunctor<LMPHostType,0>
        f(atomKK,k_buf,k_sendlist,k_copylist,datamask_exchange);
      Kokkos::parallel_for(nsend,f);
    }
  } else {
    if (size_exchange == size_exchange_default) {
      AtomVecKokkos_PackExchangeFunctor<LMPDeviceType,1>
        f(atomKK,k_buf,k_sendlist,k_copylist,datamask_exchange);
      Kokkos::parallel_for(nsend,f);
    } else {
      AtomVecKokkos_PackExchangeFunctor<LMPDeviceType,0>
        f(atomKK,k_buf,k_sendlist,k_copylist,datamask_exchange);
      Kokkos::parallel_for(nsend,f);
    }
  }

  atomKK->modified(space,datamask_exchange);

  return nsend*size_exchange;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int OUTPUT_INDICES,int DEFAULT>
struct AtomVecKokkos_UnpackExchangeFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr _x;
  typename AT::t_kkfloat_1d_3 _v;
  typename AT::t_tagint_1d _tag;
  typename AT::t_int_1d _type;
  typename AT::t_int_1d _mask;
  typename AT::t_imageint_1d _image;
  typename AT::t_kkfloat_1d _q;
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
  typename AT::t_kkfloat_1d_4 _mu;
  typename AT::t_kkfloat_1d_4 _sp;
  typename AT::t_kkfloat_1d _radius,_rmass;
  typename AT::t_kkfloat_1d_3 _omega;
  typename AT::t_kkfloat_1d _dpdTheta,_uCond,_uMech,_uChem,_uCG,_uCGnew;

  typename AT::t_double_2d_lr_um _buf;
  typename AT::t_int_1d _nlocal;
  typename AT::t_int_1d _indices;
  int _dim;
  double _lo,_hi;
  int _size_exchange;
  uint64_t _datamask;

  AtomVecKokkos_UnpackExchangeFunctor(
    const AtomKokkos* atomKK,
    const DAT::tdual_double_2d_lr buf,
    DAT::tdual_int_1d nlocal,
    DAT::tdual_int_1d indices,
    int dim, double lo, double hi,
    uint64_t datamask):
      _x(atomKK->k_x.view<DeviceType>()),
      _v(atomKK->k_v.view<DeviceType>()),
      _tag(atomKK->k_tag.view<DeviceType>()),
      _type(atomKK->k_type.view<DeviceType>()),
      _mask(atomKK->k_mask.view<DeviceType>()),
      _image(atomKK->k_image.view<DeviceType>()),
      _q(atomKK->k_q.view<DeviceType>()),
      _molecule(atomKK->k_molecule.view<DeviceType>()),
      _nspecial(atomKK->k_nspecial.view<DeviceType>()),
      _special(atomKK->k_special.view<DeviceType>()),
      _num_bond(atomKK->k_num_bond.view<DeviceType>()),
      _bond_type(atomKK->k_bond_type.view<DeviceType>()),
      _bond_atom(atomKK->k_bond_atom.view<DeviceType>()),
      _num_angle(atomKK->k_num_angle.view<DeviceType>()),
      _angle_type(atomKK->k_angle_type.view<DeviceType>()),
      _angle_atom1(atomKK->k_angle_atom1.view<DeviceType>()),
      _angle_atom2(atomKK->k_angle_atom2.view<DeviceType>()),
      _angle_atom3(atomKK->k_angle_atom3.view<DeviceType>()),
      _num_dihedral(atomKK->k_num_dihedral.view<DeviceType>()),
      _dihedral_type(atomKK->k_dihedral_type.view<DeviceType>()),
      _dihedral_atom1(atomKK->k_dihedral_atom1.view<DeviceType>()),
      _dihedral_atom2(atomKK->k_dihedral_atom2.view<DeviceType>()),
      _dihedral_atom3(atomKK->k_dihedral_atom3.view<DeviceType>()),
      _dihedral_atom4(atomKK->k_dihedral_atom4.view<DeviceType>()),
      _num_improper(atomKK->k_num_improper.view<DeviceType>()),
      _improper_type(atomKK->k_improper_type.view<DeviceType>()),
      _improper_atom1(atomKK->k_improper_atom1.view<DeviceType>()),
      _improper_atom2(atomKK->k_improper_atom2.view<DeviceType>()),
      _improper_atom3(atomKK->k_improper_atom3.view<DeviceType>()),
      _improper_atom4(atomKK->k_improper_atom4.view<DeviceType>()),
      _mu(atomKK->k_mu.view<DeviceType>()),
      _sp(atomKK->k_sp.view<DeviceType>()),
      _radius(atomKK->k_radius.view<DeviceType>()),
      _rmass(atomKK->k_rmass.view<DeviceType>()),
      _omega(atomKK->k_omega.view<DeviceType>()),
      _dpdTheta(atomKK->k_dpdTheta.view<DeviceType>()),
      _uCond(atomKK->k_uCond.view<DeviceType>()),
      _uMech(atomKK->k_uMech.view<DeviceType>()),
      _uChem(atomKK->k_uChem.view<DeviceType>()),
      _uCG(atomKK->k_uCG.view<DeviceType>()),
      _uCGnew(atomKK->k_uCGnew.view<DeviceType>()),

      _nlocal(nlocal.template view<DeviceType>()),
      _indices(indices.template view<DeviceType>()),
      _dim(dim),_lo(lo),_hi(hi),_size_exchange(atomKK->avecKK->size_exchange),
      _datamask(datamask) {
    const int maxsendlist = (buf.template view<DeviceType>().extent(0)*
                             buf.template view<DeviceType>().extent(1))/_size_exchange;
    buffer_view<DeviceType>(_buf,buf,maxsendlist,_size_exchange);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &myrecv) const {
    double x = _buf(myrecv,_dim+1);
    int i = -1;
    if (x >= _lo && x < _hi) {
      i = Kokkos::atomic_fetch_add(&_nlocal(0),1);
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

      if constexpr (!DEFAULT) {

        if (_datamask & Q_MASK)
          _q(i) = _buf(myrecv,m++);

        if (_datamask & MOLECULE_MASK)
          _molecule(i) = (tagint) d_ubuf(_buf(myrecv,m++)).i;

        if (_datamask & BOND_MASK) {
          _num_bond(i) = (int) d_ubuf(_buf(myrecv,m++)).i;
          for (int k = 0; k < _num_bond(i); k++) {
            _bond_type(i,k) = (int) d_ubuf(_buf(myrecv,m++)).i;
            _bond_atom(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
          }
        }

        if (_datamask & ANGLE_MASK) {
          _num_angle(i) = (int) d_ubuf(_buf(myrecv,m++)).i;
          for (int k = 0; k < _num_angle(i); k++) {
            _angle_type(i,k) = (int) d_ubuf(_buf(myrecv,m++)).i;
            _angle_atom1(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
            _angle_atom2(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
            _angle_atom3(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
          }
        }

        if (_datamask & DIHEDRAL_MASK) {
          _num_dihedral(i) = (int) d_ubuf(_buf(myrecv,m++)).i;
          for (int k = 0; k < _num_dihedral(i); k++) {
            _dihedral_type(i,k) = (int) d_ubuf(_buf(myrecv,m++)).i;
            _dihedral_atom1(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
            _dihedral_atom2(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
            _dihedral_atom3(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
            _dihedral_atom4(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
          }
        }

        if (_datamask & IMPROPER_MASK) {
          _num_improper(i) = (int) d_ubuf(_buf(myrecv,m++)).i;
          for (int k = 0; k < _num_improper(i); k++) {
            _improper_type(i,k) = (int) d_ubuf(_buf(myrecv,m++)).i;
            _improper_atom1(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
            _improper_atom2(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
            _improper_atom3(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
            _improper_atom4(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
          }
        }

        if (_datamask & SPECIAL_MASK) {
          _nspecial(i,0) = (int) d_ubuf(_buf(myrecv,m++)).i;
          _nspecial(i,1) = (int) d_ubuf(_buf(myrecv,m++)).i;
          _nspecial(i,2) = (int) d_ubuf(_buf(myrecv,m++)).i;
          for (int k = 0; k < _nspecial(i,2); k++)
            _special(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
        }

        if (_datamask & MU_MASK) {
          _mu(i,0) = _buf(myrecv,m++);
          _mu(i,1) = _buf(myrecv,m++);
          _mu(i,2) = _buf(myrecv,m++);
          _mu(i,3) = _buf(myrecv,m++);
        }

        if (_datamask & SP_MASK) {
          _sp(i,0) = _buf(myrecv,m++);
          _sp(i,1) = _buf(myrecv,m++);
          _sp(i,2) = _buf(myrecv,m++);
          _sp(i,3) = _buf(myrecv,m++);
        }

        if (_datamask & RADIUS_MASK)
          _radius(i) = _buf(myrecv,m++);

        if (_datamask & RMASS_MASK)
          _rmass(i) = _buf(myrecv,m++);

        if (_datamask & OMEGA_MASK) {
          _omega(i,0) = _buf(myrecv,m++);
          _omega(i,1) = _buf(myrecv,m++);
          _omega(i,2) = _buf(myrecv,m++);
        }

        // DPD-REACT package

        if (_datamask & DPDTHETA_MASK) {
          _dpdTheta(i) = _buf(myrecv,m++);
          _uCond(i) = _buf(myrecv,m++);
          _uMech(i) = _buf(myrecv,m++);
          _uChem(i) = _buf(myrecv,m++);
          _uCG(i) = _buf(myrecv,m++);
          _uCGnew(i) = _buf(myrecv,m++);
        }
      }
    }

    if constexpr (OUTPUT_INDICES)
      _indices(myrecv) = i;
  }
};

/* ---------------------------------------------------------------------- */
int AtomVecKokkos::unpack_exchange_kokkos(DAT::tdual_double_2d_lr &k_buf, int nrecv, int nlocal,
                                              int dim, double lo, double hi, ExecutionSpace space,
                                              DAT::tdual_int_1d &k_indices)
{
  while (nlocal + nrecv/size_exchange >= nmax) grow(0);

  atomKK->sync(space,datamask_exchange);

  if (space == HostKK) {
    k_count.view_host()(0) = nlocal;

    if (k_indices.view_host().data()) {
      if (size_exchange == size_exchange_default) {
        AtomVecKokkos_UnpackExchangeFunctor<LMPHostType,1,1>
          f(atomKK,k_buf,k_count,k_indices,dim,lo,hi,datamask_exchange);
        Kokkos::parallel_for(nrecv/size_exchange,f);
      } else {
        AtomVecKokkos_UnpackExchangeFunctor<LMPHostType,1,0>
          f(atomKK,k_buf,k_count,k_indices,dim,lo,hi,datamask_exchange);
        Kokkos::parallel_for(nrecv/size_exchange,f);
      }
    } else {
      if (size_exchange == size_exchange_default) {
        AtomVecKokkos_UnpackExchangeFunctor<LMPHostType,0,1>
          f(atomKK,k_buf,k_count,k_indices,dim,lo,hi,datamask_exchange);
        Kokkos::parallel_for(nrecv/size_exchange,f);
      } else {
        AtomVecKokkos_UnpackExchangeFunctor<LMPHostType,0,0>
          f(atomKK,k_buf,k_count,k_indices,dim,lo,hi,datamask_exchange);
        Kokkos::parallel_for(nrecv/size_exchange,f);
      }
    }
  } else {
    k_count.view_host()(0) = nlocal;
    k_count.modify_host();
    k_count.sync_device();

    if (k_indices.view_host().data()) {
      if (size_exchange == size_exchange_default) {
        AtomVecKokkos_UnpackExchangeFunctor<LMPDeviceType,1,1>
          f(atomKK,k_buf,k_count,k_indices,dim,lo,hi,datamask_exchange);
        Kokkos::parallel_for(nrecv/size_exchange,f);
      } else {
        AtomVecKokkos_UnpackExchangeFunctor<LMPDeviceType,1,0>
          f(atomKK,k_buf,k_count,k_indices,dim,lo,hi,datamask_exchange);
        Kokkos::parallel_for(nrecv/size_exchange,f);
      }
    } else {
      if (size_exchange == size_exchange_default) {
        AtomVecKokkos_UnpackExchangeFunctor<LMPDeviceType,0,1>
          f(atomKK,k_buf,k_count,k_indices,dim,lo,hi,datamask_exchange);
        Kokkos::parallel_for(nrecv/size_exchange,f);
      } else {
        AtomVecKokkos_UnpackExchangeFunctor<LMPDeviceType,0,0>
          f(atomKK,k_buf,k_count,k_indices,dim,lo,hi,datamask_exchange);
        Kokkos::parallel_for(nrecv/size_exchange,f);
      }
    }

    k_count.modify_device();
    k_count.sync_host();
  }

  atomKK->modified(space,datamask_exchange);

  return k_count.view_host()(0);
}

/* ---------------------------------------------------------------------- */

uint64_t AtomVecKokkos::field2mask(std::string field)
{
  if (field == "id")
    return TAG_MASK;
  else if (field == "type")
    return TYPE_MASK;
  else if (field == "mask")
    return MASK_MASK;
  else if (field == "image")
    return IMAGE_MASK;
  else if (field == "x")
    return X_MASK;
  else if (field == "v")
    return V_MASK;
  else if (field == "f")
    return F_MASK;
  else if (field == "rmass")
    return RMASS_MASK;
  else if (field == "q")
    return Q_MASK;
  else if (field == "mu")
    return MU_MASK;
  else if (field == "mu3")
    return MU_MASK;
  else if (field == "radius")
    return RADIUS_MASK;
  else if (field == "omega")
    return OMEGA_MASK;
  else if (field == "torque")
    return TORQUE_MASK;
  else if (field == "molecule")
    return MOLECULE_MASK;
  else if (field == "nspecial")
    return SPECIAL_MASK;
  else if (field == "num_bond")
    return BOND_MASK;
  else if (field == "num_angle")
    return ANGLE_MASK;
  else if (field == "num_dihedral")
    return DIHEDRAL_MASK;
  else if (field == "num_improper")
    return IMPROPER_MASK;
  else if (field == "sp")
    return SP_MASK;
  else if (field == "fm")
    return FM_MASK;
  else if (field == "fm_long")
    return FML_MASK;
  else if (field == "rho") // conflicts with SPH package "rho"
    return DPDRHO_MASK;
  else if (field == "dpdTheta")
    return DPDTHETA_MASK;
  else if (field == "uCond")
    return UCOND_MASK;
  else if (field == "uMech")
    return UMECH_MASK;
  else if (field == "uChem")
    return UCHEM_MASK;
  else if (field == "uCG")
    return UCG_MASK;
  else if (field == "uCGnew")
    return UCGNEW_MASK;
  else if (field == "duChem")
    return DUCHEM_MASK;
  else
    return EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

int AtomVecKokkos::field2size(std::string field)
{
  if (field == "id") return 1;
  else if (field == "type") return 1;
  else if (field == "mask") return 1;
  else if (field == "image") return 1;
  else if (field == "x") return 3;
  else if (field == "v") return 3;
  else if (field == "f") return 3;
  else if (field == "rmass") return 1;
  else if (field == "q") return 1;
  else if (field == "mu") return 4;
  else if (field == "mu3") return 3;
  else if (field == "radius") return 1;
  else if (field == "omega") return 3;
  else if (field == "torque") return 3;
  else if (field == "molecule") return 1;
  else if (field == "special") return 3+atom->maxspecial;
  else if (field == "num_bond") return 1+2*atom->bond_per_atom;
  else if (field == "num_angle") return 1+4*atom->angle_per_atom;
  else if (field == "num_dihedral") return 1+5*atom->dihedral_per_atom;
  else if (field == "num_improper") return 1+5*atom->dihedral_per_atom;
  else if (field == "sp") return 4;
  else if (field == "fm") return 3;
  else if (field == "fm_long") return 3;
  else if (field == "rho") return 1;
  else if (field == "dpdTheta") return 1;
  else if (field == "uCond") return 1;
  else if (field == "uMech") return 1;
  else if (field == "uChem") return 1;
  else if (field == "uCG") return 1;
  else if (field == "uCGnew") return 1;
  else if (field == "duChem") return 1;
  else return 0;
}

/* ---------------------------------------------------------------------- */

void AtomVecKokkos::set_atom_masks()
{
  datamask_grow = EMPTY_MASK;
  for (int i = 0; i < default_grow.size(); i++)
    datamask_grow |= field2mask(default_grow[i]);
  for (int i = 0; i < ngrow; i++)
    datamask_grow |= field2mask(fields_grow[i]);

  datamask_comm = EMPTY_MASK;
  for (int i = 0; i < default_comm.size(); i++)
    datamask_comm |= field2mask(default_comm[i]);
  for (int i = 0; i < ncomm; i++)
    datamask_comm |= field2mask(fields_comm[i]);

  datamask_comm_vel = EMPTY_MASK;
  for (int i = 0; i < default_comm_vel.size(); i++)
    datamask_comm_vel |= field2mask(default_comm_vel[i]);
  for (int i = 0; i < ncomm_vel; i++)
    datamask_comm_vel |= field2mask(fields_comm_vel[i]);

  datamask_reverse = EMPTY_MASK;
  for (int i = 0; i < default_reverse.size(); i++)
    datamask_reverse |= field2mask(default_reverse[i]);
  for (int i = 0; i < nreverse; i++)
    datamask_reverse |= field2mask(fields_reverse[i]);

  datamask_border = EMPTY_MASK;
  for (int i = 0; i < default_border.size(); i++)
    datamask_border |= field2mask(default_border[i]);
  for (int i = 0; i < nborder; i++)
    datamask_border |= field2mask(fields_border[i]);

  datamask_border_vel = EMPTY_MASK;
  for (int i = 0; i < default_border_vel.size(); i++)
    datamask_border_vel |= field2mask(default_border_vel[i]);
  for (int i = 0; i < nborder_vel; i++)
    datamask_border_vel |= field2mask(fields_border_vel[i]);

  datamask_exchange = EMPTY_MASK;
  for (int i = 0; i < default_exchange.size(); i++)
    datamask_exchange |= field2mask(default_exchange[i]);
  for (int i = 0; i < nexchange; i++)
    datamask_exchange |= field2mask(fields_exchange[i]);
}

/* ---------------------------------------------------------------------- */

void AtomVecKokkos::set_size_exchange()
{
  size_exchange_default = 1; // 1 to store buffer length
  for (int i = 0; i < default_exchange.size(); i++)
    size_exchange_default += field2size(default_exchange[i]);

  size_exchange = size_exchange_default;

  for (int i = 0; i < nexchange; i++)
    size_exchange += field2size(fields_exchange[i]);
}
