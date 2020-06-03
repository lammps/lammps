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

#include "atom_vec_kokkos.h"
#include "atom_kokkos.h"
#include "comm_kokkos.h"
#include "domain.h"
#include "atom_masks.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecKokkos::AtomVecKokkos(LAMMPS *lmp) : AtomVec(lmp)
{
  kokkosable = 1;
  buffer = NULL;
  buffer_size = 0;

  no_comm_vel_flag = 0;
  no_border_vel_flag = 1;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space,int PBC_FLAG,int TRICLINIC>
struct AtomVecKokkos_PackComm {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;

  typename AT::t_float_1d_3_randomread _x;
  typename AT::t_float_2d_um _buf;
  typename AT::t_int_2d_const _list;
  const int _iswap;
  KK_FLOAT _xprd,_yprd,_zprd,_xy,_xz,_yz;
  KK_FLOAT _pbc[6];

  AtomVecKokkos_PackComm(
      const typename DAT::tdual_float_1d_3 &x,
      const typename DAT::tdual_float_2d &buf,
      const typename DAT::tdual_int_2d &list,
      const int & iswap,
      const KK_FLOAT &xprd, const KK_FLOAT &yprd, const KK_FLOAT &zprd,
      const KK_FLOAT &xy, const KK_FLOAT &xz, const KK_FLOAT &yz, const int* const pbc):
      _x(DualViewHelper<Space>::view(x)),_list(DualViewHelper<Space>::view(list)),_iswap(iswap),
      _xprd(xprd),_yprd(yprd),_zprd(zprd),
      _xy(xy),_xz(xz),_yz(yz) {
        const size_t maxsend = (DualViewHelper<Space>::view(buf).extent(0)*DualViewHelper<Space>::view(buf).extent(1))/3;
        const size_t elements = 3;
        buffer_view<Space>(_buf,buf,maxsend,elements);
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

int AtomVecKokkos::pack_comm_kokkos(const int &n,
                                          const DAT::tdual_int_2d &list,
                                          const int & iswap,
                                          const DAT::tdual_float_2d &buf,
                                          const int &pbc_flag,
                                          const int* const pbc)
{
  // Check whether to always run forward communication on the host
  // Choose correct forward PackComm kernel

  if(commKK->forward_comm_on_host) {
    sync(Host,X_MASK);
    if(pbc_flag) {
      if(domain->triclinic) {
        struct AtomVecKokkos_PackComm<Host,1,1> f(atomKK->k_x,buf,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecKokkos_PackComm<Host,1,0> f(atomKK->k_x,buf,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      }
    } else {
      if(domain->triclinic) {
        struct AtomVecKokkos_PackComm<Host,0,1> f(atomKK->k_x,buf,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecKokkos_PackComm<Host,0,0> f(atomKK->k_x,buf,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      }
    }
  } else {
    sync(Device,X_MASK);
    if(pbc_flag) {
      if(domain->triclinic) {
        struct AtomVecKokkos_PackComm<Device,1,1> f(atomKK->k_x,buf,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecKokkos_PackComm<Device,1,0> f(atomKK->k_x,buf,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      }
    } else {
      if(domain->triclinic) {
        struct AtomVecKokkos_PackComm<Device,0,1> f(atomKK->k_x,buf,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecKokkos_PackComm<Device,0,0> f(atomKK->k_x,buf,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
        Kokkos::parallel_for(n,f);
      }
    }
  }

        return n*size_forward;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space,int PBC_FLAG,int TRICLINIC>
struct AtomVecKokkos_PackCommSelf {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;

  typename AT::t_float_1d_3_randomread _x;
  typename AT::t_float_1d_3 _xw;
  int _nfirst;
  typename AT::t_int_2d_const _list;
  const int _iswap;
  KK_FLOAT _xprd,_yprd,_zprd,_xy,_xz,_yz;
  KK_FLOAT _pbc[6];

  AtomVecKokkos_PackCommSelf(
      const typename DAT::tdual_float_1d_3 &x,
      const int &nfirst,
      const typename DAT::tdual_int_2d &list,
      const int & iswap,
      const KK_FLOAT &xprd, const KK_FLOAT &yprd, const KK_FLOAT &zprd,
      const KK_FLOAT &xy, const KK_FLOAT &xz, const KK_FLOAT &yz, const int* const pbc):
      _x(DualViewHelper<Space>::view(x)),_xw(DualViewHelper<Space>::view(x)),_nfirst(nfirst),_list(DualViewHelper<Space>::view(list)),_iswap(iswap),
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

int AtomVecKokkos::pack_comm_self(const int &n, const DAT::tdual_int_2d &list, const int & iswap,
                                        const int nfirst, const int &pbc_flag, const int* const pbc) {
  if(commKK->forward_comm_on_host) {
    sync(Host,X_MASK);
    modified(Host,X_MASK);
    if(pbc_flag) {
      if(domain->triclinic) {
      struct AtomVecKokkos_PackCommSelf<Host,1,1> f(atomKK->k_x,nfirst,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      } else {
      struct AtomVecKokkos_PackCommSelf<Host,1,0> f(atomKK->k_x,nfirst,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      }
    } else {
      if(domain->triclinic) {
      struct AtomVecKokkos_PackCommSelf<Host,0,1> f(atomKK->k_x,nfirst,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      } else {
      struct AtomVecKokkos_PackCommSelf<Host,0,0> f(atomKK->k_x,nfirst,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      }
    }
  } else {
    sync(Device,X_MASK);
    modified(Device,X_MASK);
    if(pbc_flag) {
      if(domain->triclinic) {
      struct AtomVecKokkos_PackCommSelf<Device,1,1> f(atomKK->k_x,nfirst,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      } else {
      struct AtomVecKokkos_PackCommSelf<Device,1,0> f(atomKK->k_x,nfirst,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      }
    } else {
      if(domain->triclinic) {
      struct AtomVecKokkos_PackCommSelf<Device,0,1> f(atomKK->k_x,nfirst,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      } else {
      struct AtomVecKokkos_PackCommSelf<Device,0,0> f(atomKK->k_x,nfirst,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc);
      Kokkos::parallel_for(n,f);
      }
    }
  }
        return n*3;
}


/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space,int TRICLINIC>
struct AtomVecKokkos_PackCommSelfFused {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;

  typename AT::t_float_1d_3_randomread _x;
  typename AT::t_float_1d_3 _xw;
  typename AT::t_int_2d_const _list;
  typename AT::t_int_2d_const _pbc;
  typename AT::t_int_1d_const _pbc_flag;
  typename AT::t_int_1d_const _firstrecv;
  typename AT::t_int_1d_const _sendnum_scan;
  typename AT::t_int_1d_const _g2l;
  KK_FLOAT _xprd,_yprd,_zprd,_xy,_xz,_yz;

  AtomVecKokkos_PackCommSelfFused(
      const typename DAT::tdual_float_1d_3 &x,
      const typename DAT::tdual_int_2d &list,
      const typename DAT::tdual_int_2d &pbc,
      const typename DAT::tdual_int_1d &pbc_flag,
      const typename DAT::tdual_int_1d &firstrecv,
      const typename DAT::tdual_int_1d &sendnum_scan,
      const typename DAT::tdual_int_1d &g2l,
      const KK_FLOAT &xprd, const KK_FLOAT &yprd, const KK_FLOAT &zprd,
      const KK_FLOAT &xy, const KK_FLOAT &xz, const KK_FLOAT &yz):
      _x(DualViewHelper<Space>::view(x)),_xw(DualViewHelper<Space>::view(x)),
      _list(DualViewHelper<Space>::view(list)),
      _pbc(DualViewHelper<Space>::view(pbc)),
      _pbc_flag(DualViewHelper<Space>::view(pbc_flag)),
      _firstrecv(DualViewHelper<Space>::view(firstrecv)),
      _sendnum_scan(DualViewHelper<Space>::view(sendnum_scan)),
      _g2l(DualViewHelper<Space>::view(g2l)),
      _xprd(xprd),_yprd(yprd),_zprd(zprd),
      _xy(xy),_xz(xz),_yz(yz) {};

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
        _xw(i+_nfirst,0) = _x(j,0);
        _xw(i+_nfirst,1) = _x(j,1);
        _xw(i+_nfirst,2) = _x(j,2);
    } else {
      if (TRICLINIC == 0) {
        _xw(i+_nfirst,0) = _x(j,0) + _pbc(ii,0)*_xprd;
        _xw(i+_nfirst,1) = _x(j,1) + _pbc(ii,1)*_yprd;
        _xw(i+_nfirst,2) = _x(j,2) + _pbc(ii,2)*_zprd;
      } else {
        _xw(i+_nfirst,0) = _x(j,0) + _pbc(ii,0)*_xprd + _pbc(ii,5)*_xy + _pbc(ii,4)*_xz;
        _xw(i+_nfirst,1) = _x(j,1) + _pbc(ii,1)*_yprd + _pbc(ii,3)*_yz;
        _xw(i+_nfirst,2) = _x(j,2) + _pbc(ii,2)*_zprd;
      }
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecKokkos::pack_comm_self_fused(const int &n, const DAT::tdual_int_2d &list, const DAT::tdual_int_1d &sendnum_scan,
                                         const DAT::tdual_int_1d &firstrecv, const DAT::tdual_int_1d &pbc_flag, const DAT::tdual_int_2d &pbc,
                                         const DAT::tdual_int_1d &g2l) {
  if(commKK->forward_comm_on_host) {
    sync(Host,X_MASK);
    modified(Host,X_MASK);
    if(domain->triclinic) {
    struct AtomVecKokkos_PackCommSelfFused<Host,1> f(atomKK->k_x,list,pbc,pbc_flag,firstrecv,sendnum_scan,g2l,
        domain->xprd,domain->yprd,domain->zprd,
        domain->xy,domain->xz,domain->yz);
    Kokkos::parallel_for(n,f);
    } else {
    struct AtomVecKokkos_PackCommSelfFused<Host,0> f(atomKK->k_x,list,pbc,pbc_flag,firstrecv,sendnum_scan,g2l,
        domain->xprd,domain->yprd,domain->zprd,
        domain->xy,domain->xz,domain->yz);
    Kokkos::parallel_for(n,f);
    }
  } else {
    sync(Device,X_MASK);
    modified(Device,X_MASK);
    if(domain->triclinic) {
    struct AtomVecKokkos_PackCommSelfFused<Device,1> f(atomKK->k_x,list,pbc,pbc_flag,firstrecv,sendnum_scan,g2l,
        domain->xprd,domain->yprd,domain->zprd,
        domain->xy,domain->xz,domain->yz);
    Kokkos::parallel_for(n,f);
    } else {
    struct AtomVecKokkos_PackCommSelfFused<Device,0> f(atomKK->k_x,list,pbc,pbc_flag,firstrecv,sendnum_scan,g2l,
        domain->xprd,domain->yprd,domain->zprd,
        domain->xy,domain->xz,domain->yz);
    Kokkos::parallel_for(n,f);
    }
  }
  return n*3;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
struct AtomVecKokkos_UnpackComm {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;

  typename AT::t_float_1d_3 _x;
  typename AT::t_float_2d_const _buf;
  int _first;

  AtomVecKokkos_UnpackComm(
      const typename DAT::tdual_float_1d_3 &x,
      const typename DAT::tdual_float_2d &buf,
      const int& first):_x(DualViewHelper<Space>::view(x)),_buf(DualViewHelper<Space>::view(buf)),
                        _first(first) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
      _x(i+_first,0) = _buf(i,0);
      _x(i+_first,1) = _buf(i,1);
      _x(i+_first,2) = _buf(i,2);
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecKokkos::unpack_comm_kokkos(const int &n, const int &first,
    const DAT::tdual_float_2d &buf ) {
  if(commKK->forward_comm_on_host) {
    sync(Host,X_MASK);
    modified(Host,X_MASK);
    struct AtomVecKokkos_UnpackComm<Host> f(atomKK->k_x,buf,first);
    Kokkos::parallel_for(n,f);
  } else {
    sync(Device,X_MASK);
    modified(Device,X_MASK);
    struct AtomVecKokkos_UnpackComm<Device> f(atomKK->k_x,buf,first);
    Kokkos::parallel_for(n,f);
  }
}


/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space,int PBC_FLAG,int TRICLINIC,int DEFORM_VREMAP>
struct AtomVecKokkos_PackCommVel {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;

  typename AT::t_float_1d_3_randomread _x;
  typename AT::t_int_1d _mask;
  typename AT::t_float_1d_3 _v;
  typename AT::t_float_2d_um _buf;
  typename AT::t_int_2d_const _list;
  const int _iswap;
  double _xprd,_yprd,_zprd,_xy,_xz,_yz;
  double _pbc[6];
  double _h_rate[6];
  const int _deform_vremap;

  AtomVecKokkos_PackCommVel(
    const typename DAT::tdual_float_1d_3 &x,
    const typename DAT::tdual_int_1d &mask,
    const typename DAT::tdual_float_1d_3 &v,
    const typename DAT::tdual_float_2d &buf,
    const typename DAT::tdual_int_2d &list,
    const int &iswap,
    const double &xprd, const double &yprd, const double &zprd,
    const double &xy, const double &xz, const double &yz, const int* const pbc,
    const double * const h_rate,
    const int &deform_vremap):
    _x(DualViewHelper<Space>::view(x)),
    _mask(DualViewHelper<Space>::view(mask)),
    _v(DualViewHelper<Space>::view(v)),
    _list(DualViewHelper<Space>::view(list)),_iswap(iswap),
    _xprd(xprd),_yprd(yprd),_zprd(zprd),
    _xy(xy),_xz(xz),_yz(yz),
    _deform_vremap(deform_vremap)
  {
    const size_t elements = 6;
    const int maxsend = (DualViewHelper<Space>::view(buf).extent(0)*DualViewHelper<Space>::view(buf).extent(1))/elements;
    _buf = typename AT::t_float_2d_um(DualViewHelper<Space>::view(buf).data(),maxsend,elements);
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
      _buf(i,3) = _v(j,0);
      _buf(i,4) = _v(j,1);
      _buf(i,5) = _v(j,2);
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
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecKokkos::pack_comm_vel_kokkos(
  const int &n,
  const DAT::tdual_int_2d &list,
  const int & iswap,
  const DAT::tdual_float_2d &buf,
  const int &pbc_flag,
  const int* const pbc)
{
  if(commKK->forward_comm_on_host) {
    sync(Host,X_MASK|V_MASK);
    if (pbc_flag) {
      if (deform_vremap) {
        if (domain->triclinic) {
          struct AtomVecKokkos_PackCommVel<Host,1,1,1> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_v,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackCommVel<Host,1,0,1> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_v,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        }
      } else {
        if(domain->triclinic) {
          struct AtomVecKokkos_PackCommVel<Host,1,1,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_v,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackCommVel<Host,1,0,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_v,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        }
      }
    } else {
      if(domain->triclinic) {
        struct AtomVecKokkos_PackCommVel<Host,0,1,0> f(
          atomKK->k_x,atomKK->k_mask,
          atomKK->k_v,
          buf,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecKokkos_PackCommVel<Host,0,0,0> f(
          atomKK->k_x,atomKK->k_mask,
          atomKK->k_v,
          buf,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
        Kokkos::parallel_for(n,f);
      }
    }
  } else {
    sync(Device,X_MASK|V_MASK);
    if(pbc_flag) {
      if(deform_vremap) {
        if(domain->triclinic) {
          struct AtomVecKokkos_PackCommVel<Device,1,1,1> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_v,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackCommVel<Device,1,0,1> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_v,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        }
      } else {
        if(domain->triclinic) {
          struct AtomVecKokkos_PackCommVel<Device,1,1,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_v,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecKokkos_PackCommVel<Device,1,0,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_v,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        }
      }
    } else {
      if(domain->triclinic) {
        struct AtomVecKokkos_PackCommVel<Device,0,1,0> f(
          atomKK->k_x,atomKK->k_mask,
          atomKK->k_v,
          buf,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
        Kokkos::parallel_for(n,f);
      } else {
        struct AtomVecKokkos_PackCommVel<Device,0,0,0> f(
          atomKK->k_x,atomKK->k_mask,
          atomKK->k_v,
          buf,list,iswap,
          domain->xprd,domain->yprd,domain->zprd,
          domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
        Kokkos::parallel_for(n,f);
      }
    }
  }
  return n*6;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
struct AtomVecKokkos_UnpackCommVel {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;

  typename AT::t_float_1d_3 _x;
  typename AT::t_float_1d_3 _v;
  typename AT::t_float_2d_const _buf;
  int _first;

  AtomVecKokkos_UnpackCommVel(
    const typename DAT::tdual_float_1d_3 &x,
    const typename DAT::tdual_float_1d_3 &v,
    const typename DAT::tdual_float_2d &buf,
    const int& first):
    _x(DualViewHelper<Space>::view(x)),
    _v(DualViewHelper<Space>::view(v)),
    _first(first)
  {
    const size_t elements = 6;
    const int maxsend = (DualViewHelper<Space>::view(buf).extent(0)*DualViewHelper<Space>::view(buf).extent(1))/elements;
    buffer_view<Space>(_buf,buf,maxsend,elements);
  };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    _x(i+_first,0) = _buf(i,0);
    _x(i+_first,1) = _buf(i,1);
    _x(i+_first,2) = _buf(i,2);
    _v(i+_first,0) = _buf(i,3);
    _v(i+_first,1) = _buf(i,4);
    _v(i+_first,2) = _buf(i,5);
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecKokkos::unpack_comm_vel_kokkos(const int &n, const int &first,
    const DAT::tdual_float_2d &buf ) {
  if(commKK->forward_comm_on_host) {
    sync(Host,X_MASK|V_MASK);
    modified(Host,X_MASK|V_MASK);
    struct AtomVecKokkos_UnpackCommVel<Host> f(atomKK->k_x,atomKK->k_v,buf,first);
    Kokkos::parallel_for(n,f);
  } else {
    sync(Device,X_MASK|V_MASK);
    modified(Device,X_MASK|V_MASK);
    struct AtomVecKokkos_UnpackCommVel<Device> f(atomKK->k_x,atomKK->k_v,buf,first);
    Kokkos::parallel_for(n,f);
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecKokkos::pack_comm(int n, int *list, double *buf,
                             int pbc_flag, int *pbc)
{
  int i,j,m;
  KK_FLOAT dx,dy,dz;

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
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecKokkos::pack_comm_vel(int n, int *list, double *buf,
                                 int pbc_flag, int *pbc)
{
  int i,j,m;
  KK_FLOAT dx,dy,dz,dvx,dvy,dvz;

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
        if (atom->mask[i] & deform_groupbit) {
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
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecKokkos::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    h_x(i,0) = buf[m++];
    h_x(i,1) = buf[m++];
    h_x(i,2) = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecKokkos::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    h_x(i,0) = buf[m++];
    h_x(i,1) = buf[m++];
    h_x(i,2) = buf[m++];
    h_v(i,0) = buf[m++];
    h_v(i,1) = buf[m++];
    h_v(i,2) = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
struct AtomVecKokkos_PackReverse {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;

  typename AT::t_float_1d_3_randomread _f;
  typename AT::t_float_2d _buf;
  int _first;

  AtomVecKokkos_PackReverse(
      const typename DAT::tdual_float_1d_3 &f,
      const typename DAT::tdual_float_2d &buf,
      const int& first):_f(DualViewHelper<Space>::view(f)),_buf(DualViewHelper<Space>::view(buf)),
                        _first(first) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    _buf(i,0) = _f(i+_first,0);
    _buf(i,1) = _f(i+_first,1);
    _buf(i,2) = _f(i+_first,2);
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecKokkos::pack_reverse_kokkos(const int &n, const int &first,
    const DAT::tdual_float_2d &buf ) {
  if(commKK->reverse_comm_on_host) {
    sync(Host,F_MASK);
    struct AtomVecKokkos_PackReverse<Host> f(atomKK->k_f,buf,first);
    Kokkos::parallel_for(n,f);
  } else {
    sync(Device,F_MASK);
    struct AtomVecKokkos_PackReverse<Device> f(atomKK->k_f,buf,first);
    Kokkos::parallel_for(n,f);
  }

  return n*size_reverse;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
struct AtomVecKokkos_UnPackReverseSelf {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;

  typename AT::t_float_1d_3_randomread _f;
  typename AT::t_float_1d_3 _fw;
  int _nfirst;
  typename AT::t_int_2d_const _list;
  const int _iswap;

  AtomVecKokkos_UnPackReverseSelf(
      const typename DAT::tdual_float_1d_3 &f,
      const int &nfirst,
      const typename DAT::tdual_int_2d &list,
      const int & iswap):
      _f(DualViewHelper<Space>::view(f)),_fw(DualViewHelper<Space>::view(f)),_nfirst(nfirst),_list(DualViewHelper<Space>::view(list)),_iswap(iswap) {
  };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(_iswap,i);
    _fw(j,0) += _f(i+_nfirst,0);
    _fw(j,1) += _f(i+_nfirst,1);
    _fw(j,2) += _f(i+_nfirst,2);
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecKokkos::unpack_reverse_self(const int &n, const DAT::tdual_int_2d &list, const int & iswap,
                                        const int nfirst) {
  if(commKK->reverse_comm_on_host) {
    sync(Host,F_MASK);
    struct AtomVecKokkos_UnPackReverseSelf<Host> f(atomKK->k_f,nfirst,list,iswap);
    Kokkos::parallel_for(n,f);
    modified(Host,F_MASK);
  } else {
    sync(Device,F_MASK);
    struct AtomVecKokkos_UnPackReverseSelf<Device> f(atomKK->k_f,nfirst,list,iswap);
    Kokkos::parallel_for(n,f);
    modified(Device,F_MASK);
  }
  return n*3;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
struct AtomVecKokkos_UnPackReverse {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;

  typename AT::t_float_1d_3 _f;
  typename AT::t_float_2d_const _buf;
  typename AT::t_int_2d_const _list;
  const int _iswap;

  AtomVecKokkos_UnPackReverse(
      const typename DAT::tdual_float_1d_3 &f,
      const typename DAT::tdual_float_2d &buf,
      const typename DAT::tdual_int_2d &list,
      const int & iswap):
      _f(DualViewHelper<Space>::view(f)),_list(DualViewHelper<Space>::view(list)),_iswap(iswap) {
        const size_t maxsend = (DualViewHelper<Space>::view(buf).extent(0)*DualViewHelper<Space>::view(buf).extent(1))/3;
        const size_t elements = 3;
        buffer_view<Space>(_buf,buf,maxsend,elements);
  };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(_iswap,i);
    _f(j,0) += _buf(i,0);
    _f(j,1) += _buf(i,1);
    _f(j,2) += _buf(i,2);
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecKokkos::unpack_reverse_kokkos(const int &n,
                                          const DAT::tdual_int_2d &list,
                                          const int & iswap,
                                          const DAT::tdual_float_2d &buf)
{
  // Check whether to always run reverse communication on the host
  // Choose correct reverse UnPackReverse kernel

  if(commKK->reverse_comm_on_host) {
    struct AtomVecKokkos_UnPackReverse<Host> f(atomKK->k_f,buf,list,iswap);
    Kokkos::parallel_for(n,f);
    modified(Host,F_MASK);
  } else {
    struct AtomVecKokkos_UnPackReverse<Device> f(atomKK->k_f,buf,list,iswap);
    Kokkos::parallel_for(n,f);
    modified(Device,F_MASK);
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecKokkos::pack_reverse(int n, int first, double *buf)
{
  if(n > 0)
    sync(Host,F_MASK);

  int m = 0;
  const int last = first + n;
  for (int i = first; i < last; i++) {
    buf[m++] = h_f(i,0);
    buf[m++] = h_f(i,1);
    buf[m++] = h_f(i,2);
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecKokkos::unpack_reverse(int n, int *list, double *buf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    const int j = list[i];
    h_f(j,0) += buf[m++];
    h_f(j,1) += buf[m++];
    h_f(j,2) += buf[m++];
  }

  if(n > 0)
    modified(Host,F_MASK);
}

/* ----------------------------------------------------------------------
 *    unpack one line from Velocities section of data file
 *    ------------------------------------------------------------------------- */

void AtomVecKokkos::data_vel(int m, char **values)
{
  double **v = atom->v;
  v[m][0] = utils::numeric(FLERR,values[0],true,lmp);
  v[m][1] = utils::numeric(FLERR,values[1],true,lmp);
  v[m][2] = utils::numeric(FLERR,values[2],true,lmp);

  modified(Host,V_MASK);
}

/* ----------------------------------------------------------------------
 *    pack velocity info for data file
 *    ------------------------------------------------------------------------- */

void AtomVecKokkos::pack_vel(double **buf)
{
  double **v = atom->v;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  sync(Host,V_MASK|TAG_MASK);

  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = v[i][0];
    buf[i][2] = v[i][1];
    buf[i][3] = v[i][2];
  }
}

/* ----------------------------------------------------------------------
 *    write velocity info to data file
 *    ------------------------------------------------------------------------- */

void AtomVecKokkos::write_vel(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT " %-1.16e %-1.16e %-1.16e\n",
            (tagint) ubuf(buf[i][0]).i,buf[i][1],buf[i][2],buf[i][3]);
}

