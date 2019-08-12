/* -*- c++ -*- ----------------------------------------------------------

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(eam/kk,PairEAMKokkos<LMPDeviceType>)
PairStyle(eam/kk/device,PairEAMKokkos<LMPDeviceType>)
PairStyle(eam/kk/host,PairEAMKokkos<LMPHostType>)

#else

#ifndef LMP_PAIR_EAM_KOKKOS_H
#define LMP_PAIR_EAM_KOKKOS_H

#include "kokkos_base.h"
#include "pair_kokkos.h"
#include "pair_eam.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

struct TagPairEAMPackForwardComm{};
struct TagPairEAMUnpackForwardComm{};
struct TagPairEAMInitialize{};

template<int NEIGHFLAG, int NEWTON_PAIR>
struct TagPairEAMKernelA{};

template<int EFLAG>
struct TagPairEAMKernelB{};

template<int EFLAG>
struct TagPairEAMKernelAB{};

template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairEAMKernelC{};

template<class DeviceType>
class PairEAMKokkos : public PairEAM, public KokkosBase {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairEAMKokkos(class LAMMPS *);
  virtual ~PairEAMKokkos();
  void compute(int, int);
  void init_style();
  void *extract(const char *, int &) { return NULL; }

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMPackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMUnpackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMInitialize, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMKernelA<NEIGHFLAG,NEWTON_PAIR>, const int&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMKernelB<EFLAG>, const int&, EV_FLOAT&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMKernelB<EFLAG>, const int&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMKernelAB<EFLAG>, const int&, EV_FLOAT&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMKernelAB<EFLAG>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEAMKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &epair, const F_FLOAT &fpair, const F_FLOAT &delx,
                  const F_FLOAT &dely, const F_FLOAT &delz) const;

  int pack_forward_comm_kokkos(int, DAT::tdual_int_2d, int, DAT::tdual_xfloat_1d&,
                       int, int *);
  void unpack_forward_comm_kokkos(int, int, DAT::tdual_xfloat_1d&);
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);

 protected:
  void cleanup_copy();

  typename AT::t_x_array_randomread x;
  typename AT::t_f_array f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_tagint_1d tag;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename ArrayTypes<DeviceType>::t_efloat_1d d_eatom;
  typename ArrayTypes<DeviceType>::t_virial_array d_vatom;

  int need_dup;
  Kokkos::Experimental::ScatterView<F_FLOAT*, typename DAT::t_ffloat_1d::array_layout,DeviceType,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterDuplicated> dup_rho;
  Kokkos::Experimental::ScatterView<F_FLOAT*[3], typename DAT::t_f_array::array_layout,DeviceType,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterDuplicated> dup_f;
  Kokkos::Experimental::ScatterView<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,DeviceType,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterDuplicated> dup_eatom;
  Kokkos::Experimental::ScatterView<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,DeviceType,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterDuplicated> dup_vatom;
  Kokkos::Experimental::ScatterView<F_FLOAT*, typename DAT::t_ffloat_1d::array_layout,DeviceType,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterNonDuplicated> ndup_rho;
  Kokkos::Experimental::ScatterView<F_FLOAT*[3], typename DAT::t_f_array::array_layout,DeviceType,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterNonDuplicated> ndup_f;
  Kokkos::Experimental::ScatterView<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,DeviceType,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterNonDuplicated> ndup_eatom;
  Kokkos::Experimental::ScatterView<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,DeviceType,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterNonDuplicated> ndup_vatom;

  DAT::tdual_ffloat_1d k_rho;
  DAT::tdual_ffloat_1d k_fp;
  typename AT::t_ffloat_1d d_rho;
  typename AT::t_ffloat_1d d_fp;
  HAT::t_ffloat_1d h_rho;
  HAT::t_ffloat_1d h_fp;

  typename AT::t_int_1d_randomread d_type2frho;
  typename AT::t_int_2d_randomread d_type2rhor;
  typename AT::t_int_2d_randomread d_type2z2r;

  typedef Kokkos::DualView<F_FLOAT**[7],Kokkos::LayoutRight,DeviceType> tdual_ffloat_2d_n7;
  typedef typename tdual_ffloat_2d_n7::t_dev_const_randomread t_ffloat_2d_n7_randomread;
  typedef typename tdual_ffloat_2d_n7::t_host t_host_ffloat_2d_n7;

  t_ffloat_2d_n7_randomread d_frho_spline;
  t_ffloat_2d_n7_randomread d_rhor_spline;
  t_ffloat_2d_n7_randomread d_z2r_spline;
  void interpolate(int, double, double *, t_host_ffloat_2d_n7, int);

  void file2array();
  void array2spline();

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;
  //NeighListKokkos<DeviceType> k_list;

  int iswap;
  int first;
  typename AT::t_int_2d d_sendlist;
  typename AT::t_xfloat_1d_um v_buf;


  int neighflag,newton_pair;
  int nlocal,nall,eflag,vflag;

  friend void pair_virial_fdotr_compute<PairEAMKokkos>(PairEAMKokkos*);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot use chosen neighbor list style with pair eam/kk

That style is not supported by Kokkos.

*/
