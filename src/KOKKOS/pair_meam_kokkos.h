/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(meam/kk,PairMEAMKokkos<LMPDeviceType>);
PairStyle(meam/kk/device,PairMEAMKokkos<LMPDeviceType>);
PairStyle(meam/kk/host,PairMEAMKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_MEAM_KOKKOS_H
#define LMP_PAIR_MEAM_KOKKOS_H

#include "kokkos_base.h"
#include "pair_kokkos.h"
#include "pair_meam.h"
#include "meam_kokkos.h"

namespace LAMMPS_NS {

struct TagPairMEAMNeighStrip{};
struct TagPairMEAMOffsets{};
struct TagPairMEAMPackForwardComm{};
struct TagPairMEAMUnpackForwardComm{};
struct TagPairMEAMPackReverseComm{};
struct TagPairMEAMUnpackReverseComm{};

template<class DeviceType>
class MEAMKokkos;

template<class DeviceType>
class PairMEAMKokkos : public PairMEAM, public KokkosBase {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef int value_type;

  PairMEAMKokkos(class LAMMPS *);
  ~PairMEAMKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMEAMPackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMEAMUnpackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMEAMPackReverseComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMEAMUnpackReverseComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMEAMNeighStrip,  const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMEAMOffsets,  const int, int&) const;

  int pack_forward_comm_kokkos(int, DAT::tdual_int_1d, DAT::tdual_xfloat_1d&,
                               int, int *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm_kokkos(int, int, DAT::tdual_xfloat_1d&) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm_kokkos(int, int, DAT::tdual_xfloat_1d&) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm_kokkos(int, DAT::tdual_int_1d,
                                  DAT::tdual_xfloat_1d&) override;
  void unpack_reverse_comm(int, int *, double *) override;

 protected:
  class MEAMKokkos<DeviceType> *meam_inst_kk;
  typename AT::t_x_array x;
  typename AT::t_f_array f;
  typename AT::t_int_1d type;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename AT::t_efloat_1d d_eatom;
  typename AT::t_virial_array d_vatom;

  typename AT::t_int_1d d_offset;

  DAT::tdual_int_1d k_map;
  typename AT::t_int_1d d_map;
  typename AT::t_int_2d d_scale;
  typename AT::t_int_1d d_ilist_half;
  typename AT::t_int_1d d_numneigh_half;
  typename AT::t_neighbors_2d d_neighbors_half;
  typename AT::t_int_1d d_numneigh_full;
  typename AT::t_neighbors_2d d_neighbors_full;
  typename AT::t_int_1d d_sendlist;
  typename AT::t_xfloat_1d_um v_buf;

  int first;
  int neighflag,nlocal,nall,eflag,vflag;

  typename ArrayTypes<DeviceType>::t_ffloat_1d d_rho, d_rho0, d_rho1, d_rho2, d_rho3, d_frhop;
  typename ArrayTypes<DeviceType>::t_ffloat_1d d_gamma, d_dgamma1, d_dgamma2, d_dgamma3, d_arho2b;
  typename ArrayTypes<DeviceType>::t_ffloat_2d d_arho1, d_arho2, d_arho3, d_arho3b, d_t_ave, d_tsq_ave;
  // msmeam params
  typename ArrayTypes<DeviceType>::t_ffloat_1d d_arho2mb;
  typename ArrayTypes<DeviceType>::t_ffloat_2d d_arho1m, d_arho2m, d_arho3m, d_arho3mb;

  void update_meam_views();

  friend void pair_virial_fdotr_compute<PairMEAMKokkos>(PairMEAMKokkos*);
};

}
#endif
#endif

