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
PairStyle(multi/lucy/rx/kk,PairMultiLucyRXKokkos<LMPDeviceType>);
PairStyle(multi/lucy/rx/kk/device,PairMultiLucyRXKokkos<LMPDeviceType>);
PairStyle(multi/lucy/rx/kk/host,PairMultiLucyRXKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_MULTI_LUCY_RX_KOKKOS_H
#define LMP_PAIR_MULTI_LUCY_RX_KOKKOS_H


#include "pair_multi_lucy_rx.h"
#include "pair_kokkos.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagPairMultiLucyRXPackForwardComm{};
struct TagPairMultiLucyRXUnpackForwardComm{};

struct TagPairMultiLucyRXgetMixingWeights{};

template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, int TABSTYLE>
struct TagPairMultiLucyRXCompute{};

struct TagPairMultiLucyRXZero{};

template<int NEIGHFLAG, int NEWTON_PAIR, bool ONE_TYPE>
struct TagPairMultiLucyRXComputeLocalDensity{};

template<class DeviceType>
class PairMultiLucyRXKokkos : public PairMultiLucyRX, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairMultiLucyRXKokkos(class LAMMPS *);
  ~PairMultiLucyRXKokkos() override;

  void compute(int, int) override;
  void settings(int, char **) override;

  template<int TABSTYLE>
  void compute_style(int, int);

  void init_style() override;
  int pack_forward_comm_kokkos(int, DAT::tdual_int_1d, DAT::tdual_double_1d&,
                               int, int *) override;
  void unpack_forward_comm_kokkos(int, int, DAT::tdual_double_1d&) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  void computeLocalDensity();

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMultiLucyRXPackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMultiLucyRXUnpackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMultiLucyRXgetMixingWeights, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, int TABSTYLE>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMultiLucyRXCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG,TABSTYLE>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, int TABSTYLE>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMultiLucyRXCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG,TABSTYLE>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMultiLucyRXZero, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, bool ONE_TYPE>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMultiLucyRXComputeLocalDensity<NEIGHFLAG,NEWTON_PAIR,ONE_TYPE>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                  const KK_FLOAT &dely, const KK_FLOAT &delz) const;

 private:
  int nlocal;
  int neighflag;
  int eflag,vflag;

  double cutsq_type11;
  double rcut_type11;
  double factor_type11;

  enum{LOOKUP,LINEAR,SPLINE,BITMAP};

  //struct Table {
  //  int ninput,rflag,fpflag,match;
  //  double rlo,rhi,fplo,fphi,cut;
  //  double *rfile,*efile,*ffile;
  //  double *e2file,*f2file;
  //  double innersq,delta,invdelta,deltasq6;
  //  double *rsq,*drsq,*e,*de,*f,*df,*e2,*f2;
  //};

  /*struct TableDeviceConst {
    typename AT::t_int_2d_lr_randomread tabindex;
    typename AT::t_double_1d_randomread innersq,invdelta;
    typename AT::t_double_2d_lr_randomread rsq,e,de,f,df;
  };*/
 //Its faster not to use texture fetch if the number of tables is less than 32!
  struct TableDeviceConst {
    typename AT::t_int_2d_lr tabindex;
    typename AT::t_double_1d innersq,invdelta;
    typename AT::t_double_2d_lr_randomread rsq,e,de,f,df;
  };

  struct TableDevice {
    typename AT::t_int_2d_lr tabindex;
    typename AT::t_double_1d innersq,invdelta;
    typename AT::t_double_2d_lr rsq,e,de,f,df;
  };

  struct TableHost {
    HAT::t_int_2d_lr tabindex;
    HAT::t_double_1d innersq,invdelta;
    HAT::t_double_2d_lr rsq,e,de,f,df;
  };

  TableDeviceConst d_table_const;
  TableDevice* d_table;
  TableHost* h_table;

  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  void allocate() override;
  int update_table;
  void create_kokkos_tables();

  KOKKOS_INLINE_FUNCTION
  void getMixingWeights(int, KK_FLOAT &, KK_FLOAT &, KK_FLOAT &, KK_FLOAT &) const;

  typename AT::t_kkfloat_1d d_mixWtSite1old,d_mixWtSite2old,d_mixWtSite1,d_mixWtSite2;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_1d rho;
  typename HAT::t_double_1d h_rho;
  typename AT::t_kkfloat_1d uCG, uCGnew;
  typename AT::t_kkfloat_2d dvector;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  DAT::tdual_int_scalar k_error_flag;

  DAT::ttransform_kkfloat_2d k_cutsq;
  typename AT::t_kkfloat_2d d_cutsq;

  int first;
  typename AT::t_int_1d d_sendlist;
  typename AT::t_double_1d_um v_buf;

  friend void pair_virial_fdotr_compute<PairMultiLucyRXKokkos>(PairMultiLucyRXKokkos*);
};

}

#endif
#endif

