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

PairStyle(lj/cut/coul/dsf/kk,PairLJCutCoulDSFKokkos<Device>)
PairStyle(lj/cut/coul/dsf/kk/device,PairLJCutCoulDSFKokkos<Device>)
PairStyle(lj/cut/coul/dsf/kk/host,PairLJCutCoulDSFKokkos<Host>)

#else

#ifndef LMP_PAIR_LJ_CUT_COUL_DSF_KOKKOS_H
#define LMP_PAIR_LJ_CUT_COUL_DSF_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_cut_coul_dsf.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<ExecutionSpace Space>
class PairLJCutCoulDSFKokkos : public PairLJCutCoulDSF {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=1};
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;
  PairLJCutCoulDSFKokkos(class LAMMPS *);
  ~PairLJCutCoulDSFKokkos();

  void compute(int, int);

  void init_style();
  double init_one(int, int);

 protected:
  void cleanup_copy();

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fpair(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fcoul(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_evdwl(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_ecoul(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const;

  Kokkos::DualView<params_lj_coul**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_lj_coul**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  // hardwired to space for 12 atom types
  params_lj_coul m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cut_ljsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cut_coulsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_float_1d_3_randomread x;
  typename AT::t_float_1d_3 c_x;
  typename AT::t_float_1d_3 f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_float_1d_randomread q;

  DAT::tdual_float_1d k_eatom;
  DAT::tdual_float_1d_6 k_vatom;
  typename AT::t_float_1d d_eatom;
  typename AT::t_float_1d_6 d_vatom;

  int newton_pair;

  DAT::tdual_float_2d k_cutsq;
  typename AT::t_float_2d d_cutsq;
  DAT::tdual_float_2d k_cut_ljsq;
  typename AT::t_float_2d d_cut_ljsq;
  DAT::tdual_float_2d k_cut_coulsq;
  typename AT::t_float_2d d_cut_coulsq;


  int neighflag;
  int nlocal,nall,eflag,vflag;

  KK_FLOAT special_coul[4];
  KK_FLOAT special_lj[4];
  KK_FLOAT qqrd2e;

  void allocate();
  friend class PairComputeFunctor<Space,PairLJCutCoulDSFKokkos,FULL,true>;
  friend class PairComputeFunctor<Space,PairLJCutCoulDSFKokkos,HALF,true>;
  friend class PairComputeFunctor<Space,PairLJCutCoulDSFKokkos,HALFTHREAD,true>;
  friend class PairComputeFunctor<Space,PairLJCutCoulDSFKokkos,FULL,false>;
  friend class PairComputeFunctor<Space,PairLJCutCoulDSFKokkos,HALF,false>;
  friend class PairComputeFunctor<Space,PairLJCutCoulDSFKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<Space,PairLJCutCoulDSFKokkos,FULL,void>(PairLJCutCoulDSFKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute_neighlist<Space,PairLJCutCoulDSFKokkos,HALF,void>(PairLJCutCoulDSFKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute_neighlist<Space,PairLJCutCoulDSFKokkos,HALFTHREAD,void>(PairLJCutCoulDSFKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute<Space,PairLJCutCoulDSFKokkos,void>(PairLJCutCoulDSFKokkos*,
                                                            NeighListKokkos<Space>*);
  friend void pair_virial_fdotr_compute<Space,PairLJCutCoulDSFKokkos>(PairLJCutCoulDSFKokkos*);

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot use Kokkos pair style with rRESPA inner/middle

Self-explanatory.

E: Cannot use chosen neighbor list style with lj/cut/coul/cut/kk

That style is not supported by Kokkos.

*/
