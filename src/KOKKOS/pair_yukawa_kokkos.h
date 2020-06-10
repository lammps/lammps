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

PairStyle(yukawa/kk,PairYukawaKokkos<Device>)
PairStyle(yukawa/kk/device,PairYukawaKokkos<Device>)
PairStyle(yukawa/kk/host,PairYukawaKokkos<Host>)

#else

#ifndef LMP_PAIR_YUKAWA_KOKKOS_H
#define LMP_PAIR_YUKAWA_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_yukawa.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<ExecutionSpace Space>
class PairYukawaKokkos : public PairYukawa {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;

  PairYukawaKokkos(class LAMMPS *);
  virtual ~PairYukawaKokkos();

  void compute(int, int);
  void init_style();
  double init_one(int,int);

  struct params_yukawa {
    KOKKOS_INLINE_FUNCTION
    params_yukawa(){ cutsq=0, a = 0; offset = 0; }
    KOKKOS_INLINE_FUNCTION
    params_yukawa(int i){ cutsq=0, a = 0; offset = 0; }
    KK_FLOAT cutsq, a, offset;
  };


 protected:
  void cleanup_copy();

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fpair(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_evdwl(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_ecoul(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const
  {
    return 0;
  }


  Kokkos::DualView<params_yukawa**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_yukawa**,Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  params_yukawa m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_float_1d_3_lr_randomread x;
  typename AT::t_float_1d_3 f;
  typename AT::t_int_1d_randomread type;

  DAT::tdual_float_1d k_eatom;
  DAT::tdual_float_1d_6 k_vatom;
  typename AT::t_float_1d d_eatom;
  typename AT::t_float_1d_6 d_vatom;
  typename AT::t_tagint_1d tag;

  int newton_pair;
  KK_FLOAT special_lj[4];

  DAT::tdual_float_2d k_cutsq;
  typename AT::t_float_2d d_cutsq;


  int neighflag;
  int nlocal,nall,eflag,vflag;

  void allocate();
  friend class PairComputeFunctor<Space,PairYukawaKokkos,FULL,true>;
  friend class PairComputeFunctor<Space,PairYukawaKokkos,HALF,true>;
  friend class PairComputeFunctor<Space,PairYukawaKokkos,HALFTHREAD,true>;
  friend class PairComputeFunctor<Space,PairYukawaKokkos,FULL,false>;
  friend class PairComputeFunctor<Space,PairYukawaKokkos,HALF,false>;
  friend class PairComputeFunctor<Space,PairYukawaKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<Space,PairYukawaKokkos,FULL,void>(
    PairYukawaKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute_neighlist<Space,PairYukawaKokkos,HALF,void>(
    PairYukawaKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute_neighlist<Space,PairYukawaKokkos,HALFTHREAD,void>(
    PairYukawaKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute<Space,PairYukawaKokkos,void>(
    PairYukawaKokkos*,NeighListKokkos<Space>*);
  friend void pair_virial_fdotr_compute<Space,PairYukawaKokkos>(PairYukawaKokkos*);

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot use Kokkos pair style with rRESPA inner/middle

UNDOCUMENTED

E: Cannot use chosen neighbor list style with yukawa/kk

That style is not supported by Kokkos.

U: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

U: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

*/
