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
PairStyle(exp6/rx/kk,PairExp6rxKokkos<LMPDeviceType>);
PairStyle(exp6/rx/kk/device,PairExp6rxKokkos<LMPDeviceType>);
PairStyle(exp6/rx/kk/host,PairExp6rxKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_EXP6_RX_KOKKOS_H
#define LMP_PAIR_EXP6_RX_KOKKOS_H

#include "pair_exp6_rx.h"
#include "kokkos_type.h"
#include "pair_kokkos.h"

namespace LAMMPS_NS {

// Create a structure to hold the parameter data for all
// local and neighbor particles. Pack inside this struct
// to avoid any name clashes.

template<class DeviceType>
struct PairExp6ParamDataTypeKokkos
{
  typedef ArrayTypes<DeviceType> AT;

   int n;
   typename AT::t_kkfloat_1d epsilon1, alpha1, rm1, mixWtSite1,
          epsilon2, alpha2, rm2, mixWtSite2,
          epsilonOld1, alphaOld1, rmOld1, mixWtSite1old,
          epsilonOld2, alphaOld2, rmOld2, mixWtSite2old;

   // Default constructor -- nullify everything.
   PairExp6ParamDataTypeKokkos()
      : n(0)
   {}
};

template<class DeviceType>
struct PairExp6ParamDataTypeKokkosVect
{
  typedef ArrayTypes<DeviceType> AT;

   typename AT::t_kkfloat_1d epsilon, rm3, alpha, xMolei, epsilon_old, rm3_old,
                           alpha_old, xMolei_old, fractionOFA, fraction1,
                           fraction2, nMoleculesOFA, nMolecules1, nMolecules2,
                           nTotal, fractionOFAold, fractionOld1, fractionOld2,
                           nMoleculesOFAold, nMoleculesOld1, nMoleculesOld2,
                           nTotalold;

   // Default constructor -- nullify everything.
   PairExp6ParamDataTypeKokkosVect()
   {}
};

struct TagPairExp6rxZeroMixingWeights{};
struct TagPairExp6rxgetMixingWeights{};

template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairExp6rxCompute{};

template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairExp6rxComputeNoAtomics{};

struct TagPairExp6rxCollapseDupViews{};
struct TagPairExp6rxZeroDupViews{};

template<class DeviceType>
class PairExp6rxKokkos : public PairExp6rx {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairExp6rxKokkos(class LAMMPS *);
  ~PairExp6rxKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairExp6rxZeroMixingWeights, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairExp6rxgetMixingWeights, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairExp6rxCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairExp6rxCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairExp6rxComputeNoAtomics<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool Site1EqSite2, bool UseAtomics, bool OneType>
  KOKKOS_INLINE_FUNCTION
  void vectorized_operator(const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairExp6rxComputeNoAtomics<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairExp6rxCollapseDupViews, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairExp6rxZeroDupViews, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                  const KK_FLOAT &dely, const KK_FLOAT &delz) const;

  KOKKOS_INLINE_FUNCTION
  int sbmask(const int& j) const;

 protected:
  int eflag,vflag;
  int nlocal,newton_pair,neighflag;
  KK_FLOAT special_lj[4];
  int nthreads,ntypes;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_1d uCG, uCGnew;
  typename AT::t_kkfloat_2d dvector;

  typedef Kokkos::View<KK_FLOAT**[3],Kokkos::LayoutRight,DeviceType> t_kkfloat_1d_3_thread;
  typedef Kokkos::View<KK_FLOAT**,Kokkos::LayoutRight,DeviceType> t_kkfloat_1d_thread;

  t_kkfloat_1d_3_thread t_f;
  t_kkfloat_1d_thread t_uCG, t_uCGnew;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  DAT::tdual_int_scalar k_error_flag;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  PairExp6ParamDataTypeKokkos<DeviceType> PairExp6ParamData;
  PairExp6ParamDataTypeKokkosVect<DeviceType> PairExp6ParamDataVect;

  void allocate() override;
  DAT::tdual_int_1d k_mol2param;               // mapping from molecule to parameters
  typename AT::t_int_1d_randomread d_mol2param;

  typedef Kokkos::DualView<Param*,Kokkos::LayoutRight,DeviceType> tdual_param_1d;
  typedef typename tdual_param_1d::t_dev_const_randomread t_param_1d_randomread;

  tdual_param_1d k_params;                // parameter set for an I-J-K interaction
  t_param_1d_randomread d_params;                // parameter set for an I-J-K interaction

  DAT::ttransform_kkfloat_2d k_cutsq;
  typename AT::t_kkfloat_2d d_cutsq;

  void read_file(char *) override;
  void setup() override;

  KOKKOS_INLINE_FUNCTION
  void getMixingWeights(int, KK_FLOAT &, KK_FLOAT &, KK_FLOAT &, KK_FLOAT &, KK_FLOAT &, KK_FLOAT &, KK_FLOAT &, KK_FLOAT &, KK_FLOAT &, KK_FLOAT &, KK_FLOAT &, KK_FLOAT &, KK_FLOAT &, KK_FLOAT &, KK_FLOAT &, KK_FLOAT &) const;

  template <class ArrayT>
  void getMixingWeightsVect(const int, int, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &) const;

  KOKKOS_INLINE_FUNCTION
  void exponentScaling(KK_FLOAT, KK_FLOAT &, KK_FLOAT &) const;

  KOKKOS_INLINE_FUNCTION
  void polynomialScaling(KK_FLOAT, KK_FLOAT &, KK_FLOAT &, KK_FLOAT &) const;

  double s_coeffAlpha[6],s_coeffEps[6],s_coeffRm[6];

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT func_rin(const KK_FLOAT &) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT expValue(const KK_FLOAT) const;

  friend void pair_virial_fdotr_compute<PairExp6rxKokkos>(PairExp6rxKokkos*);
};

}

#endif
#endif

