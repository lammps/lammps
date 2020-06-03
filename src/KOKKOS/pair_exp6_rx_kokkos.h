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

#ifdef PAIR_CLASS

PairStyle(exp6/rx/kk,PairExp6rxKokkos<Device>)
PairStyle(exp6/rx/kk/device,PairExp6rxKokkos<Device>)
PairStyle(exp6/rx/kk/host,PairExp6rxKokkos<Host>)

#else

#ifndef LMP_PAIR_EXP6_RX_KOKKOS_H
#define LMP_PAIR_EXP6_RX_KOKKOS_H

#include "pair_exp6_rx.h"
#include "kokkos_type.h"
#include "pair_kokkos.h"

namespace LAMMPS_NS {

// Create a structure to hold the parameter data for all
// local and neighbor particles. Pack inside this struct
// to avoid any name clashes.

template<ExecutionSpace Space>
struct PairExp6ParamDataTypeKokkos
{
  typedef ArrayTypes<Space> AT;

   int n;
   typename AT::t_float_1d epsilon1, alpha1, rm1, mixWtSite1,
          epsilon2, alpha2, rm2, mixWtSite2,
          epsilonOld1, alphaOld1, rmOld1, mixWtSite1old,
          epsilonOld2, alphaOld2, rmOld2, mixWtSite2old;

   // Default constructor -- nullify everything.
   PairExp6ParamDataTypeKokkos<Space>(void)
      : n(0)
   {}
};

template<ExecutionSpace Space>
struct PairExp6ParamDataTypeKokkosVect
{
  typedef ArrayTypes<Space> AT;

   typename AT::t_float_1d epsilon, rm3, alpha, xMolei, epsilon_old, rm3_old,
                           alpha_old, xMolei_old, fractionOFA, fraction1,
                           fraction2, nMoleculesOFA, nMolecules1, nMolecules2,
                           nTotal, fractionOFAold, fractionOld1, fractionOld2,
                           nMoleculesOFAold, nMoleculesOld1, nMoleculesOld2,
                           nTotalold;

   // Default constructor -- nullify everything.
   PairExp6ParamDataTypeKokkosVect<Space>(void)
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

template<ExecutionSpace Space>
class PairExp6rxKokkos : public PairExp6rx {
 public:
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;
  typedef typename GetFloatType<Space>::type SPACE_FLOAT;
  typedef EV_FLOAT value_type;

  PairExp6rxKokkos(class LAMMPS *);
  virtual ~PairExp6rxKokkos();
  void compute(int, int);
  void coeff(int, char **);
  void init_style();

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

  typename AT::t_float_1d_3_randomread x;
  typename AT::t_float_1d_3 f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_float_1d uCG, uCGnew;
  typename AT::t_float_2d dvector;

  typedef Kokkos::View<typename AT::t_float_1d::data_type*[3],Kokkos::LayoutRight,DeviceType> t_float_1d_3_thread;
  typedef Kokkos::View<typename AT::t_float_1d::data_type*,Kokkos::LayoutRight,DeviceType> t_float_1d_thread;

  t_float_1d_3_thread t_f;
  t_float_1d_thread t_uCG, t_uCGnew;

  DAT::tdual_float_1d k_eatom;
  DAT::tdual_float_1d_6 k_vatom;
  typename AT::t_float_1d d_eatom;
  typename AT::t_float_1d_6 d_vatom;

  DAT::tdual_int_scalar k_error_flag;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  PairExp6ParamDataTypeKokkos<Space> PairExp6ParamData;
  PairExp6ParamDataTypeKokkosVect<Space> PairExp6ParamDataVect;

  void allocate();
  DAT::tdual_int_1d k_mol2param;               // mapping from molecule to parameters
  typename AT::t_int_1d_randomread d_mol2param;

  typedef Kokkos::DualView<Param*,Kokkos::LayoutRight,DeviceType> tdual_param_1d;
  typedef typename tdual_param_1d::t_dev_const_randomread t_param_1d_randomread;

  tdual_param_1d k_params;                // parameter set for an I-J-K interaction
  t_param_1d_randomread d_params;                // parameter set for an I-J-K interaction

  DAT::tdual_float_2d k_cutsq;
  typename AT::t_float_2d d_cutsq;

  void read_file(char *);
  void setup();

  KOKKOS_INLINE_FUNCTION
  void getMixingWeights(int, SPACE_FLOAT &, SPACE_FLOAT &, SPACE_FLOAT &, SPACE_FLOAT &, SPACE_FLOAT &, SPACE_FLOAT &, SPACE_FLOAT &, SPACE_FLOAT &, SPACE_FLOAT &, SPACE_FLOAT &, SPACE_FLOAT &, SPACE_FLOAT &, SPACE_FLOAT &, SPACE_FLOAT &, SPACE_FLOAT &, SPACE_FLOAT &) const;

  template <class ArrayT>
  void getMixingWeightsVect(const int, int, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &, ArrayT &) const;

  KOKKOS_INLINE_FUNCTION
  void exponentScaling(KK_FLOAT, SPACE_FLOAT &, SPACE_FLOAT &) const;

  KOKKOS_INLINE_FUNCTION
  void polynomialScaling(KK_FLOAT, SPACE_FLOAT &, SPACE_FLOAT &, SPACE_FLOAT &) const;

  KK_FLOAT s_coeffAlpha[6],s_coeffEps[6],s_coeffRm[6];

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT func_rin(const SPACE_FLOAT &) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT expValue(const KK_FLOAT) const;

  friend void pair_virial_fdotr_compute<Space,PairExp6rxKokkos>(PairExp6rxKokkos*);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E:  alpha_ij is 6.0 in pair exp6

Self-explanatory

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: PairExp6rxKokkos requires a fix rx command

The fix rx command must come before the pair style command in the input file

E:  There are no rx species specified

There must be at least one species specified through the fix rx command

E:  Site1 name not recognized in pair coefficients

The site1 keyword does not match the species keywords specified throug the fix rx command

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E:  Cannot open exp6/rx potential file %s

Self-explanatory

E:  Incorrect format in exp6/rx potential file

Self-explanatory

E:  Illegal exp6/rx parameters.  Rm and Epsilon must be greater than zero.  Alpha cannot be negative.

Self-explanatory

E:  Illegal exp6/rx parameters.  Interaction potential does not exist.

Self-explanatory

E:  Potential file has duplicate entry.

Self-explanatory

E:  The number of molecules in CG particle is less than 10*DBL_EPSILON.

Self-explanatory.  Check the species concentrations have been properly set
and check the reaction kinetic solver parameters in fix rx to more for
sufficient accuracy.


*/
