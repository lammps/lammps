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

PairStyle(table/rx/kk,PairTableKokkos<LMPDeviceType>)
PairStyle(table/rx/kk/device,PairTableKokkos<LMPDeviceType>)
PairStyle(table/rx/kk/host,PairTableKokkos<LMPHostType>)

#else

#ifndef LMP_PAIR_TABLE_RX_KOKKOS_H
#define LMP_PAIR_TABLE_RX_KOKKOS_H

#include "pair_table_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairTableRXKokkos : public PairTable {
 public:

  using DAT = ArrayTypes<DeviceType>;

  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF|N2};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;

  PairTableRXKokkos(class LAMMPS *);
  virtual ~PairTableRXKokkos();

  virtual void compute(int, int);

  template<int TABSTYLE>
  void compute_style(int, int);

  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  virtual double single(int, int, int, int, double, double, double, double &);

  void init_style();


 protected:

  struct TableDeviceConst {
    typename ArrayTypes<DeviceType>::t_ffloat_2d cutsq;
    typename ArrayTypes<DeviceType>::t_int_2d tabindex;
    typename ArrayTypes<DeviceType>::t_int_1d nshiftbits,nmask;
    typename ArrayTypes<DeviceType>::t_ffloat_1d innersq,invdelta,deltasq6;
    typename ArrayTypes<DeviceType>::t_ffloat_2d_randomread rsq,drsq,e,de,f,df,e2,f2;
  };

  struct TableDevice {
    typename ArrayTypes<DeviceType>::t_ffloat_2d cutsq;
    typename ArrayTypes<DeviceType>::t_int_2d tabindex;
    typename ArrayTypes<DeviceType>::t_int_1d nshiftbits,nmask;
    typename ArrayTypes<DeviceType>::t_ffloat_1d innersq,invdelta,deltasq6;
    typename ArrayTypes<DeviceType>::t_ffloat_2d rsq,drsq,e,de,f,df,e2,f2;
  };

  struct TableHost {
    typename ArrayTypes<LMPHostType>::t_ffloat_2d cutsq;
    typename ArrayTypes<LMPHostType>::t_int_2d tabindex;
    typename ArrayTypes<LMPHostType>::t_int_1d nshiftbits,nmask;
    typename ArrayTypes<LMPHostType>::t_ffloat_1d innersq,invdelta,deltasq6;
    typename ArrayTypes<LMPHostType>::t_ffloat_2d rsq,drsq,e,de,f,df,e2,f2;
  };

  TableDeviceConst d_table_const;
  TableDevice* d_table;
  TableHost* h_table;

  F_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  typename ArrayTypes<DeviceType>::t_ffloat_2d d_cutsq;

  virtual void allocate();
  void compute_table(Table *);

  typename ArrayTypes<DeviceType>::t_x_array_randomread x;
  typename ArrayTypes<DeviceType>::t_x_array_const c_x;
  typename ArrayTypes<DeviceType>::t_f_array f;
  typename ArrayTypes<DeviceType>::t_int_1d_randomread type;
  typename ArrayTypes<DeviceType>::t_efloat_1d d_eatom;
  typename ArrayTypes<DeviceType>::t_virial_array d_vatom;

 protected:
  int nlocal,nall,eflag,vflag,neighflag,newton_pair;

  int update_table;
  void create_kokkos_tables();
  void cleanup_copy();

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_fpair(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_evdwl(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_ecoul(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const {
    return 0;
  }

  friend class PairComputeFunctor<PairTableRXKokkos,FULL,true,S_TableCompute<DeviceType,LOOKUP> >;
  friend class PairComputeFunctor<PairTableRXKokkos,HALF,true,S_TableCompute<DeviceType,LOOKUP> >;
  friend class PairComputeFunctor<PairTableRXKokkos,HALFTHREAD,true,S_TableCompute<DeviceType,LOOKUP> >;
  friend class PairComputeFunctor<PairTableRXKokkos,N2,true,S_TableCompute<DeviceType,LOOKUP> >;
  friend class PairComputeFunctor<PairTableRXKokkos,FULL,false,S_TableCompute<DeviceType,LOOKUP> >;
  friend class PairComputeFunctor<PairTableRXKokkos,HALF,false,S_TableCompute<DeviceType,LOOKUP> >;
  friend class PairComputeFunctor<PairTableRXKokkos,HALFTHREAD,false,S_TableCompute<DeviceType,LOOKUP> >;
  friend class PairComputeFunctor<PairTableRXKokkos,N2,false,S_TableCompute<DeviceType,LOOKUP> >;

  friend class PairComputeFunctor<PairTableRXKokkos,FULL,true,S_TableCompute<DeviceType,LINEAR> >;
  friend class PairComputeFunctor<PairTableRXKokkos,HALF,true,S_TableCompute<DeviceType,LINEAR> >;
  friend class PairComputeFunctor<PairTableRXKokkos,HALFTHREAD,true,S_TableCompute<DeviceType,LINEAR> >;
  friend class PairComputeFunctor<PairTableRXKokkos,N2,true,S_TableCompute<DeviceType,LINEAR> >;
  friend class PairComputeFunctor<PairTableRXKokkos,FULL,false,S_TableCompute<DeviceType,LINEAR> >;
  friend class PairComputeFunctor<PairTableRXKokkos,HALF,false,S_TableCompute<DeviceType,LINEAR> >;
  friend class PairComputeFunctor<PairTableRXKokkos,HALFTHREAD,false,S_TableCompute<DeviceType,LINEAR> >;
  friend class PairComputeFunctor<PairTableRXKokkos,N2,false,S_TableCompute<DeviceType,LINEAR> >;

  friend class PairComputeFunctor<PairTableRXKokkos,FULL,true,S_TableCompute<DeviceType,SPLINE> >;
  friend class PairComputeFunctor<PairTableRXKokkos,HALF,true,S_TableCompute<DeviceType,SPLINE> >;
  friend class PairComputeFunctor<PairTableRXKokkos,HALFTHREAD,true,S_TableCompute<DeviceType,SPLINE> >;
  friend class PairComputeFunctor<PairTableRXKokkos,N2,true,S_TableCompute<DeviceType,SPLINE> >;
  friend class PairComputeFunctor<PairTableRXKokkos,FULL,false,S_TableCompute<DeviceType,SPLINE> >;
  friend class PairComputeFunctor<PairTableRXKokkos,HALF,false,S_TableCompute<DeviceType,SPLINE> >;
  friend class PairComputeFunctor<PairTableRXKokkos,HALFTHREAD,false,S_TableCompute<DeviceType,SPLINE> >;
  friend class PairComputeFunctor<PairTableRXKokkos,N2,false,S_TableCompute<DeviceType,SPLINE> >;

  friend class PairComputeFunctor<PairTableRXKokkos,FULL,true,S_TableCompute<DeviceType,BITMAP> >;
  friend class PairComputeFunctor<PairTableRXKokkos,HALF,true,S_TableCompute<DeviceType,BITMAP> >;
  friend class PairComputeFunctor<PairTableRXKokkos,HALFTHREAD,true,S_TableCompute<DeviceType,BITMAP> >;
  friend class PairComputeFunctor<PairTableRXKokkos,N2,true,S_TableCompute<DeviceType,BITMAP> >;
  friend class PairComputeFunctor<PairTableRXKokkos,FULL,false,S_TableCompute<DeviceType,BITMAP> >;
  friend class PairComputeFunctor<PairTableRXKokkos,HALF,false,S_TableCompute<DeviceType,BITMAP> >;
  friend class PairComputeFunctor<PairTableRXKokkos,HALFTHREAD,false,S_TableCompute<DeviceType,BITMAP> >;
  friend class PairComputeFunctor<PairTableRXKokkos,N2,false,S_TableCompute<DeviceType,BITMAP> >;

  friend void pair_virial_fdotr_compute<PairTableRXKokkos>(PairTableRXKokkos*);

  /* PairTableRX members */

  int nspecies;
  char *site1, *site2;
  int isite1, isite2;
  bool fractionalWeighting;

  template <class ExecDevice>
  KOKKOS_INLINE_FUNCTION
  void getMixingWeights(typename ArrayTypes<ExecDevice>::t_float_2d_randomread,
      int, double &, double &, double &, double &);
};

}

#endif
#endif

/* ERROR/WARNING messages:

 */
