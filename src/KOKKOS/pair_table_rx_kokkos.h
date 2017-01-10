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

PairStyle(table/rx/kk,PairTableRXKokkos<LMPDeviceType>)
PairStyle(table/rx/kk/device,PairTableRXKokkos<LMPDeviceType>)
PairStyle(table/rx/kk/host,PairTableRXKokkos<LMPHostType>)

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
  typename ArrayTypes<DeviceType>::t_efloat_1d uCG;
  typename ArrayTypes<DeviceType>::t_efloat_1d uCGnew;
  typename ArrayTypes<DeviceType>::t_efloat_1d d_eatom;
  typename ArrayTypes<DeviceType>::t_virial_array d_vatom;

  int nlocal,nall,eflag,vflag,neighflag,newton_pair;

  int update_table;
  void create_kokkos_tables();
  void cleanup_copy();

  template<bool STACKPARAMS, int TABSTYLE>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_fpair(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, int TABSTYLE>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_evdwl(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

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

  Kokkos::View<double*, DeviceType> mixWtSite1old_;
  Kokkos::View<double*, DeviceType> mixWtSite2old_;
  Kokkos::View<double*, DeviceType> mixWtSite1_;
  Kokkos::View<double*, DeviceType> mixWtSite2_;

  /* a duplicate of PairComputeFunctor to deal with uCG */
  template <int NEIGHFLAG, bool STACKPARAMS, int TABSTYLE>
  struct Functor {
    using device_type = DeviceType;
    typedef EV_FLOAT value_type;
    PairTableRXKokkos c;
    // arrays are atomic for Half(Thread) neighbor style
    Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,
                 device_type,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > f;
    Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,
                 device_type,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > uCG;
    Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,
                 device_type,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > uCGnew;
    NeighListKokkos<device_type> list;
    Functor(PairTableRXKokkos* c_ptr, NeighListKokkos<device_type>* list_ptr);
    ~Functor();
    KOKKOS_INLINE_FUNCTION int sbmask(const int& j) const {
      return j >> SBBITS & 3;
    }
    template<int EVFLAG, int NEWTON_PAIR>
    KOKKOS_INLINE_FUNCTION
    EV_FLOAT compute_item(const int&) const;
    KOKKOS_INLINE_FUNCTION
    void
    ev_tally(EV_FLOAT &ev, const int &i, const int &j,
             const F_FLOAT &epair, const F_FLOAT &fpair, const F_FLOAT &delx,
             const F_FLOAT &dely, const F_FLOAT &delz) const;
    KOKKOS_INLINE_FUNCTION
    void operator()(const int) const;
    KOKKOS_INLINE_FUNCTION
    void operator()(const int, value_type&) const;
  };

};

}

#endif
#endif

/* ERROR/WARNING messages:

 */
