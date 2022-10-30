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
PairStyle(table/kk,PairTableKokkos<LMPDeviceType>);
PairStyle(table/kk/device,PairTableKokkos<LMPDeviceType>);
PairStyle(table/kk/host,PairTableKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_TABLE_KOKKOS_H
#define LMP_PAIR_TABLE_KOKKOS_H

#include "pair_table.h"
#include "pair_kokkos.h"
#include "neigh_list_kokkos.h"
#include "atom_kokkos.h"

namespace LAMMPS_NS {

template<class Device,int TABSTYLE>
struct S_TableCompute {
  static constexpr int TabStyle = TABSTYLE;
};

template <class DeviceType, int NEIGHFLAG, int TABSTYLE>
struct PairTableComputeFunctor;

template<class DeviceType>
class PairTableKokkos : public PairTable {
 public:

  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  PairTableKokkos(class LAMMPS *);
  ~PairTableKokkos() override;

  void compute(int, int) override;

  template<int TABSTYLE>
  void compute_style(int, int);

  void settings(int, char **) override;
  double init_one(int, int) override;

  void init_style() override;


 protected:

  /*struct TableDeviceConst {
    typename AT::t_ffloat_2d_randomread cutsq;
    typename AT::t_int_2d_randomread tabindex;
    typename AT::t_int_1d_randomread nshiftbits,nmask;
    typename AT::t_ffloat_1d_randomread innersq,invdelta,deltasq6;
    typename AT::t_ffloat_2d_randomread rsq,drsq,e,de,f,df,e2,f2;
  };*/
 //Its faster not to use texture fetch if the number of tables is less than 32!
  struct TableDeviceConst {
    typename AT::t_ffloat_2d cutsq;
    typename AT::t_int_2d tabindex;
    typename AT::t_int_1d nshiftbits,nmask;
    typename AT::t_ffloat_1d innersq,invdelta,deltasq6;
    typename AT::t_ffloat_2d_randomread rsq,drsq,e,de,f,df,e2,f2;
  };

  struct TableDevice {
    typename AT::t_ffloat_2d cutsq;
    typename AT::t_int_2d tabindex;
    typename AT::t_int_1d nshiftbits,nmask;
    typename AT::t_ffloat_1d innersq,invdelta,deltasq6;
    typename AT::t_ffloat_2d rsq,drsq,e,de,f,df,e2,f2;
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

  typename AT::t_ffloat_2d d_cutsq;

  void allocate() override;
  void compute_table(Table *) override;

  typename AT::t_x_array_randomread x;
  typename AT::t_x_array_const c_x;
  typename AT::t_f_array f;
  typename AT::t_int_1d_randomread type;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename AT::t_efloat_1d d_eatom;
  typename AT::t_virial_array d_vatom;

 protected:
  int nlocal,nall,eflag,vflag,neighflag,newton_pair;

  int update_table;
  void create_kokkos_tables();

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_fpair(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_evdwl(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_ecoul(const F_FLOAT& /*rsq*/, const int& /*i*/, const int& /*j*/,
                        const int& /*itype*/, const int& /*jtype*/) const { return 0; }

  friend struct PairComputeFunctor<PairTableKokkos,FULL,true,S_TableCompute<DeviceType,LOOKUP> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALF,true,S_TableCompute<DeviceType,LOOKUP> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALFTHREAD,true,S_TableCompute<DeviceType,LOOKUP> >;
  friend struct PairComputeFunctor<PairTableKokkos,FULL,false,S_TableCompute<DeviceType,LOOKUP> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALF,false,S_TableCompute<DeviceType,LOOKUP> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALFTHREAD,false,S_TableCompute<DeviceType,LOOKUP> >;

  friend struct PairComputeFunctor<PairTableKokkos,FULL,true,S_TableCompute<DeviceType,LINEAR> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALF,true,S_TableCompute<DeviceType,LINEAR> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALFTHREAD,true,S_TableCompute<DeviceType,LINEAR> >;
  friend struct PairComputeFunctor<PairTableKokkos,FULL,false,S_TableCompute<DeviceType,LINEAR> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALF,false,S_TableCompute<DeviceType,LINEAR> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALFTHREAD,false,S_TableCompute<DeviceType,LINEAR> >;

  friend struct PairComputeFunctor<PairTableKokkos,FULL,true,S_TableCompute<DeviceType,SPLINE> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALF,true,S_TableCompute<DeviceType,SPLINE> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALFTHREAD,true,S_TableCompute<DeviceType,SPLINE> >;
  friend struct PairComputeFunctor<PairTableKokkos,FULL,false,S_TableCompute<DeviceType,SPLINE> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALF,false,S_TableCompute<DeviceType,SPLINE> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALFTHREAD,false,S_TableCompute<DeviceType,SPLINE> >;

  friend struct PairComputeFunctor<PairTableKokkos,FULL,true,S_TableCompute<DeviceType,BITMAP> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALF,true,S_TableCompute<DeviceType,BITMAP> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALFTHREAD,true,S_TableCompute<DeviceType,BITMAP> >;
  friend struct PairComputeFunctor<PairTableKokkos,FULL,false,S_TableCompute<DeviceType,BITMAP> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALF,false,S_TableCompute<DeviceType,BITMAP> >;
  friend struct PairComputeFunctor<PairTableKokkos,HALFTHREAD,false,S_TableCompute<DeviceType,BITMAP> >;

  friend void pair_virial_fdotr_compute<PairTableKokkos>(PairTableKokkos*);
};






}

#endif
#endif

