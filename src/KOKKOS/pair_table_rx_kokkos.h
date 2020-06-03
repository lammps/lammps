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

PairStyle(table/rx/kk,PairTableRXKokkos<Device>)
PairStyle(table/rx/kk/device,PairTableRXKokkos<Device>)
PairStyle(table/rx/kk/host,PairTableRXKokkos<Host>)

#else

#ifndef LMP_PAIR_TABLE_RX_KOKKOS_H
#define LMP_PAIR_TABLE_RX_KOKKOS_H

#include "pair_table_kokkos.h"
#include "kokkos_few.h"

namespace LAMMPS_NS {

template<ExecutionSpace Space>
class PairTableRXKokkos : public PairTable {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;
  typedef typename GetFloatType<Space>::type SPACE_FLOAT;

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

  struct TableDual {
    DAT::tdual_float_2d k_cutsq;
    DAT::tdual_int_2d k_tabindex;
    DAT::tdual_int_1d k_nshiftbits,k_nmask;
    DAT::tdual_float_1d k_innersq,k_invdelta,k_deltasq6;
    DAT::tdual_float_2d k_rsq,k_drsq,k_e,k_de,k_f,k_df,k_e2,k_f2;
  };

  struct TableDeviceConst {
    typename AT::t_float_2d cutsq;
    typename AT::t_int_2d tabindex;
    typename AT::t_int_1d nshiftbits,nmask;
    typename AT::t_float_1d innersq,invdelta,deltasq6;
    typename AT::t_float_2d_randomread rsq,drsq,e,de,f,df,e2,f2;
  };

  struct TableDevice {
    typename AT::t_float_2d cutsq;
    typename AT::t_int_2d tabindex;
    typename AT::t_int_1d nshiftbits,nmask;
    typename AT::t_float_1d innersq,invdelta,deltasq6;
    typename AT::t_float_2d rsq,drsq,e,de,f,df,e2,f2;

    TableDevice(const TableDual *rhs) {
      cutsq = DualViewHelper<Space>::view(rhs->k_cutsq);
      tabindex = DualViewHelper<Space>::view(rhs->k_tabindex);
      nshiftbits = DualViewHelper<Space>::view(rhs->k_nshiftbits);
      nmask = DualViewHelper<Space>::view(rhs->k_nmask);
      innersq = DualViewHelper<Space>::view(rhs->k_innersq);
      invdelta = DualViewHelper<Space>::view(rhs->k_invdelta);
      deltasq6 = DualViewHelper<Space>::view(rhs->k_deltasq6);
      rsq = DualViewHelper<Space>::view(rhs->k_rsq);
      drsq = DualViewHelper<Space>::view(rhs->k_drsq);
      e = DualViewHelper<Space>::view(rhs->k_e);
      de = DualViewHelper<Space>::view(rhs->k_de);
      f = DualViewHelper<Space>::view(rhs->k_f);
      df = DualViewHelper<Space>::view(rhs->k_df);
      e2 = DualViewHelper<Space>::view(rhs->k_e2);
      f2 = DualViewHelper<Space>::view(rhs->k_f2);
    }
  };

  struct TableHost {
    HAT::t_float_2d cutsq;
    HAT::t_int_2d tabindex;
    HAT::t_int_1d nshiftbits,nmask;
    HAT::t_float_1d innersq,invdelta,deltasq6;
    HAT::t_float_2d rsq,drsq,e,de,f,df,e2,f2;

    TableHost(const TableDual *rhs) {
      cutsq = rhs->k_cutsq.h_view;
      tabindex = rhs->k_tabindex.h_view;
      nshiftbits = rhs->k_nshiftbits.h_view;
      nmask = rhs->k_nmask.h_view;
      innersq = rhs->k_innersq.h_view;
      invdelta = rhs->k_invdelta.h_view;
      deltasq6 = rhs->k_deltasq6.h_view;
      rsq = rhs->k_rsq.h_view;
      drsq = rhs->k_drsq.h_view;
      e = rhs->k_e.h_view;
      de = rhs->k_de.h_view;
      f = rhs->k_f.h_view;
      df = rhs->k_df.h_view;
      e2 = rhs->k_e2.h_view;
      f2 = rhs->k_f2.h_view;
    }
  };

  TableDual* k_table;
  TableDeviceConst d_table_const;
  TableDevice* d_table;
  TableHost* h_table;

  Few<Few<KK_FLOAT, MAX_TYPES_STACKPARAMS+1>, MAX_TYPES_STACKPARAMS+1> m_cutsq;

  typename AT::t_float_2d d_cutsq;

  virtual void allocate();
  void compute_table(Table *);

  typename AT::t_float_1d_3_randomread x;
  typename AT::t_float_1d_3 f;

  int neighflag;

  int update_table;
  void create_kokkos_tables();
  void cleanup_copy();

  friend void pair_virial_fdotr_compute<Space,PairTableRXKokkos>(PairTableRXKokkos*);

  /* PairTableRX members */

  typename AT::t_float_1d mixWtSite1old;
  typename AT::t_float_1d mixWtSite2old;
  typename AT::t_float_1d mixWtSite1;
  typename AT::t_float_1d mixWtSite2;

  int nspecies;
  char *site1, *site2;
  int isite1, isite2;
  bool fractionalWeighting;

  DAT::tdual_float_1d k_eatom;
  DAT::tdual_float_1d_6 k_vatom;
  typename AT::t_float_1d d_eatom;
  typename AT::t_float_1d_6 d_vatom;
};

}

#endif
#endif

/* ERROR/WARNING messages:

 */
