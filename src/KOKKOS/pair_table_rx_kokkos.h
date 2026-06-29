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
PairStyle(table/rx/kk,PairTableRXKokkos<LMPDeviceType>);
PairStyle(table/rx/kk/device,PairTableRXKokkos<LMPDeviceType>);
PairStyle(table/rx/kk/host,PairTableRXKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_TABLE_RX_KOKKOS_H
#define LMP_PAIR_TABLE_RX_KOKKOS_H

#include "pair_table_kokkos.h"
#include "kokkos_few.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairTableRXKokkos : public PairTable {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  PairTableRXKokkos(class LAMMPS *);
  ~PairTableRXKokkos() override;

  void compute(int, int) override;

  template<int TABSTYLE>
  void compute_style(int, int);

  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  double single(int, int, int, int, double, double, double, double &) override;

  void init_style() override;

  struct TableDeviceConst {
    typename AT::t_double_2d_lr cutsq;
    typename AT::t_int_2d_lr tabindex;
    typename AT::t_int_1d nshiftbits,nmask;
    typename AT::t_double_1d innersq,invdelta,deltasq6;
    typename AT::t_double_2d_lr_randomread rsq,drsq,e,de,f,df,e2,f2;
  };

  struct TableDevice {
    typename AT::t_double_2d_lr cutsq;
    typename AT::t_int_2d_lr tabindex;
    typename AT::t_int_1d nshiftbits,nmask;
    typename AT::t_double_1d innersq,invdelta,deltasq6;
    typename AT::t_double_2d_lr rsq,drsq,e,de,f,df,e2,f2;
  };

  struct TableHost {
    typename ArrayTypes<LMPHostType>::t_double_2d_lr cutsq;
    typename ArrayTypes<LMPHostType>::t_int_2d_lr tabindex;
    typename ArrayTypes<LMPHostType>::t_int_1d nshiftbits,nmask;
    typename ArrayTypes<LMPHostType>::t_double_1d innersq,invdelta,deltasq6;
    typename ArrayTypes<LMPHostType>::t_double_2d_lr rsq,drsq,e,de,f,df,e2,f2;
  };

  TableDeviceConst d_table_const;
  TableDevice* d_table;
  TableHost* h_table;

  Few<Few<double, MAX_TYPES_STACKPARAMS+1>, MAX_TYPES_STACKPARAMS+1> m_cutsq;

  typename AT::t_double_2d_lr d_cutsq;

  void allocate() override;
  void compute_table(Table *) override;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;

  int neighflag;

  int update_table;
  void create_kokkos_tables();

  friend void pair_virial_fdotr_compute<PairTableRXKokkos>(PairTableRXKokkos*);

  /* PairTableRX members */

  Kokkos::View<KK_FLOAT*, DeviceType> mixWtSite1old;
  Kokkos::View<KK_FLOAT*, DeviceType> mixWtSite2old;
  Kokkos::View<KK_FLOAT*, DeviceType> mixWtSite1;
  Kokkos::View<KK_FLOAT*, DeviceType> mixWtSite2;

  int nspecies;
  char *site1, *site2;
  int isite1, isite2;
  bool fractionalWeighting;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;
};

}

#endif
#endif

