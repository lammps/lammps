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
#include "kokkos_few.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairTableRXKokkos : public PairTable {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
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

  Few<Few<F_FLOAT, MAX_TYPES_STACKPARAMS+1>, MAX_TYPES_STACKPARAMS+1> m_cutsq;

  typename ArrayTypes<DeviceType>::t_ffloat_2d d_cutsq;

  virtual void allocate();
  void compute_table(Table *);

  typename ArrayTypes<DeviceType>::t_x_array_randomread x;
  typename ArrayTypes<DeviceType>::t_f_array f;

  int neighflag;

  int update_table;
  void create_kokkos_tables();
  void cleanup_copy();

  friend void pair_virial_fdotr_compute<PairTableRXKokkos>(PairTableRXKokkos*);

  /* PairTableRX members */

  Kokkos::View<double*, DeviceType> mixWtSite1old;
  Kokkos::View<double*, DeviceType> mixWtSite2old;
  Kokkos::View<double*, DeviceType> mixWtSite1;
  Kokkos::View<double*, DeviceType> mixWtSite2;

  int nspecies;
  char *site1, *site2;
  int isite1, isite2;
  bool fractionalWeighting;

  typename ArrayTypes<DeviceType>::tdual_efloat_1d k_eatom;
  typename ArrayTypes<DeviceType>::tdual_virial_array k_vatom;
  typename ArrayTypes<DeviceType>::t_efloat_1d d_eatom;
  typename ArrayTypes<DeviceType>::t_virial_array d_vatom;
};

}

#endif
#endif

/* ERROR/WARNING messages:

 */
