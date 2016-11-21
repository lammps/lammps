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

PairStyle(dpd/fdt/energy/kk,PairDPDfdtEnergyKokkos<LMPDeviceType>)
PairStyle(dpd/fdt/energy/kk/device,PairDPDfdtEnergyKokkos<LMPDeviceType>)
PairStyle(dpd/fdt/energy/kk/host,PairDPDfdtEnergyKokkos<LMPHostType>)

#else

#ifndef LMP_PAIR_DPD_FDT_ENERGY_KOKKOS_H
#define LMP_PAIR_DPD_FDT_ENERGY_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_dpd_fdt_energy.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairDPDfdtEnergyKokkos : public PairDPDfdtEnergy {
 public:
  enum {EnabledNeighFlags=HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  PairDPDfdtEnergyKokkos(class LAMMPS *);
  virtual ~PairDPDfdtEnergyKokkos();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void init_style();
  double init_one(int, int);

 protected:
  void cleanup_copy();

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_fpair(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_evdwl(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

  F_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename ArrayTypes<DeviceType>::t_x_array_randomread x;
  typename ArrayTypes<DeviceType>::t_x_array c_x;
  typename ArrayTypes<DeviceType>::t_f_array f;
  typename ArrayTypes<DeviceType>::t_int_1d_randomread type;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename ArrayTypes<DeviceType>::t_efloat_1d d_eatom;
  typename ArrayTypes<DeviceType>::t_virial_array d_vatom;
  typename ArrayTypes<DeviceType>::t_tagint_1d tag;

  int newton_pair;
  double special_lj[4];

  typename ArrayTypes<DeviceType>::tdual_ffloat_2d k_cutsq;
  typename ArrayTypes<DeviceType>::t_ffloat_2d d_cutsq;


  int neighflag;
  int nlocal,nall,eflag,vflag;

  void allocate();

  friend class PairComputeFunctor<PairDPDfdtEnergyKokkos,HALF,true>;
  friend class PairComputeFunctor<PairDPDfdtEnergyKokkos,HALFTHREAD,true>;
  friend class PairComputeFunctor<PairDPDfdtEnergyKokkos,HALF,false>;
  friend class PairComputeFunctor<PairDPDfdtEnergyKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairDPDfdtEnergyKokkos,HALF,void>(PairDPDfdtEnergyKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairDPDfdtEnergyKokkos,HALFTHREAD,void>(PairDPDfdtEnergyKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairDPDfdtEnergyKokkos,void>(PairDPDfdtEnergyKokkos*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairDPDfdtEnergyKokkos>(PairDPDfdtEnergyKokkos*);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair dpd/fdt/energy requires ghost atoms store velocity

Use the communicate vel yes command to enable this.

E: Pair dpd/fdt/energy requires newton pair on

Self-explanatory.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

*/
