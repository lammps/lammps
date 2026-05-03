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

PairStyle(coul/long/kk,        PairCoulLongKokkos<LMPDeviceType,PairCoulLong,false,false>);
PairStyle(coul/long/kk/device, PairCoulLongKokkos<LMPDeviceType,PairCoulLong,false,false>);
PairStyle(coul/long/kk/host,   PairCoulLongKokkos<LMPHostType,  PairCoulLong,false,false>);

/*
PairStyle(tip4p/long/kk,        PairCoulLongKokkos<LMPDeviceType,PairTIP4PLong,true,false>);
PairStyle(tip4p/long/kk/device, PairCoulLongKokkos<LMPDeviceType,PairTIP4PLong,true,false>);
PairStyle(tip4p/long/kk/host,   PairCoulLongKokkos<LMPHostType,  PairTIP4PLong,true,false>);


PairStyle(coul/long/soft/kk,        PairCoulLongKokkos<LMPDeviceType,PairCoulLongSoft,false,true>);
PairStyle(coul/long/soft/kk/device, PairCoulLongKokkos<LMPDeviceType,PairCoulLongSoft,false,true>);
PairStyle(coul/long/soft/kk/host,   PairCoulLongKokkos<LMPHostType,  PairCoulLongSoft,false,true>);


PairStyle(tip4p/long/soft/kk,        PairCoulLongKokkos<LMPDeviceType,PairTIP4PLongSoft,true,true>);
PairStyle(tip4p/long/soft/kk/device, PairCoulLongKokkos<LMPDeviceType,PairTIP4PLongSoft,true,true>);
PairStyle(tip4p/long/soft/kk/host,   PairCoulLongKokkos<LMPHostType,  PairTIP4PLongSoft,true,true>);
*/

// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_COUL_LONG_KOKKOS_H
#define LMP_PAIR_COUL_LONG_KOKKOS_H

#include "pair_kokkos.h"
#include "neigh_list_kokkos.h"
#include "fix.h"

#include "pair_coul_long.h"
#include "pair_tip4p_long.h"
#ifdef LMP_FEP
  #include "pair_coul_long_soft.h"
  #include "pair_tip4p_long_soft.h"
#endif

namespace LAMMPS_NS {

template<class DeviceType, class PairCoulLongBase, bool TIP4P, bool SOFT>
class PairCoulLongKokkos : public PairKokkos<DeviceType,PairCoulLongBase,false,TIP4P,SOFT>
{
 public:

  PairCoulLongKokkos(class LAMMPS *);


 protected:


  

  friend struct PairComputeFunctor<PairCoulLongKokkos,FULL,true,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,FULL,true,1,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,HALF,true,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,HALFTHREAD,true,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,FULL,false,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,FULL,false,1,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,HALF,false,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,HALFTHREAD,false,0,CoulLongTable<1>>;
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongKokkos,FULL,0,CoulLongTable<1>>(PairCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongKokkos,FULL,1,CoulLongTable<1>>(PairCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongKokkos,HALF,0,CoulLongTable<1>>(PairCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongKokkos,HALFTHREAD,0,CoulLongTable<1>>(PairCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairCoulLongKokkos,CoulLongTable<1>>(PairCoulLongKokkos*,
                                                            NeighListKokkos<DeviceType>*);
  friend struct PairComputeFunctor<PairCoulLongKokkos,FULL,true,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,FULL,true,1,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,HALF,true,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,HALFTHREAD,true,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,FULL,false,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,FULL,false,1,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,HALF,false,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairCoulLongKokkos,HALFTHREAD,false,0,CoulLongTable<0>>;
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongKokkos,FULL,0,CoulLongTable<0>>(PairCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongKokkos,FULL,1,CoulLongTable<0>>(PairCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongKokkos,HALF,0,CoulLongTable<0>>(PairCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulLongKokkos,HALFTHREAD,0,CoulLongTable<0>>(PairCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairCoulLongKokkos,CoulLongTable<0>>(PairCoulLongKokkos*,
                                                            NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairCoulLongKokkos>(PairCoulLongKokkos*);

};

}

#endif // !LMP_PAIR_COUL_LONG_KOKKOS_H
#endif

