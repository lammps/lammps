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

PairStyle(tip4p/long/kk,        PairCoulLongKokkos<LMPDeviceType,PairTIP4PLong,true,false>);
PairStyle(tip4p/long/kk/device, PairCoulLongKokkos<LMPDeviceType,PairTIP4PLong,true,false>);
PairStyle(tip4p/long/kk/host,   PairCoulLongKokkos<LMPHostType,  PairTIP4PLong,true,false>);

PairStyle(coul/long/soft/kk,        PairCoulLongKokkos<LMPDeviceType,PairCoulLongSoft,false,true>);
PairStyle(coul/long/soft/kk/device, PairCoulLongKokkos<LMPDeviceType,PairCoulLongSoft,false,true>);
PairStyle(coul/long/soft/kk/host,   PairCoulLongKokkos<LMPHostType,  PairCoulLongSoft,false,true>);

PairStyle(tip4p/long/soft/kk,        PairCoulLongKokkos<LMPDeviceType,PairTIP4PLongSoft,true,true>);
PairStyle(tip4p/long/soft/kk/device, PairCoulLongKokkos<LMPDeviceType,PairTIP4PLongSoft,true,true>);
PairStyle(tip4p/long/soft/kk/host,   PairCoulLongKokkos<LMPHostType,  PairTIP4PLongSoft,true,true>);

// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_COUL_LONG_KOKKOS_H
#define LMP_PAIR_COUL_LONG_KOKKOS_H

#include "pair_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType, class PairCoulLongBase, bool TIP4P, bool SOFT>
class PairCoulLongKokkos : public PairKokkos<DeviceType,PairCoulLongBase,false,TIP4P,SOFT>
{
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG = TIP4P ? 2 : 1};  // 2=COUL_TIP4P, 1=COUL_LONG
  PairCoulLongKokkos(class LAMMPS *);

};

}

#endif // !LMP_PAIR_COUL_LONG_KOKKOS_H
#endif

