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

PairStyle(lj/cut/coul/long/kk,        PairLJCutCoulLongKokkos<LMPDeviceType,PairLJCutCoulLong,false,false>);
PairStyle(lj/cut/coul/long/kk/device, PairLJCutCoulLongKokkos<LMPDeviceType,PairLJCutCoulLong,false,false>);
PairStyle(lj/cut/coul/long/kk/host,   PairLJCutCoulLongKokkos<LMPHostType,  PairLJCutCoulLong,false,false>);

PairStyle(lj/cut/tip4p/long/kk,        PairLJCutCoulLongKokkos<LMPDeviceType,PairLJCutTIP4PLong,true,false>);
PairStyle(lj/cut/tip4p/long/kk/device, PairLJCutCoulLongKokkos<LMPDeviceType,PairLJCutTIP4PLong,true,false>);
PairStyle(lj/cut/tip4p/long/kk/host,   PairLJCutCoulLongKokkos<LMPHostType,  PairLJCutTIP4PLong,true,false>);

#ifdef LMP_FEP

  PairStyle(lj/cut/coul/long/soft/kk,        PairLJCutCoulLongKokkos<LMPDeviceType,PairLJCutCoulLongSoft,false,true>);
  PairStyle(lj/cut/coul/long/soft/kk/device, PairLJCutCoulLongKokkos<LMPDeviceType,PairLJCutCoulLongSoft,false,true>);
  PairStyle(lj/cut/coul/long/soft/kk/host,   PairLJCutCoulLongKokkos<LMPHostType,  PairLJCutCoulLongSoft,false,true>);

  PairStyle(lj/cut/tip4p/long/soft/kk,        PairLJCutCoulLongKokkos<LMPDeviceType,PairLJCutTIP4PLongSoft,true,true>);
  PairStyle(lj/cut/tip4p/long/soft/kk/device, PairLJCutCoulLongKokkos<LMPDeviceType,PairLJCutTIP4PLongSoft,true,true>);
  PairStyle(lj/cut/tip4p/long/soft/kk/host,   PairLJCutCoulLongKokkos<LMPHostType,  PairLJCutTIP4PLongSoft,true,true>);

#endif // LMP_FEP

// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LJ_CUT_COUL_LONG_KOKKOS_H
#define LMP_PAIR_LJ_CUT_COUL_LONG_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_cut_coul_long.h"

namespace LAMMPS_NS {

template<class DeviceType, class PairBase, bool TIP4P, bool SOFT>
class PairLJCutCoulLongKokkos : public PairKokkos<DeviceType,PairBase,true,TIP4P,SOFT>
{
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG = TIP4P ? COUL_TIP4P : COUL_LONG};
  PairLJCutCoulLongKokkos(class LAMMPS *);

};

}

#endif // !LMP_PAIR_LJ_CUT_COUL_LONG_KOKKOS_H
#endif

