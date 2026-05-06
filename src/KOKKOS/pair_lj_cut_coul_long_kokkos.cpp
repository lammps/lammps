// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_lj_cut_coul_long_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include "pair_lj_cut_coul_long.h"
#include "pair_lj_cut_tip4p_long.h"
#include "pair_lj_cut_coul_long_soft.h"
#include "pair_lj_cut_tip4p_long_soft.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType, class PairBase, bool TIP4P, bool SOFT>
PairLJCutCoulLongKokkos<DeviceType,PairBase,TIP4P,SOFT>::PairLJCutCoulLongKokkos(LAMMPS *lmp):
  PairKokkos<DeviceType,PairBase,true,TIP4P,SOFT>(lmp)
{



}

namespace LAMMPS_NS {

  // lj/cut/coul/long/kk
  template class PairLJCutCoulLongKokkos<LMPDeviceType,PairLJCutCoulLong,false,false>;

  // lj/cut/tip4p/long/kk
  template class PairLJCutCoulLongKokkos<LMPDeviceType,PairLJCutTIP4PLong,true,false>;

  // lj/cut/coul/long/soft/kk
  template class PairLJCutCoulLongKokkos<LMPDeviceType,PairLJCutCoulLongSoft,false,true>;

  // lj/cut/tip4p/long/soft/kk
  template class PairLJCutCoulLongKokkos<LMPDeviceType,PairLJCutTIP4PLongSoft,true,true>;

  #ifdef LMP_KOKKOS_GPU
    //template class PairLJCutCoulLongKokkos<LMPHostType>;
  #endif

}

