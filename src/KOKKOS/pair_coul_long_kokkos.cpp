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

/* ----------------------------------------------------------------------
   Contributing author: Ray Shan (SNL)
------------------------------------------------------------------------- */

#include "pair_coul_long_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair_kokkos.h"
#include "respa.h"
#include "update.h"

#include "pair_coul_long.h"
#include "pair_tip4p_long.h"
#ifdef LMP_FEP
  #include "pair_coul_long_soft.h"
  #include "pair_tip4p_long_soft.h"
#endif

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace EwaldConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType, class PairCoulLongBase, bool TIP4P, bool SOFT>
PairCoulLongKokkos<DeviceType,PairCoulLongBase,TIP4P,SOFT>::PairCoulLongKokkos(LAMMPS *lmp):
  PairKokkos<DeviceType,PairCoulLongBase,false,TIP4P,SOFT>(lmp)
{


}

namespace LAMMPS_NS {

  // coul/long/kk
  template class PairCoulLongKokkos<LMPDeviceType,PairCoulLong,false,false>;

  // tip4p/long/kk
  template class PairCoulLongKokkos<LMPDeviceType,PairTIP4PLong,true,false>;

  #ifdef LMP_KOKKOS_GPU
    //template class PairCoulLongKokkos<LMPHostType,TIP4P,SOFT>;
  #endif
  
}



