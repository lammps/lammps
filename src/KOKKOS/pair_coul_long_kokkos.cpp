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

/*


template<class DeviceType, class PairCoulLongBase, bool TIP4P, bool SOFT>
void PairCoulLongKokkos<DeviceType,PairCoulLongBase,TIP4P,SOFT>::allocate()
{
  PairCoulLong::allocate();

  int n = atom->ntypes;

  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  d_cut_coulsq = typename AT::t_kkfloat_2d("pair:cut_coulsq",n+1,n+1);

  k_params = Kokkos::DualView<params_coul**,Kokkos::LayoutRight,DeviceType>("PairCoulLong::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}



template<class DeviceType, class PairCoulLongBase, bool TIP4P, bool SOFT>
void PairCoulLongKokkos<DeviceType,PairCoulLongBase,TIP4P,SOFT>::init_style()
{
  PairCoulLong::init_style();

  Kokkos::deep_copy(d_cut_coulsq,cut_coulsq);
  Kokkos::deep_copy(d_cut_ljsq,cut_coulsq);

  // error if rRESPA with inner levels

  if (update->whichflag == 1 && utils::strmatch(update->integrate_style,"^respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    if (respa)
      error->all(FLERR,"Cannot use Kokkos pair style with rRESPA inner/middle");
  }

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();
}

template<class DeviceType, class PairCoulLongBase, bool TIP4P, bool SOFT>
double PairCoulLongKokkos<DeviceType,PairCoulLongBase,TIP4P,SOFT>::init_one(int i, int j)
{
  double cutone = PairCoulLong::init_one(i,j);

  k_params.view_host()(i,j).cut_coulsq = cut_coulsq;

  k_params.view_host()(j,i) = k_params.view_host()(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = cutone*cutone;
    m_cut_coulsq[j][i] = m_cut_coulsq[i][j] = cut_coulsq;
    m_cut_ljsq[j][i] = m_cut_ljsq[i][j] = cut_coulsq;
  }

  k_cutsq.view_host()(i,j) = cutone*cutone;
  k_cutsq.modify_host();
  k_params.modify_host();

  return cutone;
}

*/

namespace LAMMPS_NS {

// coul/long/kk
template class PairCoulLongKokkos<LMPDeviceType,PairCoulLong,false,false>;

// tip4p/long/kk
//template class PairCoulLongKokkos<LMPDeviceType,PairTIP4PLong,true,false>;

#ifdef LMP_FEP

  // coul/long/soft/kk
  //template class PairCoulLongKokkos<LMPDeviceType,false,true>;

  // tip4p/long/soft/kk
  //template class PairCoulLongKokkos<LMPDeviceType,true,true>;

#endif // LMP_FEP

#ifdef LMP_KOKKOS_GPU
//template class PairCoulLongKokkos<LMPHostType,TIP4P,SOFT>;
#endif
}



