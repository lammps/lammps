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
------------------------------------------------------------------------- */

#include "pair_vashishta_table_kokkos.h"

#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairVashishtaTableKokkos<DeviceType>::PairVashishtaTableKokkos(LAMMPS *lmp) :
  PairVashishtaKokkos<DeviceType>(lmp)
{
  ntable = 0;
  this->tabflag = 1;
}

/* ----------------------------------------------------------------------
   global settings, same as for pair style vashishta/table
------------------------------------------------------------------------- */

template<class DeviceType>
void PairVashishtaTableKokkos<DeviceType>::settings(int narg, char **arg)
{
  if (narg != 2) this->error->all(FLERR,"Illegal pair_style command");

  ntable = utils::inumeric(FLERR,arg[0],false,this->lmp);
  this->tabinner = utils::numeric(FLERR,arg[1],false,this->lmp);

  if (this->tabinner <= 0.0)
    this->error->all(FLERR,"Illegal inner cutoff for tabulation");
}

/* ----------------------------------------------------------------------
   tabulate the two-body force and energy with the analytic expressions
   of the parent class, like PairVashishtaTable::create_tables()
------------------------------------------------------------------------- */

template<class DeviceType>
void PairVashishtaTableKokkos<DeviceType>::setup_params()
{
  PairVashishtaKokkos<DeviceType>::setup_params();

  const int nelements = this->nelements;
  const double cutmax = this->cutmax;
  this->tabinnersq = this->tabinner*this->tabinner;
  const double tabinnersq = this->tabinnersq;

  const double deltaR2 = (cutmax*cutmax - tabinnersq) / (ntable-1);
  const double oneOverDeltaR2 = 1.0/deltaR2;

  k_forceTable = DAT::tdual_kkfloat_3d("pair:vashishta:forceTable",nelements,nelements,ntable+1);
  k_potentialTable = DAT::tdual_kkfloat_3d("pair:vashishta:potentialTable",nelements,nelements,ntable+1);
  auto h_forceTable = k_forceTable.view_host();
  auto h_potentialTable = k_potentialTable.view_host();

  double rsq,fpair,eng;

  for (int i = 0; i < nelements; i++) {
    for (int j = 0; j < nelements; j++) {
      int ijparam = this->elem3param[i][j][j];
      for (int idx = 0; idx <= ntable; idx++) {
        rsq = tabinnersq + idx*deltaR2;
        this->PairVashishta::twobody(&this->params[ijparam],rsq,fpair,1,eng);
        h_forceTable(i,j,idx) = static_cast<KK_FLOAT>(fpair);
        h_potentialTable(i,j,idx) = static_cast<KK_FLOAT>(eng);
      }
    }
  }

  k_forceTable.modify_host();
  k_forceTable.template sync<DeviceType>();
  k_potentialTable.modify_host();
  k_potentialTable.template sync<DeviceType>();

  this->d_forceTable = k_forceTable.template view<DeviceType>();
  this->d_potentialTable = k_potentialTable.template view<DeviceType>();
  this->tabinnersq_kk = static_cast<KK_FLOAT>(tabinnersq);
  this->oneOverDeltaR2_kk = static_cast<KK_FLOAT>(oneOverDeltaR2);
}

/* ----------------------------------------------------------------------
   memory usage of tabulation arrays
------------------------------------------------------------------------- */

template<class DeviceType>
double PairVashishtaTableKokkos<DeviceType>::memory_usage()
{
  return (double)2*this->nelements*this->nelements*sizeof(KK_FLOAT)*(ntable+1);
}

namespace LAMMPS_NS {
template class PairVashishtaTableKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairVashishtaTableKokkos<LMPHostType>;
#endif
}
