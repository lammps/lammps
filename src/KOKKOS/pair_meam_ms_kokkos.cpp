/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_meam_ms_kokkos.h"
#include "meam.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
template<class DeviceType>
PairMEAMMSKokkos<DeviceType>::PairMEAMMSKokkos(LAMMPS *lmp) : PairMEAMKokkos<DeviceType>(lmp)
{
  this->meam_inst->msmeamflag = this->msmeamflag = 1;
  this->myname = "meam/ms/kk";
}

namespace LAMMPS_NS {
template class PairMEAMMSKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairMEAMMSKokkos<LMPHostType>;
#endif
}
