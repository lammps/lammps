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
   Contributing authors: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "pair_eam_fs_kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template <class DeviceType>
PairEAMFSKokkos<DeviceType>::PairEAMFSKokkos(LAMMPS *lmp) : PairEAMKokkos<DeviceType>(lmp)
{
  this->fileformat = PairEAM::FS;
  this->one_coeff = 1;
}

namespace LAMMPS_NS {
template class PairEAMFSKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairEAMFSKokkos<LMPHostType>;
#endif
}    // namespace LAMMPS_NS
