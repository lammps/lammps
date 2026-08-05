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
   Contributing authors: Xiaowng Zhou (Sandia)
------------------------------------------------------------------------- */

#include "pair_eam_he_kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template <class DeviceType>
PairEAMHEKokkos<DeviceType>::PairEAMHEKokkos(LAMMPS *lmp) : PairEAMKokkos<DeviceType>(lmp)
{
  this->fileformat = PairEAM::FS;
  this->one_coeff = 1;
  this->he_flag = 1;
}

namespace LAMMPS_NS {
template class PairEAMHEKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairEAMHEKokkos<LMPHostType>;
#endif
}    // namespace LAMMPS_NS
