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
   (Kokkos version) LAMMPS development team
------------------------------------------------------------------------- */

#include "pair_coul_cut_global_kokkos.h"

#include "error.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   override: enforce exactly 2 arguments (no per-pair cutoff allowed)
------------------------------------------------------------------------- */

template<class DeviceType>
void PairCoulCutGlobalKokkos<DeviceType>::coeff(int narg, char **arg)
{
  if (narg != 2)
    this->error->all(FLERR,"Incorrect args for pair coefficients" + utils::errorurl(21));

  PairCoulCut::coeff(narg,arg);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void *PairCoulCutGlobalKokkos<DeviceType>::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"cut_coul") == 0) return (void *) &this->cut_global;
  dim = 2;
  if (strcmp(str,"scale") == 0) return (void *) this->scale;
  return nullptr;
}

namespace LAMMPS_NS {
template class PairCoulCutGlobalKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairCoulCutGlobalKokkos<LMPHostType>;
#endif
}
