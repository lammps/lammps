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
   Contributing author: Claude Opus 5.5 (Anthropic), under the direction of
   Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_gw_zbl_kokkos.h"

// the implementation is in pair_gw_kokkos.h, shared with pair style gw/kk

namespace LAMMPS_NS {
template class PairGWZBLKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairGWZBLKokkos<LMPHostType>;
#endif
}
