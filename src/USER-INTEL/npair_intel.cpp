/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "npair_intel.h"
#include "nstencil.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairIntel::NPairIntel(LAMMPS *lmp) : NPair(lmp) {
  int ifix = modify->find_fix("package_intel");
  if (ifix < 0)
    error->all(FLERR,
               "The 'package intel' command is required for /intel styles");
  _fix = static_cast<FixIntel *>(modify->fix[ifix]);
  #ifdef _LMP_INTEL_OFFLOAD
  _cop = _fix->coprocessor_number();
  _off_map_stencil = 0;
  #endif
}

/* ---------------------------------------------------------------------- */

NPairIntel::~NPairIntel() {
  #ifdef _LMP_INTEL_OFFLOAD
  if (_off_map_stencil) {
    const int * stencil = this->stencil;
    #pragma offload_transfer target(mic:_cop)	\
      nocopy(stencil:alloc_if(0) free_if(1))
  }
  #endif
}

/* ---------------------------------------------------------------------- */

#ifdef _LMP_INTEL_OFFLOAD
void NPairIntel::grow_stencil()
{
  if (_off_map_stencil != stencil) {
    if (_off_map_stencil) {
      const int * stencil = _off_map_stencil;
      #pragma offload_transfer target(mic:_cop) \
        nocopy(stencil:alloc_if(0) free_if(1))
    }
    _off_map_stencil = stencil;
    const int * stencil = _off_map_stencil;
    const int maxstencil = ns->get_maxstencil();
    #pragma offload_transfer target(mic:_cop)	\
      in(stencil:length(maxstencil) alloc_if(1) free_if(0))
  }
}
#endif
