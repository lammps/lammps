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

#include <cstdlib>
#include <cstring>
#include "fix_wall_reflect_kokkos.h"
#include "atom_kokkos.h"
#include "comm.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "lattice.h"
#include "input.h"
#include "variable.h"
#include "error.h"
#include "force.h"
#include "atom_masks.h"


using namespace LAMMPS_NS;
using namespace FixConst;

enum{XLO=0,XHI=1,YLO=2,YHI=3,ZLO=4,ZHI=5};
enum{NONE=0,EDGE,CONSTANT,VARIABLE};


/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixWallReflectKokkos<DeviceType>::FixWallReflectKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixWallReflect(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | V_MASK | MASK_MASK;
  datamask_modify = X_MASK | V_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixWallReflectKokkos<DeviceType>::post_integrate()
{
  // coord = current position of wall
  // evaluate variable if necessary, wrap with clear/add

  atomKK->sync(execution_space,datamask_read);
  atomKK->modified(execution_space,datamask_modify);

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;


  if (varflag) modify->clearstep_compute();

  for (int m = 0; m < nwall; m++) {
    if (wallstyle[m] == VARIABLE) {
      coord = input->variable->compute_equal(varindex[m]);
      if (wallwhich[m] < YLO) coord *= xscale;
      else if (wallwhich[m] < ZLO) coord *= yscale;
      else coord *= zscale;
    } else coord = coord0[m];

    dim = wallwhich[m] / 2;
    side = wallwhich[m] % 2;

    copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixWallReflectPostIntegrate>(0,nlocal),*this);
    copymode = 0;
  }

  if (varflag) modify->addstep_compute(update->ntimestep + 1);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixWallReflectKokkos<DeviceType>::operator()(TagFixWallReflectPostIntegrate, const int &i) const {
  if (mask[i] & groupbit) {
    if (side == 0) {
      if (x(i,dim) < coord) {
        x(i,dim) = coord + (coord - x(i,dim));
        v(i,dim) = -v(i,dim);
      }
    } else {
      if (x(i,dim) > coord) {
        x(i,dim) = coord - (x(i,dim) - coord);
        v(i,dim) = -v(i,dim);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixWallReflectKokkos<LMPDeviceType>;
#ifdef KOKKOS_HAVE_CUDA
template class FixWallReflectKokkos<LMPHostType>;
#endif
}

