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

#include "fix_external_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain_kokkos.h"
#include "error.h"
#include "kokkos_base.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixExternalKokkos<DeviceType>::FixExternalKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixExternal(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  memory->destroy(fexternal);
  grow_arrays(atom->nmax);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixExternalKokkos<DeviceType>::~FixExternalKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_fexternal,fexternal);
  fexternal = nullptr;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixExternalKokkos<DeviceType>::init()
{
  FixExternal::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use respa with Kokkos");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixExternalKokkos<DeviceType>::post_force(int vflag)
{
  atomKK->sync(execution_space, F_MASK | MASK_MASK);

  f = atomKK->k_f.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;
  bigint ntimestep = update->ntimestep;

  int eflag = eflag_caller;
  ev_init(eflag,vflag);

  // invoke the callback in driver program
  // it will fill fexternal with forces
  // base class will move them along with atoms if ncall != napply

  if ((mode == PF_CALLBACK) && (ntimestep % ncall == 0)) {
    atomKK->sync(Host, X_MASK | TAG_MASK);
    (this->callback)(ptr_caller,update->ntimestep,
                     atom->nlocal,atom->tag,atom->x,fexternal);
  }

  // add forces from current fexternal to KOKKOS array and then to atoms in group

  if ((ntimestep % napply) == 0) {
    // transfer external force data to device
    k_fexternal.modify_host();
    k_fexternal.sync<DeviceType>();

    // apply external forces
    copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixExternal>(0,nlocal),*this);
    copymode = 0;
    atomKK->modified(execution_space, F_MASK);

    // add contribution to global virial from previously stored value

    if (vflag_global)
      for (int i = 0; i < 6; ++i)
        virial[i] = user_virial[i];
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixExternalKokkos<DeviceType>::operator()(TagFixExternal, const int &i) const {
  if (mask[i] & groupbit) {
    f(i,0) += d_fexternal(i,0);
    f(i,1) += d_fexternal(i,1);
    f(i,2) += d_fexternal(i,2);
  }
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

template<class DeviceType>
void FixExternalKokkos<DeviceType>::grow_arrays(int nmax)
{
  memoryKK->grow_kokkos(k_fexternal,fexternal,nmax,3,"external:fexternal");
  d_fexternal = k_fexternal.view<DeviceType>();
  memset(&fexternal[0][0], 0, sizeof(double)*3*nmax);
  array_atom = fexternal;
}

namespace LAMMPS_NS {
template class FixExternalKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixExternalKokkos<LMPHostType>;
#endif
}
