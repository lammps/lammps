/* -*- c++ -*- ----------------------------------------------------------
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
   Contributing authors: Stefan Paquay (Brandeis University)
------------------------------------------------------------------------- */

#include "atom_masks.h"
#include "atom_kokkos.h"
#include "fix_enforce2d_kokkos.h"

using namespace LAMMPS_NS;


template <class DeviceType>
FixEnforce2DKokkos<DeviceType>::FixEnforce2DKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixEnforce2D(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read   = X_MASK | V_MASK | F_MASK | MASK_MASK;
  datamask_modify = X_MASK | V_MASK | F_MASK;
}


template <class DeviceType>
void FixEnforce2DKokkos<DeviceType>::setup(int vflag)
{
  post_force(vflag);
}


template <class DeviceType>
void FixEnforce2DKokkos<DeviceType>::post_force(int vflag)
{
  atomKK->sync(execution_space,datamask_read);
  atomKK->modified(execution_space,datamask_modify);

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();

  mask = atomKK->k_mask.view<DeviceType>();

  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  FixEnforce2DKokkosPostForceFunctor<DeviceType> functor(this);
  Kokkos::parallel_for(nlocal,functor);

  // Probably sync here again?
  atomKK->sync(execution_space,datamask_read);
  atomKK->modified(execution_space,datamask_modify);

  for (int m = 0; m < nfixlist; m++)
    flist[m]->enforce2d();


}


template <class DeviceType>
void FixEnforce2DKokkos<DeviceType>::post_force_item( int i ) const
{

  if (mask[i] & groupbit){
    v(i,2) = 0;
    x(i,2) = 0;
    f(i,2) = 0;

    // Add for omega, angmom, torque...
  }

}


template<class DeviceType>
void FixEnforce2DKokkos<DeviceType>::cleanup_copy()
{
  id = style = NULL;
  vatom = NULL;
}


namespace LAMMPS_NS {
template class FixEnforce2DKokkos<LMPDeviceType>;
#ifdef KOKKOS_HAVE_CUDA
template class FixEnforce2DKokkos<LMPHostType>;
#endif
}
