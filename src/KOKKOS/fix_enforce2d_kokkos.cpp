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
   Contributing authors: Stefan Paquay & Matthew Peterson (Brandeis University)
------------------------------------------------------------------------- */

#include "atom_masks.h"
#include "atom_kokkos.h"
#include "comm.h"
#include "error.h"
#include "fix_enforce2d_kokkos.h"


using namespace LAMMPS_NS;


template <class DeviceType>
FixEnforce2DKokkos<DeviceType>::FixEnforce2DKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixEnforce2D(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read   = X_MASK | V_MASK | F_MASK | OMEGA_MASK | MASK_MASK
	  | TORQUE_MASK | ANGMOM_MASK; // | */ // MASK_MASK;

  datamask_modify = X_MASK | V_MASK | F_MASK | OMEGA_MASK
	  | TORQUE_MASK | ANGMOM_MASK;
}


template <class DeviceType>
void FixEnforce2DKokkos<DeviceType>::setup(int vflag)
{
  if( comm->me == 0 ){
    fprintf(screen, "omega, angmom and torque flags are %d, %d, %d\n",
            atomKK->omega_flag, atomKK->angmom_flag, atomKK->torque_flag );
  }
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

  if( atomKK->omega_flag )
    omega  = atomKK->k_omega.view<DeviceType>();

  if( atomKK->angmom_flag )
    angmom = atomKK->k_angmom.view<DeviceType>();

  if( atomKK->torque_flag )
    torque = atomKK->k_torque.view<DeviceType>();


  mask = atomKK->k_mask.view<DeviceType>();

  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  int flag_mask = 0;
  if( atomKK->omega_flag ) flag_mask  |= 1;
  if( atomKK->angmom_flag ) flag_mask |= 2;
  if( atomKK->torque_flag ) flag_mask |= 4;

  switch( flag_mask ){
    case 0:{
      FixEnforce2DKokkosPostForceFunctor<DeviceType,0,0,0> functor(this);
      Kokkos::parallel_for(nlocal,functor);
      break;
    }
    case 1:{
      FixEnforce2DKokkosPostForceFunctor<DeviceType,1,0,0> functor(this);
      Kokkos::parallel_for(nlocal,functor);
      break;
    }
    case 2:{
      FixEnforce2DKokkosPostForceFunctor<DeviceType,0,1,0> functor(this);
      Kokkos::parallel_for(nlocal,functor);
      break;
    }
    case 3:{
      FixEnforce2DKokkosPostForceFunctor<DeviceType,1,1,0> functor(this);
      Kokkos::parallel_for(nlocal,functor);
      break;
    }
    case 4:{
      FixEnforce2DKokkosPostForceFunctor<DeviceType,0,0,1> functor(this);
      Kokkos::parallel_for(nlocal,functor);
      break;
    }
    case 5:{
      FixEnforce2DKokkosPostForceFunctor<DeviceType,1,0,1> functor(this);
      Kokkos::parallel_for(nlocal,functor);
      break;
    }
    case 6:{
      FixEnforce2DKokkosPostForceFunctor<DeviceType,0,1,1> functor(this);
      Kokkos::parallel_for(nlocal,functor);
      break;
    }
    case 7:{
      FixEnforce2DKokkosPostForceFunctor<DeviceType,1,1,1> functor(this);
      Kokkos::parallel_for(nlocal,functor);
      break;
    }
    default:
      error->all(FLERR, "flag_mask outside of what it should be");
  }


  // Probably sync here again?
  atomKK->sync(execution_space,datamask_read);
  atomKK->modified(execution_space,datamask_modify);

  for (int m = 0; m < nfixlist; m++)
    flist[m]->enforce2d();
}


template <class DeviceType>
template <int omega_flag, int angmom_flag, int torque_flag>
void FixEnforce2DKokkos<DeviceType>::post_force_item( int i ) const
{
  if (mask[i] & groupbit){
    // x(i,2) = 0; // Enforce2d does not set x[2] to zero either... :/
    v(i,2) = 0.0;
    f(i,2) = 0.0;

    if(omega_flag){
      omega(i,0) = 0.0;
      omega(i,1) = 0.0;
    }

    if(angmom_flag){
      angmom(i,0) = 0.0;
      angmom(i,1) = 0.0;
    }

    if(torque_flag){
      torque(i,0) = 0.0;
      torque(i,1) = 0.0;
    }
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
