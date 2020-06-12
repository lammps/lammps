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

#include "fix_enforce2d_kokkos.h"
#include "atom_masks.h"
#include "atom_kokkos.h"
#include "comm.h"
#include "error.h"


using namespace LAMMPS_NS;


template <ExecutionSpace Space>
FixEnforce2DKokkos<Space>::FixEnforce2DKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixEnforce2D(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = Space;

  datamask_read   = V_MASK | F_MASK | OMEGA_MASK | MASK_MASK
          | TORQUE_MASK | ANGMOM_MASK;

  datamask_modify = V_MASK | F_MASK | OMEGA_MASK
          | TORQUE_MASK | ANGMOM_MASK;
}


template <ExecutionSpace Space>
void FixEnforce2DKokkos<Space>::setup(int vflag)
{
  post_force(vflag);
}


template <ExecutionSpace Space>
void FixEnforce2DKokkos<Space>::post_force(int vflag)
{
  atomKK->sync(Space,datamask_read);

  v = DualViewHelper<Space>::view(atomKK->k_v);
  f = DualViewHelper<Space>::view(atomKK->k_f);

  if (atomKK->omega_flag)
    omega  = DualViewHelper<Space>::view(atomKK->k_omega);

  if (atomKK->angmom_flag)
    angmom = DualViewHelper<Space>::view(atomKK->k_angmom);

  if (atomKK->torque_flag)
    torque = DualViewHelper<Space>::view(atomKK->k_torque);


  mask = DualViewHelper<Space>::view(atomKK->k_mask);

  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  int flag_mask = 0;
  if (atomKK->omega_flag) flag_mask  |= 1;
  if (atomKK->angmom_flag) flag_mask |= 2;
  if (atomKK->torque_flag) flag_mask |= 4;

  copymode = 1;
  switch( flag_mask ){
    case 0:{
      FixEnforce2DKokkosPostForceFunctor<Space,0,0,0> functor(this);
      Kokkos::parallel_for(nlocal,functor);
      break;
    }
    case 1:{
      FixEnforce2DKokkosPostForceFunctor<Space,1,0,0> functor(this);
      Kokkos::parallel_for(nlocal,functor);
      break;
    }
    case 2:{
      FixEnforce2DKokkosPostForceFunctor<Space,0,1,0> functor(this);
      Kokkos::parallel_for(nlocal,functor);
      break;
    }
    case 3:{
      FixEnforce2DKokkosPostForceFunctor<Space,1,1,0> functor(this);
      Kokkos::parallel_for(nlocal,functor);
      break;
    }
    case 4:{
      FixEnforce2DKokkosPostForceFunctor<Space,0,0,1> functor(this);
      Kokkos::parallel_for(nlocal,functor);
      break;
    }
    case 5:{
      FixEnforce2DKokkosPostForceFunctor<Space,1,0,1> functor(this);
      Kokkos::parallel_for(nlocal,functor);
      break;
    }
    case 6:{
      FixEnforce2DKokkosPostForceFunctor<Space,0,1,1> functor(this);
      Kokkos::parallel_for(nlocal,functor);
      break;
    }
    case 7:{
      FixEnforce2DKokkosPostForceFunctor<Space,1,1,1> functor(this);
      Kokkos::parallel_for(nlocal,functor);
      break;
    }
    default:
      error->all(FLERR, "Flag in fix_enforce2d_kokkos outside of what it should be");
  }
  copymode = 0;

  atomKK->modified(Space,datamask_modify);

  for (int m = 0; m < nfixlist; m++) {
    atomKK->sync(flist[m]->execution_space,flist[m]->datamask_read);
    flist[m]->enforce2d();
    atomKK->modified(flist[m]->execution_space,flist[m]->datamask_modify);
  }

}


template <ExecutionSpace Space>
template <int omega_flag, int angmom_flag, int torque_flag>
KOKKOS_INLINE_FUNCTION
void FixEnforce2DKokkos<Space>::post_force_item( int i ) const
{
  if (mask[i] & groupbit){
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


namespace LAMMPS_NS {
template class FixEnforce2DKokkos<Device>;
template class FixEnforce2DKokkos<Host>;
}
