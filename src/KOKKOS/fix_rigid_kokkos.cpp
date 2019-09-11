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

#include "fix_rigid_kokkos.h"
#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "domain.h"
#include "update.h"
#include "respa.h"
#include "modify.h"
#include "group.h"
#include "comm.h"
#include "random_mars.h"
#include "force.h"
#include "input.h"
#include "variable.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

#include "atom_masks.h"
#include "atom_kokkos.h"


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;


enum{SINGLE,MOLECULE,GROUP};
enum{NONE,XYZ,XY,YZ,XZ};
enum{ISO,ANISO,TRICLINIC};

#define MAXLINE 1024
#define CHUNK 1024
#define ATTRIBUTE_PERBODY 20

#define TOLERANCE 1.0e-6
#define EPSILON 1.0e-7

#define SINERTIA 0.4            // moment of inertia prefactor for sphere
#define EINERTIA 0.2            // moment of inertia prefactor for ellipsoid
#define LINERTIA (1.0/12.0)     // moment of inertia prefactor for line segment

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixRigidKokkos<DeviceType>::FixRigidKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixRigid(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | V_MASK | F_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK;
  datamask_modify = X_MASK | V_MASK;
  // Also need:
  // read:    vcm, xcm, torque, tflag
  // modify:  fcm, fflag, angmom, sum

  if (comm->me == 0) {
    fprintf(stderr, "Using fix rigid/kokkos!\n");
  }
}

template <class DeviceType>
FixRigidKokkos<DeviceType>::~FixRigidKokkos()
{
}


template <class DeviceType>
void FixRigidKokkos<DeviceType>::cleanup_copy()
{
  id = style = NULL;
  vatom = NULL;
}


/*                                            
int FixRigidKokkos::setmask();                      
void FixRigidKokkos::init();                        
void FixRigidKokkos::setup(int);                    
*/
template <class DeviceType>
void FixRigidKokkos<DeviceType>::initial_integrate(int vflag)
{
  atomKK->sync(execution_space, datamask_read);
  atomKK->modified(execution_space, datamask_modify);

  // Grab all arrays you need for initial_integrate:
  x = atomKK->k_x.view<DeviceType>();
  y = atomKK->k_y.view<DeviceType>();
  z = atomKK->k_z.view<DeviceType>();

  rmass = atomKK->k_rmass.view<DeviceType>();
  mass = atomKK->k_mass.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  if (rmass.data()) {
    FixRigidKokkosInitialIntegrateFunctor<DeviceType,1> functor(this);
    Kokkos::parallel_for(nlocal,functor);
  } else {
    FixRigidKokkosInitialIntegrateFunctor<DeviceType,0> functor(this);
    Kokkos::parallel_for(nlocal,functor);
  }




  // set_xv is called in initial_integrate so make sure everything it needs is
  // Kokkos-able?
}

/*
template <class DeviceType>
void FixRigidKokkos::post_force(int)
{
}


template <class DeviceType>
void FixRigidKokkos::final_integrate()
{
}
*/

/*
void FixRigidKokkos::initial_integrate_respa(int, int, int);
void FixRigidKokkos::final_integrate_respa(int, int);       
void FixRigidKokkos::write_restart_file(char *);            
double FixRigidKokkos::compute_scalar();            
*/


namespace LAMMPS_NS {
template class FixRigidKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class FixRigidKokkos<LMPHostType>;
#endif
}
