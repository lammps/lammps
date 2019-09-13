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


template <class DeviceType>
void FixRigidKokkos<DeviceType>::init()
{
  FixRigid::init();
  atomKK->k_mass.modify<LMPHostType>();
  atomKK->k_mass.sync<DeviceType>();
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
  double dtfm;

  
  {
    // Local block for Kokkos parallel for:
    auto l_vcm = vcm;
    auto l_xcm = xcm;
    auto l_angmom = angmom;
    auto l_masstotal = masstotal;
    auto l_fflag = fflag;
    auto l_fcm = fcm;
    auto l_torque = torque;
    auto l_tflag = tflag;
    
    Kokkos::parallel_for(nbody, LAMMPS_LAMBDA(const int& ibody) {

      const double dtfm = dtf / l_masstotal[ibody];
        
      l_vcm(ibody,0) += dtfm * l_fcm(ibody,0) * l_fflag(ibody,0);
      l_vcm(ibody,1) += dtfm * l_fcm(ibody,1) * l_fflag(ibody,1);
      l_vcm(ibody,2) += dtfm * l_fcm(ibody,2) * l_fflag(ibody,2);
      
      // update xcm by full step
      
      l_xcm(ibody,0) += dtv * l_vcm(ibody,0);
      l_xcm(ibody,1) += dtv * l_vcm(ibody,1);
      l_xcm(ibody,2) += dtv * l_vcm(ibody,2);
      
      // update angular momentum by 1/2 step
      
      l_angmom(ibody,0) += dtf * l_torque(ibody,0) * l_tflag(ibody,0);
      l_angmom(ibody,1) += dtf * l_torque(ibody,1) * l_tflag(ibody,1);
      l_angmom(ibody,2) += dtf * l_torque(ibody,2) * l_tflag(ibody,2);

      // TODO: Convert the body angmom back to per-atom velocities.
    });
  }
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
