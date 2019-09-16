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
#include "rigid_const.h"

#include "atom_masks.h"
#include "atom_kokkos.h"
#include "memory_kokkos.h"


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace RigidConst;

constexpr const bool debug_output = true;

/* ---------------------------------------------------------------------- */

// Some MathExtra functions had to be ported to Kokkos. In the future we
// probably want to collect all of them somewhere.
namespace MathExtraKokkos {

// Trigger warning about redundant inline on explicit specialization :(
KOKKOS_INLINE_FUNCTION 
template <class DeviceType>
void angmom_to_omega(typename ArrayTypes<DeviceType>::t_v_array m,
                     typename ArrayTypes<DeviceType>::t_x_array ex,
                     typename ArrayTypes<DeviceType>::t_x_array ey,
                     typename ArrayTypes<DeviceType>::t_x_array ez,
                     typename ArrayTypes<DeviceType>::t_x_array idiag,
                     typename ArrayTypes<DeviceType>::t_v_array w, int ibody)
{
  double wbody[3];

  if (idiag(ibody,0) == 0.0) wbody[0] = 0.0;
  else wbody[0] = (m(ibody,0)*ex(ibody,0) + m(ibody,1)*ex(ibody,1) + m(ibody,2)*ex(ibody,2)) / idiag(ibody,0);
  if (idiag(ibody,1) == 0.0) wbody[1] = 0.0;
  else wbody[1] = (m(ibody,0)*ey(ibody,0) + m(ibody,1)*ey(ibody,1) + m(ibody,2)*ey(ibody,2)) / idiag(ibody,1);
  if (idiag(ibody,2) == 0.0) wbody[2] = 0.0;
  else wbody[2] = (m(ibody,0)*ez(ibody,0) + m(ibody,1)*ez(ibody,1) + m(ibody,2)*ez(ibody,2)) / idiag(ibody,2);
  
  w(ibody,0) = wbody[0]*ex(ibody,0) + wbody[1]*ey(ibody,0) + wbody[2]*ez(ibody,0);
  w(ibody,1) = wbody[0]*ex(ibody,1) + wbody[1]*ey(ibody,1) + wbody[2]*ez(ibody,1);
  w(ibody,2) = wbody[0]*ex(ibody,2) + wbody[1]*ey(ibody,2) + wbody[2]*ez(ibody,2);
}


template <typename a_arr_type, typename b_arr_type>
inline void vecquat(a_arr_type a, b_arr_type b, double c[4])
{
  c[0] = -a[0]*b[1] - a[1]*b[2] - a[2]*b[3];
  c[1] = b[0]*a[0] + a[1]*b[3] - a[2]*b[2];
  c[2] = b[0]*a[1] + a[2]*b[1] - a[0]*b[3];
  c[3] = b[0]*a[2] + a[0]*b[2] - a[1]*b[1];
}

/* ----------------------------------------------------------------------
   matrix times vector
------------------------------------------------------------------------- */

template <typename arr_out_type>
inline void matvec(const double m[3][3], const double v[3],
                   arr_out_type ans)
{
  ans[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];
  ans[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];
  ans[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];
}
	

/* ----------------------------------------------------------------------
   transposed matrix times vector
------------------------------------------------------------------------- */

template <typename v_arr_type, typename arr_out_type>
inline void transpose_matvec(const double m[3][3], v_arr_type v,
                             arr_out_type ans)
{
  ans[0] = m[0][0]*v[0] + m[1][0]*v[1] + m[2][0]*v[2];
  ans[1] = m[0][1]*v[0] + m[1][1]*v[1] + m[2][1]*v[2];
  ans[2] = m[0][2]*v[0] + m[1][2]*v[1] + m[2][2]*v[2];
}



/* ----------------------------------------------------------------------
   normalize a quaternion
------------------------------------------------------------------------- */

template <typename q_arr_type>
inline void qnormalize(q_arr_type q)
{
  double norm = 1.0 / sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
  q[0] *= norm;
  q[1] *= norm;
  q[2] *= norm;
  q[3] *= norm;
}

/* ----------------------------------------------------------------------
   conjugate of a quaternion: qc = conjugate of q
   assume q is of unit length
------------------------------------------------------------------------- */

inline void qconjugate(double q[4], double qc[4])
{
  qc[0] = q[0];
  qc[1] = -q[1];
  qc[2] = -q[2];
  qc[3] = -q[3];
}

/* ----------------------------------------------------------------------
   compute rotation matrix from quaternion
   quat = [w i j k]
------------------------------------------------------------------------- */

void quat_to_mat(const double quat[4], double mat[3][3])
{
  double w2 = quat[0]*quat[0];
  double i2 = quat[1]*quat[1];
  double j2 = quat[2]*quat[2];
  double k2 = quat[3]*quat[3];
  double twoij = 2.0*quat[1]*quat[2];
  double twoik = 2.0*quat[1]*quat[3];
  double twojk = 2.0*quat[2]*quat[3];
  double twoiw = 2.0*quat[1]*quat[0];
  double twojw = 2.0*quat[2]*quat[0];
  double twokw = 2.0*quat[3]*quat[0];

  mat[0][0] = w2+i2-j2-k2;
  mat[0][1] = twoij-twokw;
  mat[0][2] = twojw+twoik;

  mat[1][0] = twoij+twokw;
  mat[1][1] = w2-i2+j2-k2;
  mat[1][2] = twojk-twoiw;

  mat[2][0] = twoik-twojw;
  mat[2][1] = twojk+twoiw;
  mat[2][2] = w2-i2-j2+k2;
}



/* ----------------------------------------------------------------------
   compute space-frame ex,ey,ez from current quaternion q
   ex,ey,ez = space-frame coords of 1st,2nd,3rd principal axis
   operation is ex = q' d q = Q d, where d is (1,0,0) = 1st axis in body frame
------------------------------------------------------------------------- */
template <typename q_arr_type>
inline void q_to_exyz(q_arr_type q, double *ex, double *ey, double *ez)
{
  double q0 = q[0];
  double q1 = q[1];
  double q2 = q[2];
  double q3 = q[3];

  ex[0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
  ex[1] = 2.0 * (q1*q2 + q0*q3);
  ex[2] = 2.0 * (q1*q3 - q0*q2);

  ey[0] = 2.0 * (q1*q2 - q0*q3);
  ey[1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
  ey[2] = 2.0 * (q2*q3 + q0*q1);

  ez[0] = 2.0 * (q1*q3 + q0*q2);
  ez[1] = 2.0 * (q2*q3 - q0*q1);
  ez[2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;
}




/* ----------------------------------------------------------------------
   compute omega from angular momentum
   w = omega = angular velocity in space frame
   wbody = angular velocity in body frame
   project space-frame angular momentum onto body axes
     and divide by principal moments
------------------------------------------------------------------------- */


template <typename v_arr_type, typename x_arr_type>
inline void mq_to_omega(v_arr_type m, double q[4],
                        x_arr_type moments, v_arr_type w)
{
  double wbody[3];
  double rot[3][3];

  MathExtraKokkos::quat_to_mat(q,rot);
  MathExtraKokkos::transpose_matvec(rot,m,wbody);
  if (moments[0] == 0.0) wbody[0] = 0.0;
  else wbody[0] /= moments[0];
  if (moments[1] == 0.0) wbody[1] = 0.0;
  else wbody[1] /= moments[1];
  if (moments[2] == 0.0) wbody[2] = 0.0;
  else wbody[2] /= moments[2];
  MathExtraKokkos::matvec(rot,wbody,w);
}



// TODO: PORT THIS FURTHER
KOKKOS_INLINE_FUNCTION
template <class DeviceType>
void richardson(typename ArrayTypes<DeviceType>::t_x_array q,
                typename ArrayTypes<DeviceType>::t_v_array m,
                typename ArrayTypes<DeviceType>::t_v_array w,
                typename ArrayTypes<DeviceType>::t_x_array moments,
                double dtq, int ibody)
{
  // full update from dq/dt = 1/2 w q

  double wq[4];
  //MathExtraKokkos::vecquat<DeviceType>(w(ibody),q(ibody),wq);

  double qfull[4];
  qfull[0] = q(ibody,0) + dtq * wq[0];
  qfull[1] = q(ibody,1) + dtq * wq[1];
  qfull[2] = q(ibody,2) + dtq * wq[2];
  qfull[3] = q(ibody,3) + dtq * wq[3];
  MathExtraKokkos::qnormalize(qfull);

  // 1st half update from dq/dt = 1/2 w q

  double qhalf[4];
  qhalf[0] = q(ibody,0) + 0.5*dtq * wq[0];
  qhalf[1] = q(ibody,1) + 0.5*dtq * wq[1];
  qhalf[2] = q(ibody,2) + 0.5*dtq * wq[2];
  qhalf[3] = q(ibody,3) + 0.5*dtq * wq[3];
  MathExtraKokkos::qnormalize(qhalf);

  // re-compute omega at 1/2 step from m at 1/2 step and q at 1/2 step
  // recompute wq

  auto m_ibody = Kokkos::subview(m, ibody, Kokkos::ALL);
  auto moments_ibody = Kokkos::subview(moments, ibody, Kokkos::ALL);
  auto w_ibody = Kokkos::subview(w, ibody, Kokkos::ALL);
  auto q_ibody = Kokkos::subview(q, ibody, Kokkos::ALL);
  
  MathExtraKokkos::mq_to_omega(m_ibody,qhalf,moments_ibody,w_ibody);
  MathExtraKokkos::vecquat(w_ibody,qhalf,wq);

  // 2nd half update from dq/dt = 1/2 w q

  qhalf[0] += 0.5*dtq * wq[0];
  qhalf[1] += 0.5*dtq * wq[1];
  qhalf[2] += 0.5*dtq * wq[2];
  qhalf[3] += 0.5*dtq * wq[3];
  MathExtraKokkos::qnormalize(qhalf);

  // corrected Richardson update

  q(ibody,0) = 2.0*qhalf[0] - qfull[0];
  q(ibody,1) = 2.0*qhalf[1] - qfull[1];
  q(ibody,2) = 2.0*qhalf[2] - qfull[2];
  q(ibody,3) = 2.0*qhalf[3] - qfull[3];
  MathExtraKokkos::qnormalize(q_ibody);
}
	

} // MathExtraKokkos



template <class DeviceType>
FixRigidKokkos<DeviceType>::FixRigidKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixRigid(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read   = (X_MASK | V_MASK | F_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK
                     | OMEGA_MASK | ANGMOM_MASK | TORQUE_MASK);
  datamask_modify = (X_MASK | V_MASK | OMEGA_MASK | ANGMOM_MASK);

  // Only some arrays are initialized in the FixRigid constructor that we need.
  // The following are all either
  //   a) set after that, so we need not worry about preserving the data, or
  //   b) Set to all zeros which is also what Kokkos does by default.

  // These are all zero-initialized in FixRigid::FixRigid(...) and
  // later initialized:
  memoryKK->create_kokkos(k_xcm,xcm,nbody,3,"rigid/kk:xcm");
  memoryKK->create_kokkos(k_vcm,vcm,nbody,3,"rigid/kk:vcm");
  memoryKK->create_kokkos(k_fcm,fcm,nbody,3,"rigid/kk:fcm");

  memoryKK->create_kokkos(k_torque,torque,nbody,3,"rigid/kk:torque");


  // These are initialized in the base c-tor but not to 0:
  memoryKK->create_kokkos(k_fflag,fflag,nbody,3,"rigid/kk:fflag");
  memoryKK->create_kokkos(k_tflag,tflag,nbody,3,"rigid/kk:tflag");

  // The initialization for these two is simple enough:
  for (int i = 0; i < nbody; i++) {
	  k_fflag.d_view(i,0) = k_fflag.d_view(i,1) = k_fflag.d_view(i,2) = 1.0;
	  k_tflag.d_view(i,0) = k_tflag.d_view(i,1) = k_tflag.d_view(i,2) = 1.0;
	  if (domain->dimension == 2) k_fflag.d_view(i,2) = k_tflag.d_view(i,2) = 0.0;
  }
  k_fflag.template modify<LMPDeviceType>();
  k_tflag.template modify<LMPDeviceType>();

  // These seem fine.
  memoryKK->create_kokkos(k_omega, omega, nbody,3,"rigid/kk:omega");
  memoryKK->create_kokkos(k_angmom,angmom,nbody,3,"rigid/kk:angmom");

  
  // These are allocated in the base c-tor but set later:
  memoryKK->create_kokkos(k_quat,   quat,    nbody, 4,"rigid/kk:quat");
  memoryKK->create_kokkos(k_inertia,inertia, nbody, 3,"rigid/kk:inertia");

  memoryKK->create_kokkos(k_ex_space,ex_space, nbody, 3,"rigid/kk:ex_space");
  memoryKK->create_kokkos(k_ey_space,ey_space, nbody, 3,"rigid/kk:ey_space");
  memoryKK->create_kokkos(k_ez_space,ez_space, nbody, 3,"rigid/kk:ez_space");


  if (debug_output && comm->me == 0) {
    fprintf(stderr, "Using fix rigid/kokkos!\n");
  }
}

template <class DeviceType>
FixRigidKokkos<DeviceType>::~FixRigidKokkos()
{
  memoryKK->destroy_kokkos(k_xcm);
  memoryKK->destroy_kokkos(k_vcm);
  memoryKK->destroy_kokkos(k_fcm);

  memoryKK->destroy_kokkos(k_tflag);
  memoryKK->destroy_kokkos(k_fflag);

  memoryKK->destroy_kokkos(k_omega);
  memoryKK->destroy_kokkos(k_angmom);
  memoryKK->destroy_kokkos(k_torque);

  memoryKK->destroy_kokkos(k_quat);
  memoryKK->destroy_kokkos(k_inertia);

  memoryKK->destroy_kokkos(k_ex_space);
  memoryKK->destroy_kokkos(k_ey_space);
  memoryKK->destroy_kokkos(k_ez_space);
}


template <class DeviceType>
void FixRigidKokkos<DeviceType>::cleanup_copy()
{
  id = style = NULL;
  vatom = NULL;


}



template <class DeviceType>
template <typename kokkos_arr, typename base_arr>
void FixRigidKokkos<DeviceType>::debug_print(kokkos_arr k_arr, base_arr arr,
                                             const char *name, int idx)
{
	if (debug_output && comm->me == 0) {
		fprintf(stderr, "  ** -->   %s is now (%g, %g, %g)\n", name,
		        arr[idx][0],arr[idx][1],arr[idx][2]);
		fprintf(stderr, "  ** --> d_%s is now (%g, %g, %g)\n", name,
		        k_arr.d_view(idx,0),k_arr.d_view(idx,1),k_arr.d_view(idx,2));
		fprintf(stderr, "  ** --> h_%s is now (%g, %g, %g)\n", name,
		        k_arr.h_view(idx,0),k_arr.h_view(idx,1),k_arr.h_view(idx,2));
	}
}


template <class DeviceType>
void FixRigidKokkos<DeviceType>::init()
{
  FixRigid::init();
  
  atomKK->k_mass.modify<LMPHostType>();
  atomKK->k_mass.sync<DeviceType>();

  
  debug_print(k_xcm, xcm, "xcm");
  debug_print(k_vcm, vcm, "vcm");
  debug_print(k_fcm, fcm, "fcm");
  debug_print(k_torque, torque, "torque");

  debug_print(k_fflag, fflag, "fflag");
  debug_print(k_tflag, tflag, "tflag");

  debug_print(k_omega, omega, "omega");
  debug_print(k_angmom, angmom, "angmom");

  debug_print(k_quat, quat, "quat");
  debug_print(k_inertia, inertia, "inertia");

  debug_print(k_ex_space, ex_space, "ex_space");
  debug_print(k_ey_space, ey_space, "ey_space");
  debug_print(k_ez_space, ez_space, "ez_space");
  
  
  
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

  debug_print(k_xcm, xcm, "xcm");
  debug_print(k_vcm, vcm, "vcm");
  debug_print(k_fcm, fcm, "fcm");
  debug_print(k_torque, torque, "torque");

  debug_print(k_fflag, fflag, "fflag");
  debug_print(k_tflag, tflag, "tflag");

  debug_print(k_omega, omega, "omega");
  debug_print(k_angmom, angmom, "angmom");

  debug_print(k_quat, quat, "quat");
  debug_print(k_inertia, inertia, "inertia");

  debug_print(k_ex_space, ex_space, "ex_space");
  debug_print(k_ey_space, ey_space, "ey_space");
  debug_print(k_ez_space, ez_space, "ez_space");
  
  if (debug_output && comm->me == 0) fprintf(stderr,"\n");


  // Grab all arrays you need for initial_integrate:
  double dtfm;
  
  {
    // Local block for Kokkos parallel for:

    // These are local arrays?
    auto l_masstotal = masstotal;
    auto l_ex_space = k_ex_space.d_view;
    auto l_ey_space = k_ey_space.d_view;
    auto l_ez_space = k_ez_space.d_view;

    auto l_torque = k_torque.d_view;
    auto l_angmom = k_angmom.d_view;
    auto l_omega  = k_omega.d_view;

    // These are handled by FixRigid itself:
    auto l_fflag = k_fflag.d_view;
    auto l_tflag = k_tflag.d_view;
    
    auto l_xcm = k_xcm.d_view;
    auto l_vcm = k_vcm.d_view;
    auto l_fcm = k_fcm.d_view;
    auto l_quat = k_quat.d_view;
    auto l_inertia = k_inertia.d_view;
    

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

      
      MathExtraKokkos::angmom_to_omega<DeviceType>(l_angmom,
                                                   l_ex_space,
                                                   l_ey_space,
                                                   l_ez_space,
                                                   l_inertia, l_omega, ibody);

      MathExtraKokkos::richardson<DeviceType>(l_quat, l_angmom, l_omega,
                                              l_inertia, dtq, ibody);
      auto q_ibody = Kokkos::subview(l_quat, ibody, Kokkos::ALL);
      double q0 = l_quat(ibody,0);
      double q1 = l_quat(ibody,1);
      double q2 = l_quat(ibody,2);
      double q3 = l_quat(ibody,3);

      
      l_ex_space(ibody,0) = q0*q0 + q1*q1 - q2*q2 - q3*q3;
      l_ex_space(ibody,1) = 2.0 * (q1*q2 + q0*q3);
      l_ex_space(ibody,2) = 2.0 * (q1*q3 - q0*q2);

      l_ey_space(ibody,0) = 2.0 * (q1*q2 - q0*q3);
      l_ey_space(ibody,1) = q0*q0 - q1*q1 + q2*q2 - q3*q3;
      l_ey_space(ibody,2) = 2.0 * (q2*q3 + q0*q1);

      l_ez_space(ibody,0) = 2.0 * (q1*q3 + q0*q2);
      l_ez_space(ibody,1) = 2.0 * (q2*q3 - q0*q1);
      l_ez_space(ibody,2) = q0*q0 - q1*q1 - q2*q2 + q3*q3;

    });
  } // Ends local block for Kokkos parallel lambda.
  // At this point, we need to set up the virial if required:
  if (vflag) v_setup(vflag);
  else evflag = 0;

  
  // Convert in-body coordinates and stuff back to per-atom quantities:
  //set_xv_kokkos(); // THIS IS A THUNK
  
  atomKK->sync(execution_space, datamask_read);
  atomKK->modified(execution_space, datamask_modify);
  set_xv();

}


template <class DeviceType>
void FixRigidKokkos<DeviceType>::set_xv_kokkos()
{

}


template <class DeviceType>
void FixRigidKokkos<DeviceType>::set_v_kokkos()
{

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
