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
#include "domain_kokkos.h"
#include "memory_kokkos.h"



using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace RigidConst;

constexpr const bool debug_output = true;

/* ----------------------------------------------------------------------
   Contributing author: Stefan Paquay (Brandeis U, stefanpaquay@gmail.com)
------------------------------------------------------------------------- */


/* ---------------------------------------------------------------------- */

// // Some MathExtra functions had to be ported to Kokkos. In the future we
// // probably want to collect all of them somewhere.

// namespace MathExtraKokkos {

// // These trigger warning about redundant inline on explicit specialization :(
// KOKKOS_INLINE_FUNCTION
// template <class DeviceType>
// void angmom_to_omega(typename ArrayTypes<DeviceType>::t_v_array m,
//                      typename ArrayTypes<DeviceType>::t_x_array ex,
//                      typename ArrayTypes<DeviceType>::t_x_array ey,
//                      typename ArrayTypes<DeviceType>::t_x_array ez,
//                      typename ArrayTypes<DeviceType>::t_x_array idiag,
//                      typename ArrayTypes<DeviceType>::t_v_array w, int ibody)
// {
//   double wbody[3];

//   if (idiag(ibody,0) == 0.0) wbody[0] = 0.0;
//   else wbody[0] = (m(ibody,0)*ex(ibody,0) + m(ibody,1)*ex(ibody,1) + m(ibody,2)*ex(ibody,2)) / idiag(ibody,0);
//   if (idiag(ibody,1) == 0.0) wbody[1] = 0.0;
//   else wbody[1] = (m(ibody,0)*ey(ibody,0) + m(ibody,1)*ey(ibody,1) + m(ibody,2)*ey(ibody,2)) / idiag(ibody,1);
//   if (idiag(ibody,2) == 0.0) wbody[2] = 0.0;
//   else wbody[2] = (m(ibody,0)*ez(ibody,0) + m(ibody,1)*ez(ibody,1) + m(ibody,2)*ez(ibody,2)) / idiag(ibody,2);

//   w(ibody,0) = wbody[0]*ex(ibody,0) + wbody[1]*ey(ibody,0) + wbody[2]*ez(ibody,0);
//   w(ibody,1) = wbody[0]*ex(ibody,1) + wbody[1]*ey(ibody,1) + wbody[2]*ez(ibody,1);
//   w(ibody,2) = wbody[0]*ex(ibody,2) + wbody[1]*ey(ibody,2) + wbody[2]*ez(ibody,2);
// }

// KOKKOS_INLINE_FUNCTION
// template <typename a_arr_type, typename b_arr_type>
// inline void vecquat(a_arr_type a, b_arr_type b, double c[4])
// {
//   c[0] = -a[0]*b[1] - a[1]*b[2] - a[2]*b[3];
//   c[1] = b[0]*a[0] + a[1]*b[3] - a[2]*b[2];
//   c[2] = b[0]*a[1] + a[2]*b[1] - a[0]*b[3];
//   c[3] = b[0]*a[2] + a[0]*b[2] - a[1]*b[1];
// }

// /* ----------------------------------------------------------------------
//    matrix times vector
// ------------------------------------------------------------------------- */
// KOKKOS_INLINE_FUNCTION
// template <typename arr_out_type>
// inline void matvec(const double m[3][3], const double v[3],
//                    arr_out_type ans)
// {
//   ans[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];
//   ans[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];
//   ans[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];
// }


// /* ----------------------------------------------------------------------
//    transposed matrix times vector
// ------------------------------------------------------------------------- */
// KOKKOS_INLINE_FUNCTION
// template <typename v_arr_type, typename arr_out_type>
// void transpose_matvec(const double m[3][3], v_arr_type v,
//                              arr_out_type ans)
// {
//   ans[0] = m[0][0]*v[0] + m[1][0]*v[1] + m[2][0]*v[2];
//   ans[1] = m[0][1]*v[0] + m[1][1]*v[1] + m[2][1]*v[2];
//   ans[2] = m[0][2]*v[0] + m[1][2]*v[1] + m[2][2]*v[2];
// }


//   /* ----------------------------------------------------------------------
//    matrix times vector
// ------------------------------------------------------------------------- */
// KOKKOS_INLINE_FUNCTION
// template <typename e_arr_type, typename v_arr_type, typename ans_arr_type>
// void matvec(e_arr_type ex, e_arr_type ey,
//                    e_arr_type ez, v_arr_type v, ans_arr_type ans)
// {
//   ans[0] = ex[0]*v[0] + ey[0]*v[1] + ez[0]*v[2];
//   ans[1] = ex[1]*v[0] + ey[1]*v[1] + ez[1]*v[2];
//   ans[2] = ex[2]*v[0] + ey[2]*v[1] + ez[2]*v[2];
// }


// /* ----------------------------------------------------------------------
//    normalize a quaternion
// ------------------------------------------------------------------------- */
// KOKKOS_INLINE_FUNCTION
// template <typename q_arr_type>
//  void qnormalize(q_arr_type q)
// {
//   double norm = 1.0 / sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
//   q[0] *= norm;
//   q[1] *= norm;
//   q[2] *= norm;
//   q[3] *= norm;
// }


// /* ----------------------------------------------------------------------
//    conjugate of a quaternion: qc = conjugate of q
//    assume q is of unit length
// ------------------------------------------------------------------------- */
// KOKKOS_INLINE_FUNCTION
//  void qconjugate(double q[4], double qc[4])
// {
//   qc[0] = q[0];
//   qc[1] = -q[1];
//   qc[2] = -q[2];
//   qc[3] = -q[3];
// }

// /* ----------------------------------------------------------------------
//    compute rotation matrix from quaternion
//    quat = [w i j k]
// ------------------------------------------------------------------------- */
// KOKKOS_INLINE_FUNCTION
// void quat_to_mat(const double quat[4], double mat[3][3])
// {
//   double w2 = quat[0]*quat[0];
//   double i2 = quat[1]*quat[1];
//   double j2 = quat[2]*quat[2];
//   double k2 = quat[3]*quat[3];
//   double twoij = 2.0*quat[1]*quat[2];
//   double twoik = 2.0*quat[1]*quat[3];
//   double twojk = 2.0*quat[2]*quat[3];
//   double twoiw = 2.0*quat[1]*quat[0];
//   double twojw = 2.0*quat[2]*quat[0];
//   double twokw = 2.0*quat[3]*quat[0];

//   mat[0][0] = w2+i2-j2-k2;
//   mat[0][1] = twoij-twokw;
//   mat[0][2] = twojw+twoik;

//   mat[1][0] = twoij+twokw;
//   mat[1][1] = w2-i2+j2-k2;
//   mat[1][2] = twojk-twoiw;

//   mat[2][0] = twoik-twojw;
//   mat[2][1] = twojk+twoiw;
//   mat[2][2] = w2-i2-j2+k2;
// }



// /* ----------------------------------------------------------------------
//    compute space-frame ex,ey,ez from current quaternion q
//    ex,ey,ez = space-frame coords of 1st,2nd,3rd principal axis
//    operation is ex = q' d q = Q d, where d is (1,0,0) = 1st axis in body frame
// ------------------------------------------------------------------------- */
// KOKKOS_INLINE_FUNCTION
// template <typename q_arr_type, typename e_arr_type>
//  void q_to_exyz(q_arr_type q,
//                       e_arr_type ex, e_arr_type ey, e_arr_type ez)
// {
//   double q0 = q[0];
//   double q1 = q[1];
//   double q2 = q[2];
//   double q3 = q[3];

//   ex[0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
//   ex[1] = 2.0 * (q1*q2 + q0*q3);
//   ex[2] = 2.0 * (q1*q3 - q0*q2);

//   ey[0] = 2.0 * (q1*q2 - q0*q3);
//   ey[1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
//   ey[2] = 2.0 * (q2*q3 + q0*q1);

//   ez[0] = 2.0 * (q1*q3 + q0*q2);
//   ez[1] = 2.0 * (q2*q3 - q0*q1);
//   ez[2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;
// }




// /* ----------------------------------------------------------------------
//    compute omega from angular momentum
//    w = omega = angular velocity in space frame
//    wbody = angular velocity in body frame
//    project space-frame angular momentum onto body axes
//      and divide by principal moments
// ------------------------------------------------------------------------- */

// KOKKOS_INLINE_FUNCTION
// template <typename v_arr_type, typename x_arr_type>
//  void mq_to_omega(v_arr_type m, double q[4],
//                         x_arr_type moments, v_arr_type w)
// {
//   double wbody[3];
//   double rot[3][3];

//   MathExtraKokkos::quat_to_mat(q,rot);
//   MathExtraKokkos::transpose_matvec(rot,m,wbody);
//   if (moments[0] == 0.0) wbody[0] = 0.0;
//   else wbody[0] /= moments[0];
//   if (moments[1] == 0.0) wbody[1] = 0.0;
//   else wbody[1] /= moments[1];
//   if (moments[2] == 0.0) wbody[2] = 0.0;
//   else wbody[2] /= moments[2];
//   MathExtraKokkos::matvec(rot,wbody,w);
// }


// KOKKOS_INLINE_FUNCTION
// template <typename q_arr_type, typename x_arr_type, typename v_arr_type>
// void richardson(q_arr_type q_ibody, v_arr_type m_ibody, v_arr_type w_ibody,
//                 x_arr_type moments_ibody, double dtq)
// {
//   // full update from dq/dt = 1/2 w q

//   double wq[4];
//   MathExtraKokkos::vecquat(w_ibody,q_ibody,wq);

//   double qfull[4];
//   qfull[0] = q_ibody[0] + dtq * wq[0];
//   qfull[1] = q_ibody[1] + dtq * wq[1];
//   qfull[2] = q_ibody[2] + dtq * wq[2];
//   qfull[3] = q_ibody[3] + dtq * wq[3];
//   MathExtraKokkos::qnormalize(qfull);

//   // 1st half update from dq/dt = 1/2 w q

//   double qhalf[4];
//   qhalf[0] = q_ibody[0] + 0.5*dtq * wq[0];
//   qhalf[1] = q_ibody[1] + 0.5*dtq * wq[1];
//   qhalf[2] = q_ibody[2] + 0.5*dtq * wq[2];
//   qhalf[3] = q_ibody[3] + 0.5*dtq * wq[3];
//   MathExtraKokkos::qnormalize(qhalf);


//   // re-compute omega at 1/2 step from m at 1/2 step and q at 1/2 step
//   // recompute wq

//   MathExtraKokkos::mq_to_omega(m_ibody,qhalf,moments_ibody,w_ibody);
//   MathExtraKokkos::vecquat(w_ibody,qhalf,wq);

//   // 2nd half update from dq/dt = 1/2 w q

//   qhalf[0] += 0.5*dtq * wq[0];
//   qhalf[1] += 0.5*dtq * wq[1];
//   qhalf[2] += 0.5*dtq * wq[2];
//   qhalf[3] += 0.5*dtq * wq[3];
//   MathExtraKokkos::qnormalize(qhalf);

//   // corrected Richardson update

//   q_ibody[0] = 2.0*qhalf[0] - qfull[0];
//   q_ibody[1] = 2.0*qhalf[1] - qfull[1];
//   q_ibody[2] = 2.0*qhalf[2] - qfull[2];
//   q_ibody[3] = 2.0*qhalf[3] - qfull[3];
//    MathExtraKokkos::qnormalize(q_ibody);
// }
//
//
// } // MathExtraKokkos



// Debug helper functions:
template <class DeviceType>
template <typename kokkos_arr, typename base_arr>
void FixRigidKokkos<DeviceType>::debug_print_vec(kokkos_arr k_arr, base_arr arr,
                                                 const char *name, int idx)
{
  if (debug_output && comm->me == 0) {
    fprintf(stderr, "  ** -->   %s is now (%g, %g, %g)\n", name,
            arr[idx][0],arr[idx][1],arr[idx][2]);
#ifndef KOKKOS_ENABLE_CUDA
    fprintf(stderr, "  ** --> d_%s is now (%g, %g, %g)\n", name,
            k_arr.d_view(idx,0),k_arr.d_view(idx,1),k_arr.d_view(idx,2));
#endif
    fprintf(stderr, "  ** --> h_%s is now (%g, %g, %g)\n", name,
            k_arr.h_view(idx,0),k_arr.h_view(idx,1),k_arr.h_view(idx,2));
  }
}



template <class DeviceType>
template <typename kokkos_arr, typename base_arr>
void FixRigidKokkos<DeviceType>::debug_print_quat(kokkos_arr k_arr, base_arr arr,
                                                  const char *name, int idx)
{
  if (debug_output && comm->me == 0) {
    fprintf(stderr, "  ** -->   %s is now (%g, %g, %g, %g)\n", name,
            arr[idx][0],arr[idx][1],arr[idx][2], arr[idx][3]);
#ifndef KOKKOS_ENABLE_CUDA
    fprintf(stderr, "  ** --> d_%s is now (%g, %g, %g, %g)\n", name,
            k_arr.d_view(idx,0),k_arr.d_view(idx,1),k_arr.d_view(idx,2),k_arr.d_view(idx,3));
#endif
    fprintf(stderr, "  ** --> h_%s is now (%g, %g, %g, %g)\n", name,
            k_arr.h_view(idx,0),k_arr.h_view(idx,1),k_arr.h_view(idx,2),k_arr.h_view(idx,3));
  }
}




// helper functions to "copy-and-swap" arrays. They function just like
// atomKK::create_kokkos, except the data in array is preserved and
// copied into data.

// 1d variant:
template <class DeviceType>
template <typename arr_type>
void FixRigidKokkos<DeviceType>::create_and_copy(arr_type &data,
                                                 typename arr_type::value_type *&array,
                                                 int n1, const char *name)
{
  typename arr_type::value_type *tmp = new typename arr_type::value_type[n1];
  if (!tmp) error->all(FLERR, "Memory allocation failed!");

  memcpy(tmp, array, n1*sizeof(typename arr_type::value_type));
  memoryKK->create_kokkos(data, array, n1, name);
  memcpy(array, tmp, n1*sizeof(typename arr_type::value_type));

  delete [] tmp;
}

// 2d variant:
template <typename DeviceType>
template <typename arr_type>
void FixRigidKokkos<DeviceType>::create_and_copy(arr_type &data,
                                                 typename arr_type::value_type **&array,
                                                 int n1, int n2,
                                                 const char *name)
{
  typename arr_type::value_type *tmp = new typename arr_type::value_type[n1*n2];
  if (!tmp) error->all(FLERR, "Memory allocation failed!");

  // This does assume the 2D array is fully contiguous...
  memcpy(tmp, array[0], n1*n2*sizeof(typename arr_type::value_type));
  memoryKK->create_kokkos(data, array, n1, n2, name);
  memcpy(array[0], tmp, n1*n2*sizeof(typename arr_type::value_type));

  delete [] tmp;
}


template <class DeviceType>
FixRigidKokkos<DeviceType>::FixRigidKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixRigid(lmp, narg, arg)
{
  kokkosable = 1;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  atomKK = (AtomKokkos *) atom;

  datamask_read   = (X_MASK | V_MASK | F_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK
                     | OMEGA_MASK | ANGMOM_MASK | TORQUE_MASK);
  datamask_modify = (X_MASK | V_MASK | OMEGA_MASK | ANGMOM_MASK);

  // Most arrays allocated in the constructor of FixRigid are either
  //   a) set after that, so we need not worry about preserving the data, or
  //   b) Set to all zeros which is also what Kokkos does by default.
  // Those that are allocated _and_set in the c-tor use create_and_copy
  // to preserve the values already in the array, the rest do not need to.
  //
  // nrigid, body, tflag and fflag are set to specific values in base c-tor.

  int nmax = atomKK->nmax;

  //HAT::t_int_1d h_body(body,atom->local);
  //HAT::t_int_1d h_nrigid(nrigid,nbody);
  //HAT::t_f_array h_tflag(tflag,nbody,3);
  //HAT::t_f_array h_fflag(fflag,nbody,3);
  /*
  k_body.h_view = h_body;
  k_body.modify<LMPHostType>();
  k_body.sync<DeviceType>();

  k_nrigid.h_view = h_nrigid;
  k_nrigid.modify<LMPHostType>();
  k_nrigid.sync<DeviceType>();

  k_tflag.h_view = h_tflag;
  k_tflag.modify<LMPHostType>();
  k_tflag.sync<DeviceType>();

  k_fflag.h_view = h_fflag;
  k_fflag.modify<LMPHostType>();
  k_fflag.sync<DeviceType>();

  */


  if (debug_output && comm->me == 0) {
    fprintf(stderr, "Before call to grow_arrays, k_body has contents:\n");
    for (int i = 0; i < nmax; ++i) {
      fprintf(stderr, "%d ", k_body.h_view(i));
    }
    fprintf(stderr, "\n");
  }

  // The call to grow_arrays has to be after the create_and_copys because
  // else the empty arrays will be overwritten and their contents lost?
  grow_arrays(nmax);
  if (debug_output && comm->me == 0) {
    fprintf(stderr, "After call to grow_arrays, k_body has contents:\n");
    for (int i = 0; i < nmax; ++i) {
      fprintf(stderr, "%d ", k_body.h_view(i));
    }
    fprintf(stderr, "\n");
  }

  memoryKK->create_kokkos(k_masstotal, masstotal, nbody, "rigid/kk:masstotal");
  memoryKK->create_kokkos(k_xcm,xcm,nbody,3,"rigid/kk:xcm");
  memoryKK->create_kokkos(k_vcm,vcm,nbody,3,"rigid/kk:vcm");
  memoryKK->create_kokkos(k_fcm,fcm,nbody,3,"rigid/kk:fcm");
  memoryKK->create_kokkos(k_torque,torque,nbody,3,"rigid/kk:torque");


  memoryKK->create_kokkos(k_omega, omega, nbody,3,"rigid/kk:omega");
  memoryKK->create_kokkos(k_angmom,angmom,nbody,3,"rigid/kk:angmom");

  memoryKK->create_kokkos(k_quat,   quat,    nbody, 4,"rigid/kk:quat");
  memoryKK->create_kokkos(k_inertia,inertia, nbody, 3,"rigid/kk:inertia");

  memoryKK->create_kokkos(k_ex_space,ex_space, nbody, 3,"rigid/kk:ex_space");
  memoryKK->create_kokkos(k_ey_space,ey_space, nbody, 3,"rigid/kk:ey_space");
  memoryKK->create_kokkos(k_ez_space,ez_space, nbody, 3,"rigid/kk:ez_space");

  memoryKK->create_kokkos(k_sum, sum, nbody, 6, "rigid/kk:sum");
  memoryKK->create_kokkos(k_all, all, nbody, 6, "rigid/kk:all");
  memoryKK->create_kokkos(k_langextra,langextra,nbody,6,"rigid/kk:langextra");

  memoryKK->create_kokkos(k_imagebody,imagebody,nbody,"rigid/kk:imagebody");
  memoryKK->create_kokkos(k_remapflag,remapflag,nbody,4,"rigid/kk:remapflag");


}

template <class DeviceType>
FixRigidKokkos<DeviceType>::~FixRigidKokkos()
{
  memoryKK->destroy_kokkos(k_nrigid,nrigid);
  memoryKK->destroy_kokkos(k_tflag,tflag);
  memoryKK->destroy_kokkos(k_fflag,fflag);
  memoryKK->destroy_kokkos(k_body,body);

  memoryKK->destroy_kokkos(k_masstotal,masstotal);
  memoryKK->destroy_kokkos(k_xcm,xcm);
  memoryKK->destroy_kokkos(k_vcm,vcm);
  memoryKK->destroy_kokkos(k_fcm,fcm);
  memoryKK->destroy_kokkos(k_torque,torque);


  memoryKK->destroy_kokkos(k_omega,omega);
  memoryKK->destroy_kokkos(k_angmom,angmom);

  memoryKK->destroy_kokkos(k_quat,quat);
  memoryKK->destroy_kokkos(k_inertia,inertia);

  memoryKK->destroy_kokkos(k_ex_space,ex_space);
  memoryKK->destroy_kokkos(k_ey_space,ey_space);
  memoryKK->destroy_kokkos(k_ez_space,ez_space);

  memoryKK->destroy_kokkos(k_sum,sum);
  memoryKK->destroy_kokkos(k_all,all);
  memoryKK->destroy_kokkos(k_langextra,langextra);

  memoryKK->destroy_kokkos(k_imagebody,imagebody);
  memoryKK->destroy_kokkos(k_remapflag,remapflag);


  memoryKK->destroy_kokkos(k_xcmimage, xcmimage);
  memoryKK->destroy_kokkos(k_displace, displace);
  memoryKK->destroy_kokkos(k_eflags,eflags);
  memoryKK->destroy_kokkos(k_orient,orient);
  memoryKK->destroy_kokkos(k_dorient,dorient);

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

  // The host code uses these in FixRigid::init():
  // They should be synched first.
  // tflag, fflag, body, mu, radius, rmass, mass, ellipsoid, line, tri, type,
  // nlocal, eflags, image, x, sum, xcmimage, inertia, {ex,ey,ez}_space, quat
  k_body.sync<LMPHostType>();
  k_quat.sync<LMPHostType>();
  k_inertia.sync<LMPHostType>();
  k_ex_space.sync<LMPHostType>();
  k_ey_space.sync<LMPHostType>();
  k_ez_space.sync<LMPHostType>();

  k_tflag.sync<LMPHostType>();
  k_fflag.sync<LMPHostType>();
  k_xcmimage.sync<LMPHostType>();
  k_all.sync<LMPHostType>();
  k_sum.sync<LMPHostType>();

  k_vcm.sync<LMPHostType>();
  k_angmom.sync<LMPHostType>();
  k_imagebody.sync<LMPHostType>();
  k_xcm.sync<LMPHostType>();
  k_xcmimage.sync<LMPHostType>();
  k_displace.sync<LMPHostType>();

  atomKK->k_image.sync<LMPHostType>();
  atomKK->k_x.sync<LMPHostType>();


  if (comm->me == 0) {
    fprintf(stderr, "  ** --> IN(start; kk): quat and displace\n");
    debug_print_quat(k_quat, quat, "quat");
    debug_print_vec(k_displace, displace, "displace");
  }

  // These are also modified:
  // eflags, inertia, quat, body

  FixRigid::init();

  atomKK->k_mass.modify<LMPHostType>();
  atomKK->k_mass.sync<DeviceType>();

  k_imagebody.modify<LMPHostType>();
  k_vcm.modify<LMPHostType>();
  k_angmom.modify<LMPHostType>();

  k_xcmimage.modify<LMPHostType>();
  k_xcm.modify<LMPHostType>();
  k_displace.modify<LMPHostType>();

  k_quat.modify<LMPHostType>();
  k_inertia.modify<LMPHostType>();
  k_ex_space.modify<LMPHostType>();
  k_ey_space.modify<LMPHostType>();
  k_ez_space.modify<LMPHostType>();

  k_all.modify<LMPHostType>();
  k_sum.modify<LMPHostType>();

  k_quat.sync<DeviceType>();
  k_displace.sync<DeviceType>();
  k_inertia.sync<DeviceType>();
  k_ex_space.sync<DeviceType>();
  k_ey_space.sync<DeviceType>();
  k_ez_space.sync<DeviceType>();

  k_all.sync<DeviceType>();
  k_sum.sync<DeviceType>();
  k_body.sync<DeviceType>();

  if (comm->me == 0) {
    fprintf(stderr, "  ** --> IN(exit; kk): quat and displace\n");
    debug_print_quat(k_quat, quat, "quat");
    debug_print_vec(k_displace, displace, "displace");
  }
}






template <class DeviceType>
void FixRigidKokkos<DeviceType>::setup(int vflag)
{
  // setup modifies the following:
  // sum, all, torque, langextra, omega
  //
  // setup uses the following:
  // sum, f, body, all, x, xcm, atom->torque, torque, langextra, omega,
  // e{x,y,z}_space, xcmimage, inertia
  //
  // also calls set_v, which modifies the following:
  // v, atom_vec->omega, atom_vec->angmom
  // uses the following:
  // e{x,y,z}_space, displace, omega, v, vcm, atom_vec->mass, xcmimage, x

  k_sum.sync<LMPHostType>();
  k_all.sync<LMPHostType>();
  k_torque.sync<LMPHostType>();
  k_langextra.sync<LMPHostType>();
  k_omega.sync<LMPHostType>();

  atomKK->k_x.sync<LMPHostType>();
  atomKK->k_v.sync<LMPHostType>();
  atomKK->k_f.sync<LMPHostType>();

  k_body.sync<LMPHostType>();
  k_xcm.sync<LMPHostType>();
  k_ex_space.sync<LMPHostType>();
  k_ey_space.sync<LMPHostType>();
  k_ez_space.sync<LMPHostType>();
  k_xcmimage.sync<LMPHostType>();
  k_inertia.sync<LMPHostType>();

  k_vcm.sync<LMPHostType>();
  k_quat.sync<LMPHostType>();

  // modifies:
  k_sum.modify<LMPHostType>();
  k_all.modify<LMPHostType>();
  k_torque.modify<LMPHostType>();
  k_langextra.modify<LMPHostType>();
  k_omega.modify<LMPHostType>();
  k_quat.modify<LMPHostType>();
  atomKK->k_v.modify<LMPHostType>();

  FixRigid::setup(vflag);

  k_sum.sync<DeviceType>();
  k_all.sync<DeviceType>();
  k_torque.sync<DeviceType>();
  k_langextra.sync<DeviceType>();
  k_omega.sync<DeviceType>();
  k_quat.sync<DeviceType>();
  atomKK->k_v.sync<DeviceType>();


  if (debug_output && comm->me == 0) {
    fprintf(stderr, "  ** --> SE (exit2): After synching, we have:\n");
    debug_print_vec(k_fcm, fcm, "fcm");
    debug_print_vec(k_omega, omega, "omega");
  }
}


template <class DeviceType>
void FixRigidKokkos<DeviceType>::pre_neighbor()
{
  // pre_neighbor modifies both xcm and imagebody
  // and xcmimage and body
  k_xcm.sync<LMPHostType>();
  k_body.sync<LMPHostType>();
  k_xcmimage.sync<LMPHostType>();

  FixRigid::pre_neighbor();

  k_xcm.modify<LMPHostType>();

}





template <class DeviceType>
void FixRigidKokkos<DeviceType>::initial_integrate(int vflag)
{
  // initial integrate uses the following:
  // vcm, fcm, fflag, angmom, torque, tflag, ex_space, ey_space,
  // ez_space, inertia, omega, quat, x, v, f, rmass, mass, xcmimage, body
  //
  // initial integrate modifies the following:
  // vcm, xcm, angmom, omega, quat, ex_space, ey_space, ez_space, x, v
  //
  // set_xv uses:
  // body, xcmimage, x, v, omega, vcm, xcm, mass, f, displace
  //
  // set_xv modifies:
  // x, v

  k_vcm.sync<LMPHostType>();
  k_xcm.sync<LMPHostType>();
  k_fcm.sync<LMPHostType>();

  k_fflag.sync<LMPHostType>();
  k_tflag.sync<LMPHostType>();
  k_body.sync<LMPHostType>();
  k_xcmimage.sync<LMPHostType>();
  k_displace.sync<LMPHostType>();

  k_angmom.sync<LMPHostType>();
  k_omega.sync<LMPHostType>();
  k_torque.sync<LMPHostType>();
  k_inertia.sync<LMPHostType>();

  k_quat.sync<LMPHostType>();

  k_ex_space.sync<LMPHostType>();
  k_ey_space.sync<LMPHostType>();
  k_ez_space.sync<LMPHostType>();
  k_xcmimage.sync<LMPHostType>();

  atomKK->k_x.sync<LMPHostType>();
  atomKK->k_v.sync<LMPHostType>();
  atomKK->k_f.sync<LMPHostType>();
  atomKK->k_mass.sync<LMPHostType>();
  atomKK->k_rmass.sync<LMPHostType>();

  FixRigid::initial_integrate(vflag);

  k_vcm.modify<LMPHostType>();
  k_xcm.modify<LMPHostType>();
  k_angmom.modify<LMPHostType>();
  k_omega.modify<LMPHostType>();
  k_quat.modify<LMPHostType>();

  k_ex_space.modify<LMPHostType>();
  k_ey_space.modify<LMPHostType>();
  k_ez_space.modify<LMPHostType>();

  atomKK->k_x.modify<LMPHostType>();
  atomKK->k_v.modify<LMPHostType>();


  k_vcm.sync<DeviceType>();
  k_xcm.sync<DeviceType>();
  k_angmom.sync<DeviceType>();
  k_omega.sync<DeviceType>();
  k_quat.sync<DeviceType>();

  k_ex_space.sync<DeviceType>();
  k_ey_space.sync<DeviceType>();
  k_ez_space.sync<DeviceType>();

  atomKK->k_x.sync<DeviceType>();
  atomKK->k_v.sync<DeviceType>();


}


template <class DeviceType>
void FixRigidKokkos<DeviceType>::grow_arrays(int nmax)
{
  // This mirrors FixRigid::grow_arrays.
  memoryKK->grow_kokkos(k_body,body, nmax,"rigid/kk:body");
  memoryKK->grow_kokkos(k_xcmimage, xcmimage,nmax,"rigid/kk:xcmimage");
  memoryKK->grow_kokkos(k_displace, displace,nmax,3,"rigid/kk:displace");
  if (extended) {
    memoryKK->grow_kokkos(k_eflags,eflags,nmax,"rigid/kk:eflags");
    if (orientflag) memoryKK->grow_kokkos(k_orient,orient, nmax, orientflag,"rigid/kk:orient");
    if (dorientflag) memoryKK->grow_kokkos(k_dorient,dorient,nmax,3,"rigid/kk:dorient");
  }

  // check for regrow of vatom
  // must be done whether per-atom virial is accumulated on this step or not
  //   b/c this is only time grow_array() may be called
  // need to regrow b/c vatom is calculated before and after atom migration

  if (nmax > maxvatom) {
    maxvatom = atomKK->nmax;
    memory->grow(vatom,maxvatom,6,"fix:vatom");
  }
}


template <class DeviceType>
void FixRigidKokkos<DeviceType>::set_xv_kokkos()
{
}

template <class DeviceType>
void FixRigidKokkos<DeviceType>::final_integrate()
{
  // final_integrate modifies (in [] only if extended)
  // vcm, angmom, omega, v, [atom_vec->angmom, atom_vec->omega]
  //
  // final_integrate uses:
  // vcm, fcm, fflag, angmom, torque, tflag, x, v, f,
  // omega, mass, xcmimage, body
  // [atom_vec->omega, atom_vec->ellipsoid, ebonus],

  k_vcm.sync<LMPHostType>();
  k_fcm.sync<LMPHostType>();
  k_fflag.sync<LMPHostType>();
  k_angmom.sync<LMPHostType>();
  k_torque.sync<LMPHostType>();
  k_tflag.sync<LMPHostType>();
  atomKK->k_x.sync<LMPHostType>();
  atomKK->k_v.sync<LMPHostType>();
  atomKK->k_f.sync<LMPHostType>();
  atomKK->k_mass.sync<LMPHostType>();
  atomKK->k_rmass.sync<LMPHostType>();

  k_omega.sync<LMPHostType>();
  k_xcmimage.sync<LMPHostType>();
  k_body.sync<LMPHostType>();

  k_vcm.modify<LMPHostType>();
  k_angmom.modify<LMPHostType>();
  k_omega.modify<LMPHostType>();
  atomKK->k_v.modify<LMPHostType>();

  FixRigid::final_integrate();

  k_vcm.sync<DeviceType>();
  k_angmom.sync<DeviceType>();
  k_omega.sync<DeviceType>();
  atomKK->k_v.sync<DeviceType>();

}


template <class DeviceType>
double FixRigidKokkos<DeviceType>::compute_scalar()
{
  k_tflag.sync<LMPHostType>();
  k_fflag.sync<LMPHostType>();
  k_inertia.sync<LMPHostType>();

  k_angmom.sync<LMPHostType>();
  k_quat.sync<LMPHostType>();
  k_vcm.sync<LMPHostType>();


  return FixRigid::compute_scalar();
}


template <class DeviceType>
void FixRigidKokkos<DeviceType>::compute_forces_and_torques_kokkos()
{
}


template <class DeviceType>
void FixRigidKokkos<DeviceType>::set_v_kokkos()
{
}


template <class DeviceType>
void FixRigidKokkos<DeviceType>::post_force(int vflag)
{
  FixRigid::post_force(vflag);
}


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
