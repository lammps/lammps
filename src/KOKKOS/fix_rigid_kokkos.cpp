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
#include "kokkos_few.h"



using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace RigidConst;

constexpr const bool debug_output = false;

/* ----------------------------------------------------------------------
   Contributing author: Stefan Paquay (Brandeis U, stefanpaquay@gmail.com)
------------------------------------------------------------------------- */


/* ---------------------------------------------------------------------- */

// Some MathExtra functions had to be ported to Kokkos. In the future we
// probably want to collect all of them somewhere.

namespace MathExtraKokkos {

KOKKOS_INLINE_FUNCTION
void angmom_to_omega(DAT::t_v_array m,
                     DAT::t_x_array ex,
                     DAT::t_x_array ey,
                     DAT::t_x_array ez,
                     DAT::t_x_array idiag,
                     DAT::t_v_array w, int ibody)
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
KOKKOS_INLINE_FUNCTION
void vecquat(a_arr_type a, b_arr_type b, double c[4])
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
KOKKOS_INLINE_FUNCTION
void matvec(const double m[3][3], const double v[3],
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
KOKKOS_INLINE_FUNCTION
void transpose_matvec(const double m[3][3], v_arr_type v,
                      arr_out_type ans)
{
  ans[0] = m[0][0]*v[0] + m[1][0]*v[1] + m[2][0]*v[2];
  ans[1] = m[0][1]*v[0] + m[1][1]*v[1] + m[2][1]*v[2];
  ans[2] = m[0][2]*v[0] + m[1][2]*v[1] + m[2][2]*v[2];
}


  /* ----------------------------------------------------------------------
   matrix times vector
------------------------------------------------------------------------- */
template <typename e_arr_type, typename v_arr_type>
KOKKOS_INLINE_FUNCTION
Few<double,3> matvec(e_arr_type ex, e_arr_type ey,
                     e_arr_type ez, v_arr_type v)
{
  Few<double,3> ans = { ex[0]*v[0] + ey[0]*v[1] + ez[0]*v[2],
                        ex[1]*v[0] + ey[1]*v[1] + ez[1]*v[2],
                        ex[2]*v[0] + ey[2]*v[1] + ez[2]*v[2] };
  return ans;
}


/* ----------------------------------------------------------------------
   normalize a quaternion
------------------------------------------------------------------------- */
template <typename q_arr_type>
KOKKOS_INLINE_FUNCTION
 void qnormalize(q_arr_type q)
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
KOKKOS_INLINE_FUNCTION
void qconjugate(double q[4], double qc[4])
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
KOKKOS_INLINE_FUNCTION
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
template <typename q_arr_type, typename e_arr_type>
KOKKOS_INLINE_FUNCTION
void q_to_exyz(q_arr_type q,
               e_arr_type ex, e_arr_type ey, e_arr_type ez)
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
KOKKOS_INLINE_FUNCTION
void mq_to_omega(v_arr_type m, double q[4],
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


template <typename q_arr_type, typename x_arr_type, typename v_arr_type>
KOKKOS_INLINE_FUNCTION
void richardson(q_arr_type q_ibody, v_arr_type m_ibody, v_arr_type w_ibody,
                x_arr_type moments_ibody, double dtq)
{
  // full update from dq/dt = 1/2 w q

  double wq[4];
  MathExtraKokkos::vecquat(w_ibody,q_ibody,wq);

  double qfull[4];
  qfull[0] = q_ibody[0] + dtq * wq[0];
  qfull[1] = q_ibody[1] + dtq * wq[1];
  qfull[2] = q_ibody[2] + dtq * wq[2];
  qfull[3] = q_ibody[3] + dtq * wq[3];
  MathExtraKokkos::qnormalize(qfull);

  // 1st half update from dq/dt = 1/2 w q

  double qhalf[4];
  qhalf[0] = q_ibody[0] + 0.5*dtq * wq[0];
  qhalf[1] = q_ibody[1] + 0.5*dtq * wq[1];
  qhalf[2] = q_ibody[2] + 0.5*dtq * wq[2];
  qhalf[3] = q_ibody[3] + 0.5*dtq * wq[3];
  MathExtraKokkos::qnormalize(qhalf);


  // re-compute omega at 1/2 step from m at 1/2 step and q at 1/2 step
  // recompute wq

  MathExtraKokkos::mq_to_omega(m_ibody,qhalf,moments_ibody,w_ibody);
  MathExtraKokkos::vecquat(w_ibody,qhalf,wq);

  // 2nd half update from dq/dt = 1/2 w q

  qhalf[0] += 0.5*dtq * wq[0];
  qhalf[1] += 0.5*dtq * wq[1];
  qhalf[2] += 0.5*dtq * wq[2];
  qhalf[3] += 0.5*dtq * wq[3];
  MathExtraKokkos::qnormalize(qhalf);

  // corrected Richardson update

  q_ibody[0] = 2.0*qhalf[0] - qfull[0];
  q_ibody[1] = 2.0*qhalf[1] - qfull[1];
  q_ibody[2] = 2.0*qhalf[2] - qfull[2];
  q_ibody[3] = 2.0*qhalf[3] - qfull[3];
  MathExtraKokkos::qnormalize(q_ibody);
}


} // MathExtraKokkos





template <class DeviceType>
FixRigidKokkos<DeviceType>::FixRigidKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixRigid(lmp, narg, arg), rand_pool(comm->me + seed)
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
  // Those that are allocated _and_set in the c-tor use create_mirror_view
  // to preserve the values already in the array, the rest do not need to.
  //
  // nrigid, body, tflag and fflag are set to specific values in base c-tor.

  int nmax = atomKK->nmax;


  // create_mirror_view creates a host_view from device, not the
  // other way around! So we need this ugly manual memcpy.

  int *body_buffer = new int[nmax];
  double *nrigid_buffer = new double[nbody];
  double *tflag_buffer = new double[3*nbody];
  double *fflag_buffer = new double[3*nbody];

  memcpy(body_buffer,   body,   nmax*sizeof(int));
  memcpy(nrigid_buffer, nrigid, nbody*sizeof(double));
  memcpy(tflag_buffer,  tflag[0],  3*nbody*sizeof(double));
  memcpy(fflag_buffer,  fflag[0],  3*nbody*sizeof(double));

  memoryKK->create_kokkos(k_body,   body,   nmax, "rigid/kk:body");
  memoryKK->create_kokkos(k_nrigid, nrigid, nbody, "rigid/kk:nrigid");
  memoryKK->create_kokkos(k_tflag,  tflag,  nbody, 3, "rigid/kk:tflag");
  memoryKK->create_kokkos(k_fflag,  fflag,  nbody, 3, "rigid/kk:fflag");

  // The call to grow_arrays has to be after the create_mirror_views because
  // else the empty device arrays will be overwritten and the contents of the
  // host will be lost.
  grow_arrays(nmax);

  memcpy(body,   body_buffer,   nmax*sizeof(int));
  memcpy(nrigid, nrigid_buffer, nbody*sizeof(double));
  memcpy(tflag[0],  tflag_buffer,  3*nbody*sizeof(double));
  memcpy(fflag[0],  fflag_buffer,  3*nbody*sizeof(double));

  k_body.modify<LMPHostType>();
  k_nrigid.modify<LMPHostType>();
  k_tflag.modify<LMPHostType>();
  k_fflag.modify<LMPHostType>();

  k_body.sync<DeviceType>();
  k_nrigid.sync<DeviceType>();
  k_tflag.sync<DeviceType>();
  k_fflag.sync<DeviceType>();

  delete [] body_buffer;
  delete [] nrigid_buffer;
  delete [] tflag_buffer;
  delete [] fflag_buffer;


  if (debug_output && comm->me == 0) {
    fprintf(screen, "Body contains:\n");
    for (int i = 0; i < nmax; ++i) {
      fprintf(screen, "%d ", body[i]);
    }
    fprintf(screen, "\n");
  }

  /*
  HAT::t_int_1d h_body(body, nmax);
  k_body.h_view = h_body;
  k_body.modify<LMPHostType>();
  k_body.d_view = create_mirror_view(DeviceType(), k_body.h_view);

  HAT::t_f_array h_tflag(*tflag, nbody, 3);
  k_tflag.h_view = h_tflag;
  k_tflag.modify<LMPHostType>();
  k_tflag.d_view = create_mirror_view(DeviceType(), k_tflag.h_view);

  HAT::t_f_array h_fflag(*fflag, nbody, 3);
  k_fflag.h_view = h_fflag;
  k_fflag.modify<LMPHostType>();
  k_fflag.d_view = create_mirror_view(DeviceType(), k_fflag.h_view);

  HAT::t_int_1d h_nrigid(nrigid, nbody);
  k_nrigid.h_view = h_nrigid;
  k_nrigid.modify<LMPHostType>();
  k_nrigid.d_view = create_mirror_view(DeviceType(), k_nrigid.h_view);

  grow_arrays(nmax);

  */

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

  if (debug_output && comm->me == 0) {
    fprintf(screen, "  ** --> Post init, ex_space = (%g, %g, %g)\n",
            k_ex_space.h_view(0,0), k_ex_space.h_view(0,1), k_ex_space.h_view(0,2));
    fprintf(screen, "  ** --> Post init, ey_space = (%g, %g, %g)\n",
            k_ey_space.h_view(0,0), k_ey_space.h_view(0,1), k_ey_space.h_view(0,2));
    fprintf(screen, "  ** --> Post init, ez_space = (%g, %g, %g)\n",
            k_ez_space.h_view(0,0), k_ez_space.h_view(0,1), k_ez_space.h_view(0,2));
  }

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

}


template <class DeviceType>
void FixRigidKokkos<DeviceType>::pre_neighbor()
{
  // pre_neighbor modifies both xcm and imagebody
  // image_shift modifies xcmimage and body
  k_xcm.sync<DeviceType>();
  k_imagebody.sync<DeviceType>();

  {
	  // Local block for Kokkos lambda.
    auto l_xcm = k_xcm.d_view;
    auto l_imagebody = k_imagebody.d_view;
    auto domainKK = static_cast<DomainKokkos *>(domain);

    int periodic_bits = 0;
    periodic_bits +=   domain->xperiodic;
    periodic_bits += 2*domain->yperiodic;
    periodic_bits += 4*domain->zperiodic;


    Few<double,3> prd{domain->prd[0], domain->prd[1], domain->prd[2]};
    Few<double,3> prd_lambda{domain->prd_lamda[0], domain->prd_lamda[1], domain->prd_lamda[2]};

    Few<double,6> h{domain->h[0], domain->h[1], domain->h[2],
                    domain->h[3], domain->h[4], domain->h[5]};
    Few<double,6> h_inv{domain->h_inv[0], domain->h_inv[1], domain->h_inv[2],
                        domain->h_inv[3], domain->h_inv[4], domain->h_inv[5]};


    Few<double,3> boxlo{domain->boxlo[0], domain->boxlo[1], domain->boxlo[2]};
    Few<double,3> boxhi{domain->boxhi[0], domain->boxhi[1], domain->boxhi[2]};

    Few<double,3> boxlo_lambda{domain->boxlo_lamda[0], domain->boxlo_lamda[1], domain->boxlo_lamda[2]};
    Few<double,3> boxhi_lambda{domain->boxhi_lamda[0], domain->boxhi_lamda[1], domain->boxhi_lamda[2]};

    Kokkos::parallel_for(nbody, LAMMPS_LAMBDA(const int &ibody){
      Few<double,3> xcm_ibody{l_xcm(ibody,0), l_xcm(ibody,1), l_xcm(ibody,2)};
      imageint imagebody_ibody = l_imagebody(ibody);

      auto new_xcm = domainKK->remap(prd, h, h_inv, triclinic, boxlo, boxhi,
                                     boxlo_lambda, boxhi_lambda, prd_lambda,
                                     xcm_ibody, imagebody_ibody, periodic_bits);

      l_imagebody(ibody) = imagebody_ibody;
      l_xcm(ibody,0) = new_xcm[0];
      l_xcm(ibody,1) = new_xcm[1];
      l_xcm(ibody,2) = new_xcm[2];
    });
  }

  k_xcm.modify<DeviceType>();
  k_imagebody.modify<DeviceType>();

  image_shift_kokkos();

}


template <class DeviceType>
void FixRigidKokkos<DeviceType>::image_shift_kokkos()
{
  int nlocal = atomKK->nlocal;

  k_body.sync<DeviceType>();
  k_xcmimage.sync<DeviceType>();

  {
    // Local block for Kokkos::parallel_for

    imageint tdim,bdim,xdim[3];
    auto l_image = atomKK->k_image.d_view;
    auto l_imagebody = k_imagebody.d_view;
    auto l_body = k_body.d_view;
    auto l_xcmimage = k_xcmimage.d_view;

    Kokkos::parallel_for(nlocal, LAMMPS_LAMBDA(const int &i) {
        if (l_body[i] < 0) return;
        int ibody = body[i];
        imageint xdim[3];
        imageint tdim = l_image[i] & IMGMASK;
        imageint bdim = l_imagebody[ibody] & IMGMASK;

        xdim[0] = IMGMAX + tdim - bdim;
        tdim = (l_image[i] >> IMGBITS) & IMGMASK;
        bdim = (l_imagebody[ibody] >> IMGBITS) & IMGMASK;
        xdim[1] = IMGMAX + tdim - bdim;
        tdim = l_image[i] >> IMG2BITS;
        bdim = l_imagebody[ibody] >> IMG2BITS;
        xdim[2] = IMGMAX + tdim - bdim;

        l_xcmimage[i] = (xdim[2] << IMG2BITS) | (xdim[1] << IMGBITS) | xdim[0];
    });
  }
  k_xcmimage.modify<DeviceType>();
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

  k_vcm.sync<DeviceType>();
  k_xcm.sync<DeviceType>();
  k_fcm.sync<DeviceType>();

  k_fflag.sync<DeviceType>();
  k_tflag.sync<DeviceType>();
  k_body.sync<DeviceType>();
  k_xcmimage.sync<DeviceType>();
  k_displace.sync<DeviceType>();

  k_angmom.sync<DeviceType>();
  k_omega.sync<DeviceType>();
  k_torque.sync<DeviceType>();
  k_inertia.sync<DeviceType>();

  k_quat.sync<DeviceType>();

  k_ex_space.sync<DeviceType>();
  k_ey_space.sync<DeviceType>();
  k_ez_space.sync<DeviceType>();
  k_xcmimage.sync<DeviceType>();
  k_masstotal.sync<DeviceType>();

  atomKK->k_x.sync<DeviceType>();
  atomKK->k_v.sync<DeviceType>();
  atomKK->k_f.sync<DeviceType>();
  atomKK->k_mass.sync<DeviceType>();
  atomKK->k_rmass.sync<DeviceType>();

  atomKK->sync(execution_space, datamask_read);

  // Grab all arrays you need for initial_integrate:
  double dtfm;

  {
    // Local block for Kokkos parallel for:

    auto l_masstotal = k_masstotal.d_view;
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

      MathExtraKokkos::angmom_to_omega(l_angmom,
                                       l_ex_space,
                                       l_ey_space,
                                       l_ez_space,
                                       l_inertia, l_omega, ibody);

      MathExtraKokkos::richardson(q_ibody, angmom_ibody,
                                  omega_ibody, inertia_ibody, dtq);

      MathExtraKokkos::q_to_exyz(q_ibody, ex_ibody, ey_ibody, ez_ibody);
    });

  } // Ends local block for Kokkos parallel lambda.
  // At this point, we need to set up the virial if required:
  if (vflag) v_setup(vflag);
  else evflag = 0;

  if (debug_output && comm->me == 0) {
    atomKK->k_x.sync<LMPHostType>();
    atomKK->k_v.sync<LMPHostType>();
    atomKK->k_type.sync<LMPHostType>();
    k_omega.sync<LMPHostType>();
    k_ex_space.sync<LMPHostType>();
    k_ey_space.sync<LMPHostType>();
    k_ez_space.sync<LMPHostType>();
    k_displace.sync<LMPHostType>();
    k_vcm.sync<LMPHostType>();
    k_xcm.sync<LMPHostType>();

    k_angmom.sync<LMPHostType>();
    int i;
    for (i = 0; i < atomKK->nlocal; ++i) {
      if (atomKK->k_tag.h_view(i) == 4) break;
    }
    int ibody = k_body.h_view(i);
    if (debug_output && comm->me == 0) {
      fprintf(screen, "\ni = %d, ibody = %d\n", i, ibody);
      fprintf(screen, "before set_xv().\n");
      fprintf(screen, "vcm is (%g, %g, %g)\n",
              k_vcm.h_view(ibody,0),
              k_vcm.h_view(ibody,1),
              k_vcm.h_view(ibody,2));
      fprintf(screen, "xcm is (%g, %g, %g)\n",
              k_xcm.h_view(ibody,0),
              k_xcm.h_view(ibody,1),
              k_xcm.h_view(ibody,2));
      fprintf(screen, "angmom is (%g, %g, %g)\n",
              k_angmom.h_view(ibody,0),
              k_angmom.h_view(ibody,1),
              k_angmom.h_view(ibody,2));
      fprintf(screen, "omega before is (%g, %g, %g)\n",
              k_omega.h_view(ibody,0),
              k_omega.h_view(ibody,1),
              k_omega.h_view(ibody,2));
      fprintf(screen, "=========================\n");

      fprintf(screen, "x before is (%g, %g, %g)\n",
              atomKK->k_x.h_view(i,0),
              atomKK->k_x.h_view(i,1),
              atomKK->k_x.h_view(i,2));
      fprintf(screen, "v before is (%g, %g, %g)\n",
              atomKK->k_v.h_view(i,0),
              atomKK->k_v.h_view(i,1),
              atomKK->k_v.h_view(i,2));
      fprintf(screen, "ex_space before is (%g, %g, %g)\n",
              k_ex_space.h_view(ibody,0),
              k_ex_space.h_view(ibody,1),
              k_ex_space.h_view(ibody,2));
      fprintf(screen, "ey_space before is (%g, %g, %g)\n",
              k_ey_space.h_view(ibody,0),
              k_ey_space.h_view(ibody,1),
              k_ey_space.h_view(ibody,2));
      fprintf(screen, "ez_space before is (%g, %g, %g)\n",
              k_ez_space.h_view(ibody,0),
              k_ez_space.h_view(ibody,1),
              k_ez_space.h_view(ibody,2));
      fprintf(screen, "displace before is (%g, %g, %g)\n",
              k_displace.h_view(ibody,0),
              k_displace.h_view(ibody,1),
              k_displace.h_view(ibody,2));
    }
  }

  // Convert in-body coordinates and etc. back to per-atom quantities:
  set_xv_kokkos();

  k_vcm.modify<DeviceType>();
  k_xcm.modify<DeviceType>();
  k_angmom.modify<DeviceType>();
  k_omega.modify<DeviceType>();
  k_quat.modify<DeviceType>();

  k_ex_space.modify<DeviceType>();
  k_ey_space.modify<DeviceType>();
  k_ez_space.modify<DeviceType>();

  atomKK->k_x.modify<DeviceType>();
  atomKK->k_v.modify<DeviceType>();

  if (debug_output && comm->me == 0) {
    atomKK->k_x.sync<LMPHostType>();
    atomKK->k_v.sync<LMPHostType>();
    k_omega.sync<LMPHostType>();
    k_ex_space.sync<LMPHostType>();
    k_ex_space.sync<LMPHostType>();
    k_ex_space.sync<LMPHostType>();
    k_displace.sync<LMPHostType>();
    atomKK->k_tag.sync<LMPHostType>();

    int i;
    for (i = 0; i < atomKK->nlocal; ++i) {
      if (atomKK->k_tag.h_view(i) == 4) break;
    }
    int ibody = k_body.h_view(i);

    fprintf(screen, "End of initial integrate:\n");
    fprintf(screen, "\ni = %d, ibody = %d\n", i, ibody);
    fprintf(screen, "vcm is (%g, %g, %g)\n",
            k_vcm.h_view(ibody,0),
            k_vcm.h_view(ibody,1),
            k_vcm.h_view(ibody,2));
    fprintf(screen, "xcm is (%g, %g, %g)\n",
            k_xcm.h_view(ibody,0),
            k_xcm.h_view(ibody,1),
            k_xcm.h_view(ibody,2));
    fprintf(screen, "angmom is (%g, %g, %g)\n",
            k_angmom.h_view(ibody,0),
            k_angmom.h_view(ibody,1),
            k_angmom.h_view(ibody,2));
    fprintf(screen, "omega  is (%g, %g, %g)\n",
            k_omega.h_view(ibody,0),
            k_omega.h_view(ibody,1),
            k_omega.h_view(ibody,2));
    fprintf(screen, "=========================\n");

    fprintf(screen, "x  is (%g, %g, %g)\n",
            atomKK->k_x.h_view(i,0),
            atomKK->k_x.h_view(i,1),
            atomKK->k_x.h_view(i,2));
    fprintf(screen, "v  is (%g, %g, %g)\n",
            atomKK->k_v.h_view(i,0),
            atomKK->k_v.h_view(i,1),
            atomKK->k_v.h_view(i,2));
    fprintf(screen, "ex_space  is (%g, %g, %g)\n",
            k_ex_space.h_view(ibody,0),
            k_ex_space.h_view(ibody,1),
            k_ex_space.h_view(ibody,2));
    fprintf(screen, "ey_space  is (%g, %g, %g)\n",
            k_ey_space.h_view(ibody,0),
            k_ey_space.h_view(ibody,1),
            k_ey_space.h_view(ibody,2));
    fprintf(screen, "ez_space  is (%g, %g, %g)\n",
            k_ez_space.h_view(ibody,0),
            k_ez_space.h_view(ibody,1),
            k_ez_space.h_view(ibody,2));
    fprintf(screen, "displace  is (%g, %g, %g)\n",
            k_displace.h_view(ibody,0),
            k_displace.h_view(ibody,1),
            k_displace.h_view(ibody,2));

    fprintf(screen, "Done in initial_integrate!\n");
  }
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
    k_eflags.modify<LMPHostType>();
    if (orientflag) {
      memoryKK->grow_kokkos(k_orient,orient, nmax, orientflag,"rigid/kk:orient");
      k_orient.modify<LMPHostType>();
    }
    if (dorientflag) {
      memoryKK->grow_kokkos(k_dorient,dorient,nmax,3,"rigid/kk:dorient");
      k_dorient.modify<LMPHostType>();
    }
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
  double xy,xz,yz;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  int nlocal = atomKK->nlocal;

  if (triclinic) {
    xy = domain->xy;
    xz = domain->xz;
    yz = domain->yz;
  }

  // set x and v of each atom
  {
    // Local block for Kokkos parallel for:

    auto l_x = atomKK->k_x.d_view;
    auto l_v = atomKK->k_v.d_view;
    auto l_f = atomKK->k_f.d_view;

    auto l_type = atomKK->k_type.d_view;
    auto l_rmass = atomKK->k_rmass.d_view;
    auto l_mass  = atomKK->k_mass.d_view;

    auto l_ex_space = k_ex_space.d_view;
    auto l_ey_space = k_ey_space.d_view;
    auto l_ez_space = k_ez_space.d_view;

    auto l_xcm = k_xcm.d_view;
    auto l_vcm = k_vcm.d_view;
    auto l_omega = k_omega.d_view;
    auto l_displace = k_displace.d_view;

    auto l_xcmimage = k_xcmimage.d_view;
    auto l_body = k_body.d_view;

    Kokkos::parallel_for(nlocal, LAMMPS_LAMBDA(const int& i) {
      if (l_body[i] < 0) return;
      int ibody = l_body[i];

      double xbox = (l_xcmimage[i] & IMGMASK) - IMGMAX;
      double ybox = (l_xcmimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
      double zbox = (l_xcmimage[i] >> IMG2BITS) - IMGMAX;
      double x0, x1, x2;
      double v0, v1, v2;

      // save old positions and velocities for virial
      if (evflag) {
        if (triclinic == 0) {
          x0 = l_x(i,0) + xbox*xprd;
          x1 = l_x(i,1) + ybox*yprd;
          x2 = l_x(i,2) + zbox*zprd;
        } else {
          x0 = l_x(i,0) + xbox*xprd + ybox*xy + zbox*xz;
          x1 = l_x(i,1) + ybox*yprd + zbox*yz;
          x2 = l_x(i,2) + zbox*zprd;
        }
        v0 = l_v(i,0);
        v1 = l_v(i,1);
        v2 = l_v(i,2);
      }

      // x = displacement from center-of-mass, based on body orientation
      // v = vcm + omega around center-of-mass
      /*
      auto ex_space_ibody = Kokkos::subview(l_ex_space, ibody, Kokkos::ALL);
      auto ey_space_ibody = Kokkos::subview(l_ey_space, ibody, Kokkos::ALL);
      auto ez_space_ibody = Kokkos::subview(l_ez_space, ibody, Kokkos::ALL);

      auto omega_ibody = Kokkos::subview(l_omega, ibody, Kokkos::ALL);
      auto vcm_ibody   = Kokkos::subview(l_vcm,   ibody, Kokkos::ALL);
      auto xcm_ibody   = Kokkos::subview(l_xcm,   ibody, Kokkos::ALL);

      auto l_displace_i = Kokkos::subview(l_displace, i, Kokkos::ALL);

      auto ans = MathExtraKokkos::matvec(ex_space_ibody,ey_space_ibody,
                                         ez_space_ibody,l_displace_i);

      l_x(i,0) = ans[0];
      l_x(i,1) = ans[1];
      l_x(i,2) = ans[2];

      l_v(i,0) = l_omega(ibody,1)*l_x(i,2) - l_omega(ibody,2)*l_x(i,1) + l_vcm(ibody,0);
      l_v(i,1) = l_omega(ibody,2)*l_x(i,0) - l_omega(ibody,0)*l_x(i,2) + l_vcm(ibody,1);
      l_v(i,2) = l_omega(ibody,0)*l_x(i,1) - l_omega(ibody,1)*l_x(i,0) + l_vcm(ibody,2);
      */

      l_x(i,0) = l_ex_space(ibody,0)*l_displace(i,0) +
	      l_ey_space(ibody,0)*l_displace(i,1) +
	      l_ez_space(ibody,0)*l_displace(i,2);
      l_x(i,1) = l_ex_space(ibody,1)*l_displace(i,0) +
	      l_ey_space(ibody,1)*l_displace(i,1) +
	      l_ez_space(ibody,1)*l_displace(i,2);
      l_x(i,2) = l_ex_space(ibody,2)*l_displace(i,0) +
	      l_ey_space(ibody,2)*l_displace(i,1) +
	      l_ez_space(ibody,2)*l_displace(i,2);

      l_v(i,0) = l_omega(ibody,1)*l_x(i,2) - l_omega(ibody,2)*l_x(i,1) + l_vcm(ibody,0);
      l_v(i,1) = l_omega(ibody,2)*l_x(i,0) - l_omega(ibody,0)*l_x(i,2) + l_vcm(ibody,1);
      l_v(i,2) = l_omega(ibody,0)*l_x(i,1) - l_omega(ibody,1)*l_x(i,0) + l_vcm(ibody,2);

      // add center of mass to displacement
      // map back into periodic box via xbox,ybox,zbox
      // for triclinic, add in box tilt factors as well

      if (triclinic == 0) {
        l_x(i,0) += l_xcm(ibody,0) - xbox*xprd;
        l_x(i,1) += l_xcm(ibody,1) - ybox*yprd;
        l_x(i,2) += l_xcm(ibody,2) - zbox*zprd;
      } else {
        l_x(i,0) += l_xcm(ibody,0) - xbox*xprd - ybox*xy - zbox*xz;
        l_x(i,1) += l_xcm(ibody,1) - ybox*yprd - zbox*yz;
        l_x(i,2) += l_xcm(ibody,2) - zbox*zprd;
      }

      // virial = unwrapped coords dotted into body constraint force
      // body constraint force = implied force due to v change minus f external
      // assume f does not include forces internal to body
      // 1/2 factor b/c final_integrate contributes other half
      // assume per-atom contribution is due to constraint force on that atom

      if (evflag) {
        double massone;
        if (l_rmass.data()) massone = l_rmass[i];
        else massone = l_mass[l_type[i]];
        double fc0 = massone*(l_v(i,0) - v0)/dtf - l_f(i,0);
        double fc1 = massone*(l_v(i,1) - v1)/dtf - l_f(i,1);
        double fc2 = massone*(l_v(i,2) - v2)/dtf - l_f(i,2);
        double vr[6];
        vr[0] = 0.5*x0*fc0;
        vr[1] = 0.5*x1*fc1;
        vr[2] = 0.5*x2*fc2;
        vr[3] = 0.5*x0*fc1;
        vr[4] = 0.5*x0*fc2;
        vr[5] = 0.5*x1*fc2;

        int j = i;
        v_tally(1,&j,1.0,vr);
      }
    }); // end Kokkos::parallel_for
  } // end local block.


  // set orientation, omega, angmom of each extended particle
  /*
  if (extended) {
    double theta_body,theta;
    double *shape,*quatatom,*inertiaatom;
    AtomVecEllipsoid::Bonus *ebonus;
    if (avec_ellipsoid) ebonus = avec_ellipsoid->bonus;
    AtomVecLine::Bonus *lbonus;
    if (avec_line) lbonus = avec_line->bonus;
    AtomVecTri::Bonus *tbonus;
    if (avec_tri) tbonus = avec_tri->bonus;
    double **omega_one = atom->omega;
    double **angmom_one = atom->angmom;
    double **mu = atom->mu;
    int *ellipsoid = atom->ellipsoid;
    int *line = atom->line;
    int *tri = atom->tri;
    for (int i = 0; i < nlocal; i++) {
      if (body[i] < 0) continue;
      ibody = body[i];
      if (eflags[i] & SPHERE) {
        omega_one[i][0] = omega[ibody][0];
        omega_one[i][1] = omega[ibody][1];
        omega_one[i][2] = omega[ibody][2];
      } else if (eflags[i] & ELLIPSOID) {
        shape = ebonus[ellipsoid[i]].shape;
        quatatom = ebonus[ellipsoid[i]].quat;
        MathExtra::quatquat(quat[ibody],orient[i],quatatom);
        MathExtra::qnormalize(quatatom);
        ione[0] = EINERTIA*rmass[i] * (shape[1]*shape[1] + shape[2]*shape[2]);
        ione[1] = EINERTIA*rmass[i] * (shape[0]*shape[0] + shape[2]*shape[2]);
        ione[2] = EINERTIA*rmass[i] * (shape[0]*shape[0] + shape[1]*shape[1]);
        MathExtra::q_to_exyz(quatatom,exone,eyone,ezone);
        MathExtra::omega_to_angmom(omega[ibody],exone,eyone,ezone,ione,
                                   angmom_one[i]);
      } else if (eflags[i] & LINE) {
        if (quat[ibody][3] >= 0.0) theta_body = 2.0*acos(quat[ibody][0]);
        else theta_body = -2.0*acos(quat[ibody][0]);
        theta = orient[i][0] + theta_body;
        while (theta <= -MY_PI) theta += MY_2PI;
        while (theta > MY_PI) theta -= MY_2PI;
        lbonus[line[i]].theta = theta;
        omega_one[i][0] = omega[ibody][0];
        omega_one[i][1] = omega[ibody][1];
        omega_one[i][2] = omega[ibody][2];
      } else if (eflags[i] & TRIANGLE) {
        inertiaatom = tbonus[tri[i]].inertia;
        quatatom = tbonus[tri[i]].quat;
        MathExtra::quatquat(quat[ibody],orient[i],quatatom);
        MathExtra::qnormalize(quatatom);
        MathExtra::q_to_exyz(quatatom,exone,eyone,ezone);
        MathExtra::omega_to_angmom(omega[ibody],exone,eyone,ezone,
                                   inertiaatom,angmom_one[i]);
      }
      if (eflags[i] & DIPOLE) {
        MathExtra::quat_to_mat(quat[ibody],p);
        MathExtra::matvec(p,dorient[i],mu[i]);
        MathExtra::snormalize3(mu[i][3],mu[i],mu[i]);
      }
    }
  }
  */
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

  k_vcm.sync<DeviceType>();
  k_fcm.sync<DeviceType>();
  k_fflag.sync<DeviceType>();
  k_tflag.sync<DeviceType>();
  k_angmom.sync<DeviceType>();
  k_omega.sync<DeviceType>();
  k_torque.sync<DeviceType>();
  k_xcmimage.sync<DeviceType>();
  k_body.sync<DeviceType>();

  atomKK->k_x.sync<DeviceType>();
  atomKK->k_f.sync<DeviceType>();
  atomKK->k_v.sync<DeviceType>();

  k_langextra.sync<DeviceType>();

  k_sum.sync<DeviceType>();
  k_all.sync<DeviceType>();

  atomKK->k_type.sync<DeviceType>();
  atomKK->k_rmass.sync<DeviceType>();
  atomKK->k_mass.sync<DeviceType>();

  k_ex_space.sync<DeviceType>();
  k_ey_space.sync<DeviceType>();
  k_ez_space.sync<DeviceType>();

  k_xcm.sync<DeviceType>();
  k_displace.sync<DeviceType>();

  if (!earlyflag) compute_forces_and_torques();

  // update vcm and angmom
  // fflag,tflag = 0 for some dimensions in 2d

  {
    // Local block for Kokkos lambda:
    auto l_vcm = k_vcm.d_view;
    auto l_fcm = k_fcm.d_view;
    auto l_fflag = k_fflag.d_view;
    auto l_tflag = k_tflag.d_view;
    auto l_angmom = k_angmom.d_view;
    auto l_torque = k_torque.d_view;

    auto l_ex_space = k_ex_space.d_view;
    auto l_ey_space = k_ey_space.d_view;
    auto l_ez_space = k_ez_space.d_view;
    auto l_masstotal = k_masstotal.d_view;
    auto l_inertia = k_inertia.d_view;
    auto l_omega = k_omega.d_view;

    Kokkos::parallel_for(nbody, LAMMPS_LAMBDA(const int &ibody) {
      // update vcm by 1/2 step
      double dtfm = dtf / l_masstotal(ibody);

      l_vcm(ibody,0) += dtfm * l_fcm(ibody,0) * l_fflag(ibody,0);
      l_vcm(ibody,1) += dtfm * l_fcm(ibody,1) * l_fflag(ibody,1);
      l_vcm(ibody,2) += dtfm * l_fcm(ibody,2) * l_fflag(ibody,2);

      // update angular momentum by 1/2 step

      l_angmom(ibody,0) += dtf * l_torque(ibody,0) * l_tflag(ibody,0);
      l_angmom(ibody,1) += dtf * l_torque(ibody,1) * l_tflag(ibody,1);
      l_angmom(ibody,2) += dtf * l_torque(ibody,2) * l_tflag(ibody,2);

      MathExtraKokkos::angmom_to_omega(l_angmom,
                                       l_ex_space,
                                       l_ey_space,
                                       l_ez_space,
                                       l_inertia, l_omega, ibody);
    });
  }
  // set velocity/rotation of atoms in rigid bodies
  // virial is already setup from initial_integrate

  set_v_kokkos();

  // Modifies:
  // atom->v, [atom->angmom, atom->omega]
  // Needs:
  // omega, vcm, atom->v, vr, atom->omega, atom->angmom, displace, body

  // This performs set_v() on host:
  /*
  atomKK->k_v.sync<LMPHostType>();
  atomKK->k_angmom.sync<LMPHostType>();
  atomKK->k_omega.sync<LMPHostType>();
  k_omega.sync<LMPHostType>();
  k_vcm.sync<LMPHostType>();
  k_displace.sync<LMPHostType>();
  k_body.sync<LMPHostType>();

  set_v();
  atomKK->k_v.modify<LMPHostType>();
  atomKK->k_angmom.modify<LMPHostType>();
  atomKK->k_omega.modify<LMPHostType>();

  k_vcm.sync<DeviceType>();
  k_displace.sync<DeviceType>();
  k_body.sync<DeviceType>();
  */
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
  // needs:
  // sum, xcm, f, body, all, fcm, langextra, torque
  // [atom->torque]
  // modifies:
  // fcm, torque, sum,

  // sum over atoms to get force and torque on rigid body
  int nlocal = atomKK->nlocal;

  {
    // Local block for Kokkos parallel for:

    auto l_x = atomKK->k_x.d_view;
    auto l_f = atomKK->k_f.d_view;
    auto l_image = atomKK->k_image.d_view;

    auto l_body = k_body.d_view;
    auto l_xcm = k_xcm.d_view;
    auto l_fcm = k_fcm.d_view;
    auto l_torque = k_torque.d_view;

    // For the torque of extended particles:
    auto l_eflags = k_eflags.d_view;
    auto l_langextra = k_langextra.d_view;

    auto l_ex_space = k_ex_space.d_view;
    auto l_ey_space = k_ey_space.d_view;
    auto l_ez_space = k_ez_space.d_view;

    auto l_sum = k_sum.d_view;
    auto l_all = k_all.d_view;

    DomainKokkos *domainKK = static_cast<DomainKokkos*>(domain);

    Few<double,3> prd{domain->prd[0], domain->prd[1], domain->prd[2]};
    Few<double,6> h{domain->h[0], domain->h[1], domain->h[2],
                    domain->h[3], domain->h[4], domain->h[5]};

    Few<double,3> boxlo{domain->boxlo[0], domain->boxlo[1], domain->boxlo[2]};
    Few<double,3> boxhi{domain->boxhi[0], domain->boxhi[1], domain->boxhi[2]};

    Few<double,3> boxlo_lambda{domain->boxlo_lamda[0], domain->boxlo_lamda[1],
                               domain->boxlo_lamda[2]};
    Few<double,3> boxhi_lambda{domain->boxhi_lamda[0], domain->boxhi_lamda[1],
                               domain->boxhi_lamda[2]};
    Few<double,6> h_inv{domain->h_inv[0], domain->h_inv[1], domain->h_inv[2],
                        domain->h_inv[3], domain->h_inv[4], domain->h_inv[5]};
    Few<double,3> prd_lambda{domain->prd_lamda[0], domain->prd_lamda[1],
                             domain->prd_lamda[2]};
    int periodic_bits = 0;
    periodic_bits +=   domain->xperiodic;
    periodic_bits += 2*domain->yperiodic;
    periodic_bits += 4*domain->zperiodic;

    Few<double, 3> xi, unwrap;


    Kokkos::parallel_for(nlocal, LAMMPS_LAMBDA(const int &i) {
      if (l_body[i] < 0) return;
      int ibody = l_body[i];
      l_sum(ibody,0) += l_f(i,0);
      l_sum(ibody,1) += l_f(i,1);
      l_sum(ibody,2) += l_f(i,2);

      xi[0] = l_x(i,0);
      xi[1] = l_x(i,1);
      xi[2] = l_x(i,2);

      unwrap = domainKK->unmap(prd, h, triclinic,
                               xi, l_image(i));

      double dx = unwrap[0] - l_xcm(ibody,0);
      double dy = unwrap[1] - l_xcm(ibody,1);
      double dz = unwrap[2] - l_xcm(ibody,2);

      l_sum(ibody,3) += dy*l_f(i,2) - dz*l_f(i,1);
      l_sum(ibody,4) += dy*l_f(i,0) - dz*l_f(i,2);
      l_sum(ibody,5) += dy*l_f(i,1) - dz*l_f(i,0);
    });

    // If particles are extended, second loop over their torques:
    if (extended) {
      Kokkos::parallel_for(nlocal, LAMMPS_LAMBDA(const int &i) {
        if (body[i] < 0) return;
        int ibody = l_body[i];

        if (l_eflags[i] & TORQUE) {
          l_sum(ibody,3) += l_torque(i,0);
          l_sum(ibody,4) += l_torque(i,1);
          l_sum(ibody,5) += l_torque(i,2);
        }
      });
    }

    // TODO: Don't know how to rewrite this in KOKKOS parallel_reduce.
    // This is probably a  significant communication overhead:
    k_sum.sync<LMPHostType>();
    k_all.sync<LMPHostType>();
    MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);
    k_all.modify<LMPHostType>();
    k_all.sync<DeviceType>();

    // include Langevin thermostat forces
    Kokkos::parallel_for(nbody, LAMMPS_LAMBDA(const int &ibody) {
      l_fcm(ibody,0) = l_all(ibody,0) + l_langextra(ibody,0);
      l_fcm(ibody,1) = l_all(ibody,1) + l_langextra(ibody,1);
      l_fcm(ibody,2) = l_all(ibody,2) + l_langextra(ibody,2);
      l_torque(ibody,0) = l_all(ibody,3) + l_langextra(ibody,3);
      l_torque(ibody,1) = l_all(ibody,4) + l_langextra(ibody,4);
      l_torque(ibody,2) = l_all(ibody,5) + l_langextra(ibody,5);
    });
  } // ends local block for Kokkos parallel_reduces and parallel_fors.

}


template <class DeviceType>
void FixRigidKokkos<DeviceType>::set_v_kokkos()
{
  // Modifies:
  // atom->v, vr, [atom->angmom, atom->omega]
  // Needs:
  // omega, vcm, atom->v, vr, atom->omega, atom->angmom, displace, body
  //

  atomKK->k_v.sync<DeviceType>();
  atomKK->k_f.sync<DeviceType>();
  atomKK->k_x.sync<DeviceType>();
  atomKK->k_type.sync<DeviceType>();

  atomKK->k_omega.sync<DeviceType>();
  k_omega.sync<DeviceType>();
  k_vcm.sync<DeviceType>();
  k_displace.sync<DeviceType>();
  k_body.sync<DeviceType>();

  k_ex_space.sync<DeviceType>();
  k_ey_space.sync<DeviceType>();
  k_ez_space.sync<DeviceType>();
  k_displace.sync<DeviceType>();
  k_xcmimage.sync<DeviceType>();

  atomKK->k_v.modify<DeviceType>();
  atomKK->sync(execution_space, datamask_read);

  int nlocal = atomKK->nlocal;
  double xy, xz, yz;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  if (triclinic) {
    xy = domain->xy;
    xz = domain->xz;
    yz = domain->yz;
  }

  {
    // Local block for Kokkos::parallel_for.
    auto l_v = atomKK->k_v.d_view;
    auto l_x = atomKK->k_x.d_view;
    auto l_f = atomKK->k_f.d_view;

    auto l_omega = k_omega.d_view;
    auto l_body = k_body.d_view;
    auto l_rmass = atomKK->k_rmass.d_view;
    auto l_mass = atomKK->k_mass.d_view;
    auto l_type = atomKK->k_type.d_view;

    auto l_ex_space = k_ex_space.d_view;
    auto l_ey_space = k_ey_space.d_view;
    auto l_ez_space = k_ez_space.d_view;

    auto l_displace = k_displace.d_view;

    auto l_vcm = k_vcm.d_view;
    auto l_xcmimage = k_xcmimage.d_view;


    Kokkos::parallel_for(nlocal, LAMMPS_LAMBDA(const int &i) {
      if (l_body[i] < 0) return;
      const int ibody = l_body[i];

      auto l_ex_space_ibody = Kokkos::subview(l_ex_space, ibody, Kokkos::ALL);
      auto l_ey_space_ibody = Kokkos::subview(l_ey_space, ibody, Kokkos::ALL);
      auto l_ez_space_ibody = Kokkos::subview(l_ez_space, ibody, Kokkos::ALL);
      auto l_displace_i = Kokkos::subview(l_displace, i, Kokkos::ALL);

      auto delta = MathExtraKokkos::matvec(l_ex_space_ibody,
                                           l_ey_space_ibody,
                                           l_ez_space_ibody,
                                           l_displace_i);
      /*
      double delta[3];
      delta[0] = l_ex_space(ibody,0)*l_displace(i,0) +
        l_ey_space(ibody,0)*l_displace(i,1) +
        l_ez_space(ibody,0)*l_displace(i,2);
      delta[1] = l_ex_space(ibody,1)*l_displace(i,0) +
        l_ey_space(ibody,1)*l_displace(i,1) +
        l_ez_space(ibody,1)*l_displace(i,2);
      delta[2] = l_ex_space(ibody,2)*l_displace(i,0) +
        l_ey_space(ibody,2)*l_displace(i,1) +
        l_ez_space(ibody,2)*l_displace(i,2);
      */
      double v0, v1, v2;
      if (evflag) {
        v0 = l_v(i,0);
        v1 = l_v(i,1);
        v2 = l_v(i,2);
      }

      l_v(i,0) = l_omega(ibody,1)*delta[2] - l_omega(ibody,2)*delta[1]
        + l_vcm(ibody,0);
      l_v(i,1) = l_omega(ibody,2)*delta[0] - l_omega(ibody,0)*delta[2]
        + l_vcm(ibody,1);
      l_v(i,2) = l_omega(ibody,0)*delta[1] - l_omega(ibody,1)*delta[0]
        + l_vcm(ibody,2);

      double massone = 0;
      if (evflag) {
        if (l_rmass.data()) massone = l_rmass[i];
        else massone = l_mass[l_type[i]];

        double fc0 = massone*(l_v(i,0) - v0)/dtf - l_f(i,0);
        double fc1 = massone*(l_v(i,1) - v1)/dtf - l_f(i,1);
        double fc2 = massone*(l_v(i,2) - v2)/dtf - l_f(i,2);

        imageint xbox = (l_xcmimage[i] & IMGMASK) - IMGMAX;
        imageint ybox = (l_xcmimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
        imageint zbox = (l_xcmimage[i] >> IMG2BITS) - IMGMAX;

        double x0, x1, x2;
        if (triclinic == 0) {
          x0 = l_x(i,0) + xbox*xprd;
          x1 = l_x(i,1) + ybox*yprd;
          x2 = l_x(i,2) + zbox*zprd;
        } else {
          x0 = l_x(i,0) + xbox*xprd + ybox*xy + zbox*xz;
          x1 = l_x(i,1) + ybox*yprd + zbox*yz;
          x2 = l_x(i,2) + zbox*zprd;
        }

        double vr[6];
        vr[0] = 0.5*x0*fc0;
        vr[1] = 0.5*x1*fc1;
        vr[2] = 0.5*x2*fc2;
        vr[3] = 0.5*x0*fc1;
        vr[4] = 0.5*x0*fc2;
        vr[5] = 0.6*x1*fc2;

        int j = i;
        v_tally(1,&j, 1.0, vr);
      }
    });
  }

  /*
  if (extended) {
    double *shape,*quatatom,*inertiaatom;

    AtomVecEllipsoid::Bonus *ebonus;
    if (avec_ellipsoid) ebonus = avec_ellipsoid->bonus;
    AtomVecTri::Bonus *tbonus;
    if (avec_tri) tbonus = avec_tri->bonus;
    double **omega_one = atom->omega;
    double **angmom_one = atom->angmom;
    int *ellipsoid = atom->ellipsoid;
    int *tri = atom->tri;

    for (int i = 0; i < nlocal; i++) {
      if (body[i] < 0) continue;
      const int ibody = body[i];

      if (eflags[i] & SPHERE) {
        omega_one[i][0] = omega[ibody][0];
        omega_one[i][1] = omega[ibody][1];
        omega_one[i][2] = omega[ibody][2];
      } else if (eflags[i] & ELLIPSOID) {
        shape = ebonus[ellipsoid[i]].shape;
        quatatom = ebonus[ellipsoid[i]].quat;
        ione[0] = EINERTIA*rmass[i] * (shape[1]*shape[1] + shape[2]*shape[2]);
        ione[1] = EINERTIA*rmass[i] * (shape[0]*shape[0] + shape[2]*shape[2]);
        ione[2] = EINERTIA*rmass[i] * (shape[0]*shape[0] + shape[1]*shape[1]);
        MathExtra::q_to_exyz(quatatom,exone,eyone,ezone);
        MathExtra::omega_to_angmom(omega[ibody],exone,eyone,ezone,ione,
                                   angmom_one[i]);
      } else if (eflags[i] & LINE) {
        omega_one[i][0] = omega[ibody][0];
        omega_one[i][1] = omega[ibody][1];
        omega_one[i][2] = omega[ibody][2];
      } else if (eflags[i] & TRIANGLE) {
        inertiaatom = tbonus[tri[i]].inertia;
        quatatom = tbonus[tri[i]].quat;
        MathExtra::q_to_exyz(quatatom,exone,eyone,ezone);
        MathExtra::omega_to_angmom(omega[ibody],exone,eyone,ezone,
                                   inertiaatom,angmom_one[i]);
      }
    }
  }
  */
}


template <class DeviceType>
void FixRigidKokkos<DeviceType>::post_force(int vflag)
{
  if (langflag) apply_langevin_thermostat_kokkos();
  if (earlyflag) compute_forces_and_torques_kokkos();
}


// Closely mirrors FixRigid::apply_langevin_thermostat.
template <class DeviceType>
void FixRigidKokkos<DeviceType>::apply_langevin_thermostat_kokkos()
{

  k_masstotal.sync<DeviceType>();
  k_langextra.sync<DeviceType>();
  k_vcm.sync<DeviceType>();
  k_inertia.sync<DeviceType>();
  k_omega.sync<DeviceType>();

  if (comm->me == 0) {
    // Local block for Kokkos lambda coincides with if-block:

    double delta = update->ntimestep - update->beginstep;
    if (delta != 0.0) delta /= update->endstep - update->beginstep;
    t_target = t_start + delta * (t_stop-t_start);
    double tsqrt = sqrt(t_target);

    double boltz = force->boltz;
    double dt = update->dt;
    double mvv2e = force->mvv2e;
    double ftm2v = force->ftm2v;

    auto l_masstotal = k_masstotal.d_view;
    auto l_langextra = k_langextra.d_view;
    auto l_vcm = k_vcm.d_view;
    auto l_inertia = k_inertia.d_view;
    auto l_omega = k_omega.d_view;

    // Aped from fix_langevin_kokkos:
    Kokkos::parallel_for(nbody, LAMMPS_LAMBDA(const int&ibody) {

      rand_type rand_gen = rand_pool.get_state();

      double gamma1 = -l_masstotal[ibody] / t_period / ftm2v;
      double gamma2 = sqrt(l_masstotal[ibody]) * tsqrt *
        sqrt(24.0*boltz/t_period/dt/mvv2e) / ftm2v;
      l_langextra(ibody,0) = gamma1*l_vcm(ibody,0) + gamma2*(rand_gen.drand()-0.5);
      l_langextra(ibody,1) = gamma1*l_vcm(ibody,1) + gamma2*(rand_gen.drand()-0.5);
      l_langextra(ibody,2) = gamma1*l_vcm(ibody,2) + gamma2*(rand_gen.drand()-0.5);

      gamma1 = -1.0 / t_period / ftm2v;
      gamma2 = tsqrt * sqrt(24.0*boltz/t_period/dt/mvv2e) / ftm2v;
      l_langextra(ibody,3) = l_inertia(ibody,0)*gamma1*l_omega(ibody,0) +
        sqrt(l_inertia(ibody,0))*gamma2*(random->uniform()-0.5);
      l_langextra(ibody,4) = l_inertia(ibody,1)*gamma1*l_omega(ibody,2) +
        sqrt(l_inertia(ibody,1))*gamma2*(random->uniform()-0.5);
      l_langextra(ibody,5) = l_inertia(ibody,2)*gamma1*l_omega(ibody,2) +
        sqrt(l_inertia(ibody,2))*gamma2*(random->uniform()-0.5);

      rand_pool.free_state(rand_gen);
    });
  }

  k_langextra.modify<DeviceType>();
  MPI_Bcast(&langextra[0][0],6*nbody,MPI_DOUBLE,0,world);
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
