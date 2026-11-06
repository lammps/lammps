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

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "fix_rigid_omp.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "rigid_const.h"

#include <cmath>
#include <cstring>

#if defined(_OPENMP)
#include <omp.h>
#endif
#include "omp_compat.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace RigidConst;

namespace {
using dbl3_t = struct {
  double x, y, z;
};
}

/* ---------------------------------------------------------------------- */

void FixRigidOMP::initial_integrate(int vflag)
{
#if defined(_OPENMP)
#pragma omp parallel for LMP_DEFAULT_NONE schedule(static)
#endif
  for (int ibody = 0; ibody < nbody; ibody++) {

    // update vcm by 1/2 step

    const double dtfm = dtf / masstotal[ibody];
    vcm[ibody][0] += dtfm * fcm[ibody][0] * fflag[ibody][0];
    vcm[ibody][1] += dtfm * fcm[ibody][1] * fflag[ibody][1];
    vcm[ibody][2] += dtfm * fcm[ibody][2] * fflag[ibody][2];

    // update xcm by full step

    xcm[ibody][0] += dtv * vcm[ibody][0];
    xcm[ibody][1] += dtv * vcm[ibody][1];
    xcm[ibody][2] += dtv * vcm[ibody][2];

    // update angular momentum by 1/2 step

    angmom[ibody][0] += dtf * torque[ibody][0] * tflag[ibody][0];
    angmom[ibody][1] += dtf * torque[ibody][1] * tflag[ibody][1];
    angmom[ibody][2] += dtf * torque[ibody][2] * tflag[ibody][2];

    // compute omega at 1/2 step from angmom at 1/2 step and current q
    // update quaternion a full step via Richardson iteration
    // returns new normalized quaternion, also updated omega at 1/2 step
    // update ex,ey,ez to reflect new quaternion

    MathExtra::angmom_to_omega(angmom[ibody], ex_space[ibody], ey_space[ibody], ez_space[ibody],
                               inertia[ibody], omega[ibody]);
    MathExtra::richardson(quat[ibody], angmom[ibody], omega[ibody], inertia[ibody], dtq);
    MathExtra::q_to_exyz(quat[ibody], ex_space[ibody], ey_space[ibody], ez_space[ibody]);
  }    // end of omp parallel for

  // virial setup before call to set_xv

  v_init(vflag);

  // set coords/orient and velocity/rotation of atoms in rigid bodies
  // from quarternion and omega

  if (domain->dimension == 2) {
    if (triclinic) {
      if (evflag)
        set_xv_thr<1,1,2>();
      else
        set_xv_thr<1,0,2>();
    } else {
      if (evflag)
        set_xv_thr<0,1,2>();
      else
        set_xv_thr<0,0,2>();
    }
  } else {

    if (triclinic) {
      if (evflag)
        set_xv_thr<1,1,3>();
      else
        set_xv_thr<1,0,3>();
    } else {
      if (evflag)
        set_xv_thr<0,1,3>();
      else
        set_xv_thr<0,0,3>();
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixRigidOMP::compute_forces_and_torques()
{
  double * const * _noalias const x = atom->x;
  const auto * _noalias const f = (dbl3_t *) atom->f[0];
  const double * const * const torque_one = atom->torque;
  const int nlocal = atom->nlocal;

  // sum over atoms to get force and torque on rigid body
  // we have 3 different strategies for multi-threading this.

   if (rstyle == SINGLE) {
     // we have just one rigid body. use OpenMP reduction to get sum[]
     double s0=0.0,s1=0.0,s2=0.0,s3=0.0,s4=0.0,s5=0.0;

#if defined(_OPENMP)
#pragma omp parallel for LMP_DEFAULT_NONE reduction(+:s0,s1,s2,s3,s4,s5)
#endif
     for (int i = 0; i < nlocal; i++) {
       const int ibody = body[i];
       if (ibody < 0) continue;

       double unwrap[3];
       domain->unmap(x[i],xcmimage[i],unwrap);
       const double dx = unwrap[0] - xcm[0][0];
       const double dy = unwrap[1] - xcm[0][1];
       const double dz = unwrap[2] - xcm[0][2];

       s0 += f[i].x;
       s1 += f[i].y;
       s2 += f[i].z;

       s3 += dy*f[i].z - dz*f[i].y;
       s4 += dz*f[i].x - dx*f[i].z;
       s5 += dx*f[i].y - dy*f[i].x;

       if (extended && (eflags[i] & TORQUE)) {
         s3 += torque_one[i][0];
         s4 += torque_one[i][1];
         s5 += torque_one[i][2];
       }
     }
     sum[0][0]=s0; sum[0][1]=s1; sum[0][2]=s2;
     sum[0][3]=s3; sum[0][4]=s4; sum[0][5]=s5;

  } else if (rstyle == GROUP) {

     // we likely have a rather small number of groups so we loop
     // over bodies and thread over all atoms for each of them.

     for (int ib = 0; ib < nbody; ++ib) {
       double s0=0.0,s1=0.0,s2=0.0,s3=0.0,s4=0.0,s5=0.0;

#if defined(_OPENMP)
#pragma omp parallel for LMP_DEFAULT_NONE LMP_SHARED(ib) reduction(+:s0,s1,s2,s3,s4,s5)
#endif
       for (int i = 0; i < nlocal; i++) {
         const int ibody = body[i];
         if (ibody != ib) continue;

         s0 += f[i].x;
         s1 += f[i].y;
         s2 += f[i].z;

         double unwrap[3];
         domain->unmap(x[i],xcmimage[i],unwrap);
         const double dx = unwrap[0] - xcm[ibody][0];
         const double dy = unwrap[1] - xcm[ibody][1];
         const double dz = unwrap[2] - xcm[ibody][2];

         s3 += dy*f[i].z - dz*f[i].y;
         s4 += dz*f[i].x - dx*f[i].z;
         s5 += dx*f[i].y - dy*f[i].x;

         if (extended && (eflags[i] & TORQUE)) {
           s3 += torque_one[i][0];
           s4 += torque_one[i][1];
           s5 += torque_one[i][2];
         }
       }

       sum[ib][0]=s0; sum[ib][1]=s1; sum[ib][2]=s2;
       sum[ib][3]=s3; sum[ib][4]=s4; sum[ib][5]=s5;
     }

  } else if (rstyle == MOLECULE) {

     // we likely have a large number of rigid objects with only a
     // a few atoms each. so we loop over all atoms for all threads
     // and then each thread only processes some bodies.

     memset(&sum[0][0],0,6*nbody*sizeof(double));

#if defined(_OPENMP)
     const int nthreads=comm->nthreads;
#pragma omp parallel LMP_DEFAULT_NONE
#else
     const int nthreads=1;
#endif
     {
#if defined(_OPENMP)
       const int tid = omp_get_thread_num();
#else
       const int tid = 0;
#endif

       for (int i = 0; i < nlocal; i++) {
         const int ibody = body[i];
         if ((ibody < 0) || (ibody % nthreads != tid)) continue;

         double unwrap[3];
         domain->unmap(x[i],xcmimage[i],unwrap);
         const double dx = unwrap[0] - xcm[ibody][0];
         const double dy = unwrap[1] - xcm[ibody][1];
         const double dz = unwrap[2] - xcm[ibody][2];

         const double s0 = f[i].x;
         const double s1 = f[i].y;
         const double s2 = f[i].z;

         double s3 = dy*s2 - dz*s1;
         double s4 = dz*s0 - dx*s2;
         double s5 = dx*s1 - dy*s0;

         if (extended && (eflags[i] & TORQUE)) {
           s3 += torque_one[i][0];
           s4 += torque_one[i][1];
           s5 += torque_one[i][2];
         }

         sum[ibody][0] += s0; sum[ibody][1] += s1; sum[ibody][2] += s2;
         sum[ibody][3] += s3; sum[ibody][4] += s4; sum[ibody][5] += s5;
       }
     }
   } else
     error->all(FLERR,"rigid style is unsupported by fix rigid/omp");

  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);

  // update vcm and angmom
  // include Langevin thermostat forces
  // fflag,tflag = 0 for some dimensions in 2d

#if defined(_OPENMP)
#pragma omp parallel for LMP_DEFAULT_NONE schedule(static)
#endif
  for (int ibody = 0; ibody < nbody; ibody++) {
    fcm[ibody][0] = all[ibody][0];
    fcm[ibody][1] = all[ibody][1];
    fcm[ibody][2] = all[ibody][2];
    torque[ibody][0] = all[ibody][3];
    torque[ibody][1] = all[ibody][4];
    torque[ibody][2] = all[ibody][5];
  }

  // add langevin friction to force and torque of each body

  if (langflag) {
#if defined(_OPENMP)
#pragma omp parallel for LMP_DEFAULT_NONE schedule(static)
#endif
    for (int ibody = 0; ibody < nbody; ibody++) {
      fcm[ibody][0] += fflag[ibody][0]*langextra[ibody][0];
      fcm[ibody][1] += fflag[ibody][1]*langextra[ibody][1];
      fcm[ibody][2] += fflag[ibody][2]*langextra[ibody][2];
      torque[ibody][0] += tflag[ibody][0]*langextra[ibody][3];
      torque[ibody][1] += tflag[ibody][1]*langextra[ibody][4];
      torque[ibody][2] += tflag[ibody][2]*langextra[ibody][5];
    }
  }

  // add gravity force to COM of each body

  if (id_gravity) {
#if defined(_OPENMP)
#pragma omp parallel for LMP_DEFAULT_NONE schedule(static)
#endif
    for (int ibody = 0; ibody < nbody; ibody++) {
      fcm[ibody][0] += gvec[0]*masstotal[ibody];
      fcm[ibody][1] += gvec[1]*masstotal[ibody];
      fcm[ibody][2] += gvec[2]*masstotal[ibody];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixRigidOMP::final_integrate()
{
  if (!earlyflag) compute_forces_and_torques();
  if (domain->dimension == 2) enforce2d();

  // update vcm and angmom

#if defined(_OPENMP)
#pragma omp parallel for LMP_DEFAULT_NONE schedule(static)
#endif
  for (int ibody = 0; ibody < nbody; ibody++) {

    // update vcm by 1/2 step

    const double dtfm = dtf / masstotal[ibody];
    vcm[ibody][0] += dtfm * fcm[ibody][0] * fflag[ibody][0];
    vcm[ibody][1] += dtfm * fcm[ibody][1] * fflag[ibody][1];
    vcm[ibody][2] += dtfm * fcm[ibody][2] * fflag[ibody][2];

    // update angular momentum by 1/2 step

    angmom[ibody][0] += dtf * torque[ibody][0] * tflag[ibody][0];
    angmom[ibody][1] += dtf * torque[ibody][1] * tflag[ibody][1];
    angmom[ibody][2] += dtf * torque[ibody][2] * tflag[ibody][2];

    MathExtra::angmom_to_omega(angmom[ibody],ex_space[ibody],ey_space[ibody],
                               ez_space[ibody],inertia[ibody],omega[ibody]);
  }

  // set velocity/rotation of atoms in rigid bodies
  // virial is already setup from initial_integrate
  // triclinic only matters for virial calculation.

#if defined(_OPENMP)
  if (domain->dimension == 2) {
    if (evflag)
      if (triclinic)
        set_v_thr<1,1,2>();
      else
        set_v_thr<0,1,2>();
    else
      set_v_thr<0,0,2>();
  } else {
    if (evflag)
      if (triclinic)
        set_v_thr<1,1,3>();
      else
        set_v_thr<0,1,3>();
    else
      set_v_thr<0,0,3>();
  }
#else
  set_v();
#endif
}

/* ---------------------------------------------------------------------- */

void FixRigidOMP::compute_accelerations()
{
#if defined(_OPENMP)
#pragma omp parallel for LMP_DEFAULT_NONE schedule(static)
#endif
  for (int ibody = 0; ibody < nbody; ibody++) {
    double omegadot_body[3], wbody[3], tbody[3], tspace[3];
    double *ex, *ey, *ez, *langone;

    if(langflag) {
      langone = langextra[ibody];
      acc_vir[ibody][0] = fflag[ibody][0] * (fcm[ibody][0] - langone[0]) / masstotal[ibody];
      acc_vir[ibody][1] = fflag[ibody][1] * (fcm[ibody][1] - langone[1]) / masstotal[ibody];
      if (domain->dimension == 2) acc_vir[ibody][2] = 0.0;
      else acc_vir[ibody][2] = fflag[ibody][2] * (fcm[ibody][2] - langone[2]) / masstotal[ibody];
    } else {
      acc_vir[ibody][0] = fflag[ibody][0] * fcm[ibody][0] / masstotal[ibody];
      acc_vir[ibody][1] = fflag[ibody][1] * fcm[ibody][1] / masstotal[ibody];
      if (domain->dimension == 2) acc_vir[ibody][2] = 0.0;
      else acc_vir[ibody][2] = fflag[ibody][2] * fcm[ibody][2] / masstotal[ibody];
    }

    ex = ex_space[ibody], ey = ey_space[ibody], ez = ez_space[ibody];
    wbody[0] = omega[ibody][0]*ex[0] + omega[ibody][1]*ex[1] + omega[ibody][2]*ex[2];
    wbody[1] = omega[ibody][0]*ey[0] + omega[ibody][1]*ey[1] + omega[ibody][2]*ey[2];
    wbody[2] = omega[ibody][0]*ez[0] + omega[ibody][1]*ez[1] + omega[ibody][2]*ez[2];
    if(langflag) {
      tspace[0] = tflag[ibody][0] * (torque[ibody][0] - langone[3]);
      tspace[1] = tflag[ibody][1] * (torque[ibody][1] - langone[4]);
      tspace[2] = tflag[ibody][2] * (torque[ibody][2] - langone[5]);
    } else {
      tspace[0] = tflag[ibody][0] * torque[ibody][0];
      tspace[1] = tflag[ibody][1] * torque[ibody][1];
      tspace[2] = tflag[ibody][2] * torque[ibody][2];
    }
    tbody[0] = tspace[0]*ex[0] + tspace[1]*ex[1] + tspace[2]*ex[2];
    tbody[1] = tspace[0]*ey[0] + tspace[1]*ey[1] + tspace[2]*ey[2];
    tbody[2] = tspace[0]*ez[0] + tspace[1]*ez[1] + tspace[2]*ez[2];
    if (inertia[ibody][0] == 0.0) omegadot_body[0] = 0.0;
    else omegadot_body[0] = (force->ftm2v*tbody[0] + (inertia[ibody][1] - inertia[ibody][2]) * wbody[1] * wbody[2]) / inertia[ibody][0];
    if (inertia[ibody][1] == 0.0) omegadot_body[1] = 0.0;
    else omegadot_body[1] = (force->ftm2v*tbody[1] + (inertia[ibody][2] - inertia[ibody][0]) * wbody[2] * wbody[0]) / inertia[ibody][1];
    if (inertia[ibody][2] == 0.0) omegadot_body[2] = 0.0;
    else omegadot_body[2] = (force->ftm2v*tbody[2] + (inertia[ibody][0] - inertia[ibody][1]) * wbody[0] * wbody[1]) / inertia[ibody][2];
    if (domain->dimension == 2) {
      acc_vir[ibody][3] = 0.0;
      acc_vir[ibody][4] = 0.0;
    } else {
      acc_vir[ibody][3] = omegadot_body[0]*ex[0] + omegadot_body[1]*ey[0] + omegadot_body[2]*ez[0];
      acc_vir[ibody][4] = omegadot_body[0]*ex[1] + omegadot_body[1]*ey[1] + omegadot_body[2]*ez[1];
    }
    acc_vir[ibody][5] = omegadot_body[0]*ex[2] + omegadot_body[1]*ey[2] + omegadot_body[2]*ez[2];
  }    // end of omp parallel for
}

/* ----------------------------------------------------------------------
   set space-frame coords and velocity of each atom in each rigid body
   set orientation and rotation of extended particles
   x = Q displace + Xcm, mapped back to periodic box
   v = Vcm + (W cross (x - Xcm))

   NOTE: this needs to be kept in sync with FixRigidNHOMP
------------------------------------------------------------------------- */
template <int TRICLINIC, int EVFLAG, int DIMENSION>
void FixRigidOMP::set_xv_thr()
{
  auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const v = (dbl3_t *) atom->v[0];
  const double * _noalias const rmass = atom->rmass;

  const double xprd = domain->xprd;
  const double yprd = domain->yprd;
  const double zprd = domain->zprd;
  const double xy = domain->xy;
  const double xz = domain->xz;
  const double yz = domain->yz;

  // set x and v of each atom

  const int nlocal = atom->nlocal;

#if defined(_OPENMP)
#pragma omp parallel for LMP_DEFAULT_NONE schedule(static)
#endif
  for (int i = 0; i < nlocal; i++) {
    const int ibody = body[i];
    if (ibody < 0) continue;

    const auto &xcmi = * ((dbl3_t *) xcm[ibody]);
    const auto &vcmi = * ((dbl3_t *) vcm[ibody]);
    const auto &omegai = * ((dbl3_t *) omega[ibody]);

    const int xbox = (xcmimage[i] & IMGMASK) - IMGMAX;
    const int ybox = (xcmimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
    const int zbox = (xcmimage[i] >> IMG2BITS) - IMGMAX;
    const double deltax = xbox*xprd + (TRICLINIC ? ybox*xy + zbox*xz : 0.0);
    const double deltay = ybox*yprd + (TRICLINIC ? zbox*yz : 0.0);
    const double deltaz = zbox*zprd;

    // x = displacement from center-of-mass, based on body orientation
    // v = vcm + omega around center-of-mass

    MathExtra::matvec(ex_space[ibody],ey_space[ibody],
                      ez_space[ibody],displace[i],&x[i].x);

    v[i].x = omegai.y*x[i].z - omegai.z*x[i].y + vcmi.x;
    v[i].y = omegai.z*x[i].x - omegai.x*x[i].z + vcmi.y;
    v[i].z = omegai.x*x[i].y - omegai.y*x[i].x + vcmi.z;

    if (DIMENSION == 2) x[i].z = v[i].z = 0.0;

    // add center of mass to displacement
    // map back into periodic box via xbox,ybox,zbox
    // for triclinic, add in box tilt factors as well

    x[i].x += xcmi.x - deltax;
    x[i].y += xcmi.y - deltay;
    x[i].z += xcmi.z - deltaz;
  }

  // set orientation, omega, angmom of each extended particle
  // XXX: extended particle info not yet multi-threaded

  if (extended) {
    double *shape,*quatatom,*inertiaatom;
    double theta_body,theta;
    double ione[3],exone[3],eyone[3],ezone[3],p[3][3];

    AtomVecEllipsoid::Bonus *ebonus = nullptr;
    if (avec_ellipsoid) ebonus = avec_ellipsoid->bonus;
    AtomVecLine::Bonus *lbonus = nullptr;
    if (avec_line) lbonus = avec_line->bonus;
    AtomVecTri::Bonus *tbonus = nullptr;
    if (avec_tri) tbonus = avec_tri->bonus;
    double **omega_one = atom->omega;
    double **angmom_one = atom->angmom;
    double **mu = atom->mu;
    int *ellipsoid = atom->ellipsoid;
    int *line = atom->line;
    int *tri = atom->tri;

    for (int i = 0; i < nlocal; i++) {
      const int ibody = body[i];
      if (ibody < 0) continue;

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
}

/* ----------------------------------------------------------------------
   set space-frame velocity of each atom in a rigid body
   set omega and angmom of extended particles
   v = Vcm + (W cross (x - Xcm))

   NOTE: this needs to be kept in sync with FixRigidNHOMP
------------------------------------------------------------------------- */
template <int TRICLINIC, int EVFLAG, int DIMENSION>
void FixRigidOMP::set_v_thr()
{
  auto * _noalias const v = (dbl3_t *) atom->v[0];
  const auto * _noalias const f = (dbl3_t *) atom->f[0];
  const double * _noalias const rmass = atom->rmass;
  const double * _noalias const mass = atom->mass;
  const int * _noalias const type = atom->type;

  double v0=0.0,v1=0.0,v2=0.0,v3=0.0,v4=0.0,v5=0.0;

  // set acm and omegadot to compute constraints virial

  if (EVFLAG) compute_accelerations();

  // set v of each atom

  const int nlocal = atom->nlocal;

#if defined(_OPENMP)
#pragma omp parallel for LMP_DEFAULT_NONE reduction(+:v0,v1,v2,v3,v4,v5)
#endif
  for (int i = 0; i < nlocal; i++) {
    const int ibody = body[i];
    if (ibody < 0) continue;

    const auto &vcmi = * ((dbl3_t *) vcm[ibody]);
    const double *omegai = omega[ibody];
    double delta[3];

    MathExtra::matvec(ex_space[ibody],ey_space[ibody],
                      ez_space[ibody],displace[i],delta);

    v[i].x = omegai[1]*delta[2] - omegai[2]*delta[1] + vcmi.x;
    v[i].y = omegai[2]*delta[0] - omegai[0]*delta[2] + vcmi.y;
    v[i].z = omegai[0]*delta[1] - omegai[1]*delta[0] + vcmi.z;

    if (DIMENSION == 2) v[i].z = 0.0;


    // virial = unwrapped coords dotted into body constraint force
    // body constraint force = implied force from acceleration minus f external
    // assume f does not include forces internal to body
    // total acceleration does not include eventual langevin contributions
    // assume per-atom contribution is due to constraint force on that atom

    if (EVFLAG) {
      const double *acmi = acc_vir[ibody];
      const double *omegadoti = acc_vir[ibody] + 3;
      double massone, acc_rot[3], v_rot[3], acc_centr[3], fc[3], vr[6];

      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];

      MathExtra::cross3(omegadoti, delta, acc_rot);
      MathExtra::cross3(omegai, delta, v_rot) ;
      MathExtra::cross3(omegai, v_rot, acc_centr) ;
      fc[0] = massone*(acmi[0] + (acc_rot[0] + acc_centr[0])/force->ftm2v) - f[i].x;
      fc[1] = massone*(acmi[1] + (acc_rot[1] + acc_centr[1])/force->ftm2v) - f[i].y;
      if (DIMENSION == 2) fc[2] = 0.0;
      else fc[2] = massone*(acmi[2] + (acc_rot[2] + acc_centr[2])/force->ftm2v) - f[i].z;

      // if id_gravity=1 fc will also contain the gravitational field contribution

      const double x0 = delta[0] + xcm[ibody][0];
      const double x1 = delta[1] + xcm[ibody][1];
      const double x2 = delta[2] + xcm[ibody][2];

      vr[0] = x0*fc[0]; vr[1] = x1*fc[1]; vr[2] = x2*fc[2];
      vr[3] = x0*fc[1]; vr[4] = x0*fc[2]; vr[5] = x1*fc[2];

      // Fix::v_tally() is not thread safe, so we do this manually here
      // accumulate global virial into thread-local variables and reduce them later
      if (vflag_global) {
        v0 += vr[0];
        v1 += vr[1];
        v2 += vr[2];
        v3 += vr[3];
        v4 += vr[4];
        v5 += vr[5];
      }

      // accumulate per atom virial directly since we parallelize over atoms.
      if (vflag_atom) {
        vatom[i][0] += vr[0];
        vatom[i][1] += vr[1];
        vatom[i][2] += vr[2];
        vatom[i][3] += vr[3];
        vatom[i][4] += vr[4];
        vatom[i][5] += vr[5];
      }
    }
  } // end of parallel for

  // second part of thread safe virial accumulation
  // add global virial component after it was reduced across all threads
  if (EVFLAG) {
    if (vflag_global) {
      virial[0] += v0;
      virial[1] += v1;
      virial[2] += v2;
      virial[3] += v3;
      virial[4] += v4;
      virial[5] += v5;
    }
  }

  // set omega, angmom of each extended particle
  // XXX: extended particle info not yet multi-threaded

  if (extended) {
    double *shape,*quatatom,*inertiaatom;
    double ione[3],exone[3],eyone[3],ezone[3];

    AtomVecEllipsoid::Bonus *ebonus = nullptr;
    if (avec_ellipsoid) ebonus = avec_ellipsoid->bonus;
    AtomVecTri::Bonus *tbonus = nullptr;
    if (avec_tri) tbonus = avec_tri->bonus;
    double **omega_one = atom->omega;
    double **angmom_one = atom->angmom;
    int *ellipsoid = atom->ellipsoid;
    int *tri = atom->tri;

    for (int i = 0; i < nlocal; i++) {
      const int ibody = body[i];
      if (ibody < 0) continue;

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
}
