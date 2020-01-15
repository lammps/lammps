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

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "fix_rigid_omp.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "comm.h"
#include "error.h"
#include "domain.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "math_extra.h"
#include "math_const.h"
#include "rigid_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace RigidConst;

typedef struct { double x,y,z; } dbl3_t;

/* ---------------------------------------------------------------------- */

void FixRigidOMP::initial_integrate(int vflag)
{
#if defined(_OPENMP)
#pragma omp parallel for default(none) schedule(static)
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

    MathExtra::angmom_to_omega(angmom[ibody],ex_space[ibody],ey_space[ibody],
                               ez_space[ibody],inertia[ibody],omega[ibody]);
    MathExtra::richardson(quat[ibody],angmom[ibody],omega[ibody],
                          inertia[ibody],dtq);
    MathExtra::q_to_exyz(quat[ibody],
                         ex_space[ibody],ey_space[ibody],ez_space[ibody]);
  } // end of omp parallel for

  // virial setup before call to set_xv

  if (vflag) v_setup(vflag);
  else evflag = 0;

  // set coords/orient and velocity/rotation of atoms in rigid bodies
  // from quarternion and omega

  if (triclinic)
    if (evflag)
      set_xv_thr<1,1>();
    else
      set_xv_thr<1,0>();
  else
    if (evflag)
      set_xv_thr<0,1>();
    else
      set_xv_thr<0,0>();
}

/* ---------------------------------------------------------------------- */

void FixRigidOMP::compute_forces_and_torques()
{
  double * const * _noalias const x = atom->x;
  const dbl3_t * _noalias const f = (dbl3_t *) atom->f[0];
  const double * const * const torque_one = atom->torque;
  const int nlocal = atom->nlocal;

  // sum over atoms to get force and torque on rigid body
  // we have 3 different strategies for multi-threading this.

   if (rstyle == SINGLE) {
     // we have just one rigid body. use OpenMP reduction to get sum[]
     double s0=0.0,s1=0.0,s2=0.0,s3=0.0,s4=0.0,s5=0.0;

#if defined(_OPENMP)
#pragma omp parallel for default(none) reduction(+:s0,s1,s2,s3,s4,s5)
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

     // we likely have only a rather number of groups so we loops
     // over bodies and thread over all atoms for each of them.

     for (int ib = 0; ib < nbody; ++ib) {
       double s0=0.0,s1=0.0,s2=0.0,s3=0.0,s4=0.0,s5=0.0;

#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(ib) reduction(+:s0,s1,s2,s3,s4,s5)
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

     const int nthreads=comm->nthreads;
     memset(&sum[0][0],0,6*nbody*sizeof(double));

#if defined(_OPENMP)
#pragma omp parallel default(none)
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
#pragma omp parallel for default(none) schedule(static)
#endif
  for (int ibody = 0; ibody < nbody; ibody++) {
    fcm[ibody][0] = all[ibody][0] + langextra[ibody][0];
    fcm[ibody][1] = all[ibody][1] + langextra[ibody][1];
    fcm[ibody][2] = all[ibody][2] + langextra[ibody][2];
    torque[ibody][0] = all[ibody][3] + langextra[ibody][3];
    torque[ibody][1] = all[ibody][4] + langextra[ibody][4];
    torque[ibody][2] = all[ibody][5] + langextra[ibody][5];
  }

  // add gravity force to COM of each body

  if (id_gravity) {
#if defined(_OPENMP)
#pragma omp parallel for default(none) schedule(static)
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

  // update vcm and angmom

#if defined(_OPENMP)
#pragma omp parallel for default(none) schedule(static)
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

  if (evflag)
    if (triclinic)
      set_v_thr<1,1>();
    else
      set_v_thr<0,1>();
  else
    set_v_thr<0,0>();
}

/* ----------------------------------------------------------------------
   set space-frame coords and velocity of each atom in each rigid body
   set orientation and rotation of extended particles
   x = Q displace + Xcm, mapped back to periodic box
   v = Vcm + (W cross (x - Xcm))

   NOTE: this needs to be kept in sync with FixRigidNHOMP
------------------------------------------------------------------------- */
template <int TRICLINIC, int EVFLAG>
void FixRigidOMP::set_xv_thr()
{
  dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const v = (dbl3_t *) atom->v[0];
  const dbl3_t * _noalias const f = (dbl3_t *) atom->f[0];
  const double * _noalias const rmass = atom->rmass;
  const double * _noalias const mass = atom->mass;
  const int * _noalias const type = atom->type;

  double v0=0.0,v1=0.0,v2=0.0,v3=0.0,v4=0.0,v5=0.0;

  const double xprd = domain->xprd;
  const double yprd = domain->yprd;
  const double zprd = domain->zprd;
  const double xy = domain->xy;
  const double xz = domain->xz;
  const double yz = domain->yz;

  // set x and v of each atom

  const int nlocal = atom->nlocal;

#if defined(_OPENMP)
#pragma omp parallel for default(none) reduction(+:v0,v1,v2,v3,v4,v5)
#endif
  for (int i = 0; i < nlocal; i++) {
    const int ibody = body[i];
    if (ibody < 0) continue;

    const dbl3_t &xcmi = * ((dbl3_t *) xcm[ibody]);
    const dbl3_t &vcmi = * ((dbl3_t *) vcm[ibody]);
    const dbl3_t &omegai = * ((dbl3_t *) omega[ibody]);

    const int xbox = (xcmimage[i] & IMGMASK) - IMGMAX;
    const int ybox = (xcmimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
    const int zbox = (xcmimage[i] >> IMG2BITS) - IMGMAX;
    const double deltax = xbox*xprd + (TRICLINIC ? ybox*xy + zbox*xz : 0.0);
    const double deltay = ybox*yprd + (TRICLINIC ? zbox*yz : 0.0);
    const double deltaz = zbox*zprd;

    // save old positions and velocities for virial
    double x0,x1,x2,vx,vy,vz;
    if (EVFLAG) {
      x0 = x[i].x + deltax;
      x1 = x[i].y + deltay;
      x2 = x[i].z + deltaz;
      vx = v[i].x;
      vy = v[i].y;
      vz = v[i].z;
    }

    // x = displacement from center-of-mass, based on body orientation
    // v = vcm + omega around center-of-mass

    MathExtra::matvec(ex_space[ibody],ey_space[ibody],
                      ez_space[ibody],displace[i],&x[i].x);

    v[i].x = omegai.y*x[i].z - omegai.z*x[i].y + vcmi.x;
    v[i].y = omegai.z*x[i].x - omegai.x*x[i].z + vcmi.y;
    v[i].z = omegai.x*x[i].y - omegai.y*x[i].x + vcmi.z;

    // add center of mass to displacement
    // map back into periodic box via xbox,ybox,zbox
    // for triclinic, add in box tilt factors as well

    x[i].x += xcmi.x - deltax;
    x[i].y += xcmi.y - deltay;
    x[i].z += xcmi.z - deltaz;

    // virial = unwrapped coords dotted into body constraint force
    // body constraint force = implied force due to v change minus f external
    // assume f does not include forces internal to body
    // 1/2 factor b/c final_integrate contributes other half
    // assume per-atom contribution is due to constraint force on that atom

    if (EVFLAG) {
      double massone,vr[6];

      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];

      const double fc0 = 0.5*(massone*(v[i].x - vx)/dtf - f[i].x);
      const double fc1 = 0.5*(massone*(v[i].y - vy)/dtf - f[i].y);
      const double fc2 = 0.5*(massone*(v[i].z - vz)/dtf - f[i].z);

      vr[0] = x0*fc0; vr[1] = x1*fc1; vr[2] = x2*fc2;
      vr[3] = x0*fc1; vr[4] = x0*fc2; vr[5] = x1*fc2;

      // Fix::v_tally() is not thread safe, so we do this manually here
      // accumulate global virial into thread-local variables for reduction
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
  }

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

  // set orientation, omega, angmom of each extended particle
  // XXX: extended particle info not yet multi-threaded

  if (extended) {
    double *shape,*quatatom,*inertiaatom;
    double theta_body,theta;
    double ione[3],exone[3],eyone[3],ezone[3],p[3][3];

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
template <int TRICLINIC, int EVFLAG>
void FixRigidOMP::set_v_thr()
{
  dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const v = (dbl3_t *) atom->v[0];
  const dbl3_t * _noalias const f = (dbl3_t *) atom->f[0];
  const double * _noalias const rmass = atom->rmass;
  const double * _noalias const mass = atom->mass;
  const int * _noalias const type = atom->type;

  const double xprd = domain->xprd;
  const double yprd = domain->yprd;
  const double zprd = domain->zprd;
  const double xy = domain->xy;
  const double xz = domain->xz;
  const double yz = domain->yz;

  double v0=0.0,v1=0.0,v2=0.0,v3=0.0,v4=0.0,v5=0.0;

  // set v of each atom

  const int nlocal = atom->nlocal;

#if defined(_OPENMP)
#pragma omp parallel for default(none) reduction(+:v0,v1,v2,v3,v4,v5)
#endif
  for (int i = 0; i < nlocal; i++) {
    const int ibody = body[i];
    if (ibody < 0) continue;

    const dbl3_t &vcmi = * ((dbl3_t *) vcm[ibody]);
    const dbl3_t &omegai = * ((dbl3_t *) omega[ibody]);
    double delta[3],vx,vy,vz;

    MathExtra::matvec(ex_space[ibody],ey_space[ibody],
                      ez_space[ibody],displace[i],delta);

    // save old velocities for virial

    if (EVFLAG) {
      vx = v[i].x;
      vy = v[i].y;
      vz = v[i].z;
    }

    v[i].x = omegai.y*delta[2] - omegai.z*delta[1] + vcmi.x;
    v[i].y = omegai.z*delta[0] - omegai.x*delta[2] + vcmi.y;
    v[i].z = omegai.x*delta[1] - omegai.y*delta[0] + vcmi.z;

    // virial = unwrapped coords dotted into body constraint force
    // body constraint force = implied force due to v change minus f external
    // assume f does not include forces internal to body
    // 1/2 factor b/c initial_integrate contributes other half
    // assume per-atom contribution is due to constraint force on that atom

    if (EVFLAG) {
      double massone, vr[6];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];

      const int xbox = (xcmimage[i] & IMGMASK) - IMGMAX;
      const int ybox = (xcmimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
      const int zbox = (xcmimage[i] >> IMG2BITS) - IMGMAX;
      const double deltax = xbox*xprd + (TRICLINIC ? ybox*xy + zbox*xz : 0.0);
      const double deltay = ybox*yprd + (TRICLINIC ? zbox*yz : 0.0);
      const double deltaz = zbox*zprd;

      const double fc0 = 0.5*(massone*(v[i].x - vx)/dtf - f[i].x);
      const double fc1 = 0.5*(massone*(v[i].y - vy)/dtf - f[i].y);
      const double fc2 = 0.5*(massone*(v[i].z - vz)/dtf - f[i].z);

      const double x0 = x[i].x + deltax;
      const double x1 = x[i].y + deltay;
      const double x2 = x[i].z + deltaz;

      vr[0] = x0*fc0; vr[1] = x1*fc1; vr[2] = x2*fc2;
      vr[3] = x0*fc1; vr[4] = x0*fc2; vr[5] = x1*fc2;

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

    AtomVecEllipsoid::Bonus *ebonus;
    if (avec_ellipsoid) ebonus = avec_ellipsoid->bonus;
    AtomVecTri::Bonus *tbonus;
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
