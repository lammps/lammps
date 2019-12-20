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

#include "fix_rigid_small_omp.h"
#include <cmath>
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "comm.h"
#include "domain.h"
#include "timer.h"

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

void FixRigidSmallOMP::initial_integrate(int vflag)
{

#if defined(_OPENMP)
#pragma omp parallel for default(none) schedule(static)
#endif
  for (int ibody = 0; ibody < nlocal_body; ibody++) {

    Body &b = body[ibody];

    // update vcm by 1/2 step

    const double dtfm = dtf / b.mass;
    b.vcm[0] += dtfm * b.fcm[0];
    b.vcm[1] += dtfm * b.fcm[1];
    b.vcm[2] += dtfm * b.fcm[2];

    // update xcm by full step

    b.xcm[0] += dtv * b.vcm[0];
    b.xcm[1] += dtv * b.vcm[1];
    b.xcm[2] += dtv * b.vcm[2];

    // update angular momentum by 1/2 step

    b.angmom[0] += dtf * b.torque[0];
    b.angmom[1] += dtf * b.torque[1];
    b.angmom[2] += dtf * b.torque[2];

    // compute omega at 1/2 step from angmom at 1/2 step and current q
    // update quaternion a full step via Richardson iteration
    // returns new normalized quaternion, also updated omega at 1/2 step
    // update ex,ey,ez to reflect new quaternion

    MathExtra::angmom_to_omega(b.angmom,b.ex_space,b.ey_space,
                               b.ez_space,b.inertia,b.omega);
    MathExtra::richardson(b.quat,b.angmom,b.omega,b.inertia,dtq);
    MathExtra::q_to_exyz(b.quat,b.ex_space,b.ey_space,b.ez_space);
  } // end of omp parallel for

  // virial setup before call to set_xv

  if (vflag) v_setup(vflag);
  else evflag = 0;

  // forward communicate updated info of all bodies

  commflag = INITIAL;
  comm->forward_comm_fix(this,26);

  // set coords/orient and velocity/rotation of atoms in rigid bodies

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

void FixRigidSmallOMP::compute_forces_and_torques()
{
  double * const * _noalias const x = atom->x;
  const dbl3_t * _noalias const f = (dbl3_t *) atom->f[0];
  const double * const * const torque_one = atom->torque;
  const int nlocal = atom->nlocal;
  const int nthreads=comm->nthreads;

#if defined(_OPENMP)
#pragma omp parallel for default(none) schedule(static)
#endif
  for (int ibody = 0; ibody < nlocal_body+nghost_body; ibody++) {
    double * _noalias const fcm = body[ibody].fcm;
    fcm[0] = fcm[1] = fcm[2] = 0.0;
    double * _noalias const tcm = body[ibody].torque;
    tcm[0] = tcm[1] = tcm[2] = 0.0;
  }

  // sum over atoms to get force and torque on rigid body
  // we likely have a large number of rigid objects with only a
  // a few atoms each. so we loop over all atoms for all threads
  // and then each thread only processes some bodies.

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
      int ibody = atom2body[i];
      if ((ibody < 0) || (ibody % nthreads != tid)) continue;

      Body &b = body[ibody];

      double unwrap[3];
      domain->unmap(x[i],xcmimage[i],unwrap);

      double * _noalias const fcm = b.fcm;
      double * _noalias const xcm = b.xcm;
      double * _noalias const tcm = b.torque;

      const double dx = unwrap[0] - xcm[0];
      const double dy = unwrap[1] - xcm[1];
      const double dz = unwrap[2] - xcm[2];

      fcm[0] += f[i].x;
      fcm[1] += f[i].y;
      fcm[2] += f[i].z;

      tcm[0] += dy*f[i].z - dz*f[i].y;
      tcm[1] += dz*f[i].x - dx*f[i].z;
      tcm[2] += dx*f[i].y - dy*f[i].x;

      if (extended && (eflags[i] & TORQUE)) {
        tcm[0] += torque_one[i][0];
        tcm[1] += torque_one[i][1];
        tcm[2] += torque_one[i][2];
      }
    }
  } // end of omp parallel region

  // reverse communicate fcm, torque of all bodies

  commflag = FORCE_TORQUE;
  comm->reverse_comm_fix(this,6);

  // include Langevin thermostat forces and torques

  if (langflag) {
#if defined(_OPENMP)
#pragma omp parallel for default(none) schedule(static)
#endif
    for (int ibody = 0; ibody < nlocal_body; ibody++) {
      double * _noalias const fcm = body[ibody].fcm;
      fcm[0] += langextra[ibody][0];
      fcm[1] += langextra[ibody][1];
      fcm[2] += langextra[ibody][2];
      double * _noalias const tcm = body[ibody].torque;
      tcm[0] += langextra[ibody][3];
      tcm[1] += langextra[ibody][4];
      tcm[2] += langextra[ibody][5];
    }
  }

  // add gravity force to COM of each body

  if (id_gravity) {
#if defined(_OPENMP)
#pragma omp parallel for default(none) private(ibody) schedule(static)
#endif
    for (int ibody = 0; ibody < nbody; ibody++) {
      double * _noalias const fcm = body[ibody].fcm;
      const double mass = body[ibody].mass;
      fcm[0] += gvec[0]*mass;
      fcm[1] += gvec[1]*mass;
      fcm[2] += gvec[2]*mass;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixRigidSmallOMP::final_integrate()
{
  if (!earlyflag) compute_forces_and_torques();

  // update vcm and angmom, recompute omega

#if defined(_OPENMP)
#pragma omp parallel for default(none) schedule(static)
#endif
  for (int ibody = 0; ibody < nlocal_body; ibody++) {
    Body &b = body[ibody];

    // update vcm by 1/2 step

    const double dtfm = dtf / b.mass;
    b.vcm[0] += dtfm * b.fcm[0];
    b.vcm[1] += dtfm * b.fcm[1];
    b.vcm[2] += dtfm * b.fcm[2];

    // update angular momentum by 1/2 step

    b.angmom[0] += dtf * b.torque[0];
    b.angmom[1] += dtf * b.torque[1];
    b.angmom[2] += dtf * b.torque[2];

    MathExtra::angmom_to_omega(b.angmom,b.ex_space,b.ey_space,
                               b.ez_space,b.inertia,b.omega);
  }

  // forward communicate updated info of all bodies

  commflag = FINAL;
  comm->forward_comm_fix(this,10);

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
------------------------------------------------------------------------- */

template <int TRICLINIC, int EVFLAG>
void FixRigidSmallOMP::set_xv_thr()
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
    const int ibody = atom2body[i];
    if (ibody  < 0) continue;

    Body &b = body[ibody];

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

    MathExtra::matvec(b.ex_space,b.ey_space,b.ez_space,displace[i],&x[i].x);

    v[i].x = b.omega[1]*x[i].z - b.omega[2]*x[i].y + b.vcm[0];
    v[i].y = b.omega[2]*x[i].x - b.omega[0]*x[i].z + b.vcm[1];
    v[i].z = b.omega[0]*x[i].y - b.omega[1]*x[i].x + b.vcm[2];

    // add center of mass to displacement
    // map back into periodic box via xbox,ybox,zbox
    // for triclinic, add in box tilt factors as well

    x[i].x += b.xcm[0] - deltax;
    x[i].y += b.xcm[1] - deltay;
    x[i].z += b.xcm[2] - deltaz;

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
    double ione[3],exone[3],eyone[3],ezone[3],p[3][3];
    double theta_body,theta;
    double *shape,*quatatom,*inertiaatom;

    AtomVecEllipsoid::Bonus *ebonus;
    if (avec_ellipsoid) ebonus = avec_ellipsoid->bonus;
    AtomVecLine::Bonus *lbonus;
    if (avec_line) lbonus = avec_line->bonus;
    AtomVecTri::Bonus *tbonus;
    if (avec_tri) tbonus = avec_tri->bonus;
    double **omega = atom->omega;
    double **angmom = atom->angmom;
    double **mu = atom->mu;
    int *ellipsoid = atom->ellipsoid;
    int *line = atom->line;
    int *tri = atom->tri;

    for (int i = 0; i < nlocal; i++) {
      if (atom2body[i] < 0) continue;
      Body &b = body[atom2body[i]];

      if (eflags[i] & SPHERE) {
        omega[i][0] = b.omega[0];
        omega[i][1] = b.omega[1];
        omega[i][2] = b.omega[2];
      } else if (eflags[i] & ELLIPSOID) {
        shape = ebonus[ellipsoid[i]].shape;
        quatatom = ebonus[ellipsoid[i]].quat;
        MathExtra::quatquat(b.quat,orient[i],quatatom);
        MathExtra::qnormalize(quatatom);
        ione[0] = EINERTIA*rmass[i] * (shape[1]*shape[1] + shape[2]*shape[2]);
        ione[1] = EINERTIA*rmass[i] * (shape[0]*shape[0] + shape[2]*shape[2]);
        ione[2] = EINERTIA*rmass[i] * (shape[0]*shape[0] + shape[1]*shape[1]);
        MathExtra::q_to_exyz(quatatom,exone,eyone,ezone);
        MathExtra::omega_to_angmom(b.omega,exone,eyone,ezone,ione,angmom[i]);
      } else if (eflags[i] & LINE) {
        if (b.quat[3] >= 0.0) theta_body = 2.0*acos(b.quat[0]);
        else theta_body = -2.0*acos(b.quat[0]);
        theta = orient[i][0] + theta_body;
        while (theta <= -MY_PI) theta += MY_2PI;
        while (theta > MY_PI) theta -= MY_2PI;
        lbonus[line[i]].theta = theta;
        omega[i][0] = b.omega[0];
        omega[i][1] = b.omega[1];
        omega[i][2] = b.omega[2];
      } else if (eflags[i] & TRIANGLE) {
        inertiaatom = tbonus[tri[i]].inertia;
        quatatom = tbonus[tri[i]].quat;
        MathExtra::quatquat(b.quat,orient[i],quatatom);
        MathExtra::qnormalize(quatatom);
        MathExtra::q_to_exyz(quatatom,exone,eyone,ezone);
        MathExtra::omega_to_angmom(b.omega,exone,eyone,ezone,
                                   inertiaatom,angmom[i]);
      }
      if (eflags[i] & DIPOLE) {
        MathExtra::quat_to_mat(b.quat,p);
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
------------------------------------------------------------------------- */

template <int TRICLINIC, int EVFLAG>
void FixRigidSmallOMP::set_v_thr()
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
    const int ibody = atom2body[i];
    if (ibody < 0) continue;

    Body &b = body[atom2body[i]];
    double delta[3],vx,vy,vz;

    MathExtra::matvec(b.ex_space,b.ey_space,b.ez_space,displace[i],delta);

    // save old velocities for virial

    if (EVFLAG) {
      vx = v[i].x;
      vy = v[i].y;
      vz = v[i].z;
    }

    v[i].x = b.omega[1]*delta[2] - b.omega[2]*delta[1] + b.vcm[0];
    v[i].y = b.omega[2]*delta[0] - b.omega[0]*delta[2] + b.vcm[1];
    v[i].z = b.omega[0]*delta[1] - b.omega[1]*delta[0] + b.vcm[2];

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
    double ione[3],exone[3],eyone[3],ezone[3];
    double *shape,*quatatom,*inertiaatom;

    AtomVecEllipsoid::Bonus *ebonus;
    if (avec_ellipsoid) ebonus = avec_ellipsoid->bonus;
    AtomVecTri::Bonus *tbonus;
    if (avec_tri) tbonus = avec_tri->bonus;
    double **omega = atom->omega;
    double **angmom = atom->angmom;
    int *ellipsoid = atom->ellipsoid;
    int *tri = atom->tri;

    for (int i = 0; i < nlocal; i++) {
      if (atom2body[i] < 0) continue;
      Body &b = body[atom2body[i]];

      if (eflags[i] & SPHERE) {
        omega[i][0] = b.omega[0];
        omega[i][1] = b.omega[1];
        omega[i][2] = b.omega[2];
      } else if (eflags[i] & ELLIPSOID) {
        shape = ebonus[ellipsoid[i]].shape;
        quatatom = ebonus[ellipsoid[i]].quat;
        ione[0] = EINERTIA*rmass[i] * (shape[1]*shape[1] + shape[2]*shape[2]);
        ione[1] = EINERTIA*rmass[i] * (shape[0]*shape[0] + shape[2]*shape[2]);
        ione[2] = EINERTIA*rmass[i] * (shape[0]*shape[0] + shape[1]*shape[1]);
        MathExtra::q_to_exyz(quatatom,exone,eyone,ezone);
        MathExtra::omega_to_angmom(b.omega,exone,eyone,ezone,ione,
                                   angmom[i]);
      } else if (eflags[i] & LINE) {
        omega[i][0] = b.omega[0];
        omega[i][1] = b.omega[1];
        omega[i][2] = b.omega[2];
      } else if (eflags[i] & TRIANGLE) {
        inertiaatom = tbonus[tri[i]].inertia;
        quatatom = tbonus[tri[i]].quat;
        MathExtra::q_to_exyz(quatatom,exone,eyone,ezone);
        MathExtra::omega_to_angmom(b.omega,exone,eyone,ezone,
                                   inertiaatom,angmom[i]);
      }
    }
  }
}

