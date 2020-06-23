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
   Contributing author: Tony Sheh (U Michigan), Trung Dac Nguyen (U Michigan)
   references: Kamberaj et al., J. Chem. Phys. 122, 224114 (2005)
               Miller et al., J Chem Phys. 116, 8649-8659 (2002)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:
      Morteza Jalalvand (IASBS)  jalalvand.m AT gmail.com

    This is an extension of fix/rigid/nve to SPH/SDPD particles
    You can see the original copyright notice of fix/rigid authors above
    Note that the Kamberaj paper was related to the nvt variant
    and all codes relevant to that has been removed
------------------------------------------------------------------------- */

#include "fix_rigid_meso.h"
#include "math_extra.h"
#include "atom.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixRigidMeso::FixRigidMeso (LAMMPS *lmp, int narg, char **arg) :
FixRigid (lmp, narg, arg) {
  scalar_flag = 0;
  size_array_cols = 28;
  if ((atom->esph_flag != 1) || (atom->rho_flag != 1))
    error->all (FLERR, "fix rigid/meso command requires atom_style with"
                " both energy and density");

  if (langflag || tstat_flag)
    error->all (FLERR,"Can not use thermostat with fix rigid/meso");

  if (pstat_flag)
    error->all (FLERR,"Can not use barostat with fix rigid/meso");

  // memory allocation and initialization

  memory->create(conjqm,nbody,4,"rigid_nh:conjqm");
}

/* ---------------------------------------------------------------------- */

FixRigidMeso::~FixRigidMeso () {
  memory->destroy(conjqm);
}

/* ---------------------------------------------------------------------- */

int FixRigidMeso::setmask () {
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= PRE_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRigidMeso::setup (int vflag) {
  FixRigid::setup(vflag);

  double mbody[3];
  for (int ibody = 0; ibody < nbody; ibody++) {
    MathExtra::transpose_matvec (ex_space[ibody],ey_space[ibody],ez_space[ibody],
                                 angmom[ibody],mbody);
    MathExtra::quatvec (quat[ibody],mbody,conjqm[ibody]);
    conjqm[ibody][0] *= 2.0;
    conjqm[ibody][1] *= 2.0;
    conjqm[ibody][2] *= 2.0;
    conjqm[ibody][3] *= 2.0;
  }
}

/* ----------------------------------------------------------------------
   perform preforce velocity Verlet integration
   see Kamberaj paper for step references
------------------------------------------------------------------------- */

void FixRigidMeso::initial_integrate (int vflag) {
  double dtfm,mbody[3],tbody[3],fquat[4];
  double dtf2 = dtf * 2.0;

  // update xcm, vcm, quat, conjqm and angmom

  for (int ibody = 0; ibody < nbody; ibody++) {

    // step 1.1 - update vcm by 1/2 step

    dtfm = dtf / masstotal[ibody];
    vcm[ibody][0] += dtfm * fcm[ibody][0] * fflag[ibody][0];
    vcm[ibody][1] += dtfm * fcm[ibody][1] * fflag[ibody][1];
    vcm[ibody][2] += dtfm * fcm[ibody][2] * fflag[ibody][2];

    // step 1.2 - update xcm by full step

    xcm[ibody][0] += dtv * vcm[ibody][0];
    xcm[ibody][1] += dtv * vcm[ibody][1];
    xcm[ibody][2] += dtv * vcm[ibody][2];

    // step 1.3 - apply torque (body coords) to quaternion momentum

    torque[ibody][0] *= tflag[ibody][0];
    torque[ibody][1] *= tflag[ibody][1];
    torque[ibody][2] *= tflag[ibody][2];

    MathExtra::transpose_matvec (ex_space[ibody],ey_space[ibody],ez_space[ibody],
                                 torque[ibody],tbody);
    MathExtra::quatvec (quat[ibody],tbody,fquat);

    conjqm[ibody][0] += dtf2 * fquat[0];
    conjqm[ibody][1] += dtf2 * fquat[1];
    conjqm[ibody][2] += dtf2 * fquat[2];
    conjqm[ibody][3] += dtf2 * fquat[3];

    // step 1.4 to 1.13 - use no_squish rotate to update p and q

    MathExtra::no_squish_rotate (3,conjqm[ibody],quat[ibody],inertia[ibody],dtq);
    MathExtra::no_squish_rotate (2,conjqm[ibody],quat[ibody],inertia[ibody],dtq);
    MathExtra::no_squish_rotate (1,conjqm[ibody],quat[ibody],inertia[ibody],dtv);
    MathExtra::no_squish_rotate (2,conjqm[ibody],quat[ibody],inertia[ibody],dtq);
    MathExtra::no_squish_rotate (3,conjqm[ibody],quat[ibody],inertia[ibody],dtq);

    // update exyz_space
    // transform p back to angmom
    // update angular velocity

    MathExtra::q_to_exyz (quat[ibody],ex_space[ibody],ey_space[ibody],
                          ez_space[ibody]);
    MathExtra::invquatvec (quat[ibody],conjqm[ibody],mbody);
    MathExtra::matvec (ex_space[ibody],ey_space[ibody],ez_space[ibody],
                       mbody,angmom[ibody]);

    angmom[ibody][0] *= 0.5;
    angmom[ibody][1] *= 0.5;
    angmom[ibody][2] *= 0.5;

    MathExtra::angmom_to_omega (angmom[ibody],ex_space[ibody],ey_space[ibody],
                                ez_space[ibody],inertia[ibody],omega[ibody]);
  }

  // virial setup before call to set_xv

  if (vflag) v_setup(vflag);
  else evflag = 0;

  // set coords/orient and velocity/rotation of atoms in rigid bodies
  // from quarternion and omega

  set_xv();
}

/* ---------------------------------------------------------------------- */

void FixRigidMeso::final_integrate () {
  int ibody;
  double dtfm;
  double mbody[3],tbody[3],fquat[4];

  double dtf2 = dtf * 2.0;

  // late calculation of forces and torques (if requested)

  if (!earlyflag) compute_forces_and_torques();

  // update vcm and angmom
  // fflag,tflag = 0 for some dimensions in 2d

  for (ibody = 0; ibody < nbody; ibody++) {

    // update vcm by 1/2 step

    dtfm = dtf / masstotal[ibody];

    vcm[ibody][0] += dtfm * fcm[ibody][0] * fflag[ibody][0];
    vcm[ibody][1] += dtfm * fcm[ibody][1] * fflag[ibody][1];
    vcm[ibody][2] += dtfm * fcm[ibody][2] * fflag[ibody][2];

    // update conjqm, then transform to angmom, set velocity again
    // virial is already setup from initial_integrate

    torque[ibody][0] *= tflag[ibody][0];
    torque[ibody][1] *= tflag[ibody][1];
    torque[ibody][2] *= tflag[ibody][2];

    MathExtra::transpose_matvec (ex_space[ibody],ey_space[ibody],
                                 ez_space[ibody],torque[ibody],tbody);
    MathExtra::quatvec (quat[ibody],tbody,fquat);

    conjqm[ibody][0] += dtf2 * fquat[0];
    conjqm[ibody][1] += dtf2 * fquat[1];
    conjqm[ibody][2] += dtf2 * fquat[2];
    conjqm[ibody][3] += dtf2 * fquat[3];

    MathExtra::invquatvec (quat[ibody],conjqm[ibody],mbody);
    MathExtra::matvec (ex_space[ibody],ey_space[ibody],ez_space[ibody],
                       mbody,angmom[ibody]);

    angmom[ibody][0] *= 0.5;
    angmom[ibody][1] *= 0.5;
    angmom[ibody][2] *= 0.5;

    MathExtra::angmom_to_omega (angmom[ibody],ex_space[ibody],ey_space[ibody],
                                ez_space[ibody],inertia[ibody],omega[ibody]);
  }

  // set velocity/rotation of atoms in rigid bodies
  // virial is already setup from initial_integrate

  set_v();
}

/* ----------------------------------------------------------------------
   set space-frame coords and velocity of each atom in each rigid body
   set orientation and rotation of extended particles
   x = Q displace + Xcm, mapped back to periodic box
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

void FixRigidMeso::set_xv () {
  int ibody;
  int xbox,ybox,zbox;
  double x0,x1,x2,v0,v1,v2,fc0,fc1,fc2,massone;
  double xy,xz,yz;
  double vr[6];

  double **x = atom->x;
  double **v = atom->v;
  double **vest = atom->vest;
  double **f = atom->f;
  double *esph = atom->esph;
  double *desph = atom->desph;
  double *rho = atom->rho;
  double *drho = atom->drho;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  if (triclinic) {
    xy = domain->xy;
    xz = domain->xz;
    yz = domain->yz;
  }

  // set x and v of each atom

  for (int i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;

    // half-step update of particle internal energy and density
    esph[i] += dtf * desph[i];
    rho[i] += dtf * drho[i];

    ibody = body[i];

    xbox = (xcmimage[i] & IMGMASK) - IMGMAX;
    ybox = (xcmimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (xcmimage[i] >> IMG2BITS) - IMGMAX;

    // save old positions and velocities for virial

    if (evflag) {
      if (triclinic == 0) {
        x0 = x[i][0] + xbox*xprd;
        x1 = x[i][1] + ybox*yprd;
        x2 = x[i][2] + zbox*zprd;
      } else {
        x0 = x[i][0] + xbox*xprd + ybox*xy + zbox*xz;
        x1 = x[i][1] + ybox*yprd + zbox*yz;
        x2 = x[i][2] + zbox*zprd;
      }
    }

    v0 = v[i][0];
    v1 = v[i][1];
    v2 = v[i][2];

    // x = displacement from center-of-mass, based on body orientation
    // v = vcm + omega around center-of-mass
    // vest = 2*v - v_old

    MathExtra::matvec (ex_space[ibody],ey_space[ibody],
                       ez_space[ibody],displace[i],x[i]);

    v[i][0] = omega[ibody][1]*x[i][2] - omega[ibody][2]*x[i][1] +
      vcm[ibody][0];
    v[i][1] = omega[ibody][2]*x[i][0] - omega[ibody][0]*x[i][2] +
      vcm[ibody][1];
    v[i][2] = omega[ibody][0]*x[i][1] - omega[ibody][1]*x[i][0] +
      vcm[ibody][2];

    vest[i][0] = 2*v[i][0] - v0;
    vest[i][1] = 2*v[i][1] - v1;
    vest[i][2] = 2*v[i][2] - v2;

    // add center of mass to displacement
    // map back into periodic box via xbox,ybox,zbox
    // for triclinic, add in box tilt factors as well

    if (triclinic == 0) {
      x[i][0] += xcm[ibody][0] - xbox*xprd;
      x[i][1] += xcm[ibody][1] - ybox*yprd;
      x[i][2] += xcm[ibody][2] - zbox*zprd;
    } else {
      x[i][0] += xcm[ibody][0] - xbox*xprd - ybox*xy - zbox*xz;
      x[i][1] += xcm[ibody][1] - ybox*yprd - zbox*yz;
      x[i][2] += xcm[ibody][2] - zbox*zprd;
    }

    // virial = unwrapped coords dotted into body constraint force
    // body constraint force = implied force due to v change minus f external
    // assume f does not include forces internal to body
    // 1/2 factor b/c final_integrate contributes other half
    // assume per-atom contribution is due to constraint force on that atom

    if (evflag) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      fc0 = massone*(v[i][0] - v0)/dtf - f[i][0];
      fc1 = massone*(v[i][1] - v1)/dtf - f[i][1];
      fc2 = massone*(v[i][2] - v2)/dtf - f[i][2];

      vr[0] = 0.5*x0*fc0;
      vr[1] = 0.5*x1*fc1;
      vr[2] = 0.5*x2*fc2;
      vr[3] = 0.5*x0*fc1;
      vr[4] = 0.5*x0*fc2;
      vr[5] = 0.5*x1*fc2;

      v_tally(1,&i,1.0,vr);
    }
  }

  // set orientation, omega, angmom of each extended particle

  if (extended) {
    // TBD
  }
}

/* ----------------------------------------------------------------------
   set space-frame velocity of each atom in a rigid body
   set omega and angmom of extended particles
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

void FixRigidMeso::set_v () {
  int xbox,ybox,zbox;
  double x0,x1,x2,v0,v1,v2,fc0,fc1,fc2,massone;
  double xy,xz,yz;
  double delta[3],vr[6];

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *esph = atom->esph;
  double *desph = atom->desph;
  double *rho = atom->rho;
  double *drho = atom->drho;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  if (triclinic) {
    xy = domain->xy;
    xz = domain->xz;
    yz = domain->yz;
  }

  // set v of each atom

  for (int i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;

    // half-step update of particle internal energy and density
    esph[i] += dtf * desph[i];
    rho[i] += dtf * drho[i];

    const int ibody = body[i];

    MathExtra::matvec (ex_space[ibody],ey_space[ibody],
                       ez_space[ibody],displace[i],delta);

    // save old velocities for virial

    if (evflag) {
      v0 = v[i][0];
      v1 = v[i][1];
      v2 = v[i][2];
    }

    v[i][0] = omega[ibody][1]*delta[2] - omega[ibody][2]*delta[1] +
      vcm[ibody][0];
    v[i][1] = omega[ibody][2]*delta[0] - omega[ibody][0]*delta[2] +
      vcm[ibody][1];
    v[i][2] = omega[ibody][0]*delta[1] - omega[ibody][1]*delta[0] +
      vcm[ibody][2];

    // virial = unwrapped coords dotted into body constraint force
    // body constraint force = implied force due to v change minus f external
    // assume f does not include forces internal to body
    // 1/2 factor b/c initial_integrate contributes other half
    // assume per-atom contribution is due to constraint force on that atom

    if (evflag) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      fc0 = massone*(v[i][0] - v0)/dtf - f[i][0];
      fc1 = massone*(v[i][1] - v1)/dtf - f[i][1];
      fc2 = massone*(v[i][2] - v2)/dtf - f[i][2];

      xbox = (xcmimage[i] & IMGMASK) - IMGMAX;
      ybox = (xcmimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
      zbox = (xcmimage[i] >> IMG2BITS) - IMGMAX;

      if (triclinic == 0) {
        x0 = x[i][0] + xbox*xprd;
        x1 = x[i][1] + ybox*yprd;
        x2 = x[i][2] + zbox*zprd;
      } else {
        x0 = x[i][0] + xbox*xprd + ybox*xy + zbox*xz;
        x1 = x[i][1] + ybox*yprd + zbox*yz;
        x2 = x[i][2] + zbox*zprd;
      }

      vr[0] = 0.5*x0*fc0;
      vr[1] = 0.5*x1*fc1;
      vr[2] = 0.5*x2*fc2;
      vr[3] = 0.5*x0*fc1;
      vr[4] = 0.5*x0*fc2;
      vr[5] = 0.5*x1*fc2;

      v_tally(1,&i,1.0,vr);
    }
  }

  // set omega, angmom of each extended particle

  if (extended) {
    // TBD
  }
}

/* ----------------------------------------------------------------------
   return attributes of a rigid body
   19 values per body
   xcm = 0,1,2; vcm = 3,4,5; fcm = 6,7,8;
   quat = 9,10,11,12; omega = 13,14,15; torque = 16,17,18;
   inertia = 19,20,21; angmom = 22,23,24;
   image = 25,26,27
------------------------------------------------------------------------- */

double FixRigidMeso::compute_array (int i, int j) {
  if (j < 3) return xcm[i][j];
  if (j < 6) return vcm[i][j-3];
  if (j < 9) return fcm[i][j-6];
  if (j < 13) return quat[i][j-9];
  if (j < 16) return omega[i][j-13];
  if (j < 19) return torque[i][j-16];
  if (j < 22) return inertia[i][j-19];
  if (j < 25) return angmom[i][j-22];
  if (j == 25) return (imagebody[i] & IMGMASK) - IMGMAX;
  if (j == 26) return (imagebody[i] >> IMGBITS & IMGMASK) - IMGMAX;
  return (imagebody[i] >> IMG2BITS) - IMGMAX;
}
