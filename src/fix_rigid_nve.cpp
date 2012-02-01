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
               Miller et al., J Chem Phys. 116, 8649-8659 
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "string.h"
#include "fix_rigid_nve.h"
#include "math_extra.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "modify.h"
#include "group.h"
#include "comm.h"
#include "force.h"
#include "output.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixRigidNVE::FixRigidNVE(LAMMPS *lmp, int narg, char **arg) :
  FixRigid(lmp, narg, arg)
{ 
  memory->create(conjqm,nbody,4,"rigid/nve:conjqm");
}

/* ---------------------------------------------------------------------- */

FixRigidNVE::~FixRigidNVE()
{
  memory->destroy(conjqm);
}

/* ---------------------------------------------------------------------- */

void FixRigidNVE::setup(int vflag)
{
  FixRigid::setup(vflag);
  
  double mbody[3];
  for (int ibody = 0; ibody < nbody; ibody++) {
    MathExtra::transpose_matvec(ex_space[ibody],ey_space[ibody],ez_space[ibody],
				angmom[ibody],mbody);
    MathExtra::quatvec(quat[ibody],mbody,conjqm[ibody]);
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

void FixRigidNVE::initial_integrate(int vflag)
{
  double dtfm,mbody[3],tbody[3],fquat[4];
  double dtf2 = dtf * 2.0;
  
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
    
    MathExtra::transpose_matvec(ex_space[ibody],ey_space[ibody],ez_space[ibody],
				torque[ibody],tbody);
    MathExtra::quatvec(quat[ibody],tbody,fquat);
    
    conjqm[ibody][0] += dtf2 * fquat[0];
    conjqm[ibody][1] += dtf2 * fquat[1];
    conjqm[ibody][2] += dtf2 * fquat[2];
    conjqm[ibody][3] += dtf2 * fquat[3];
  
    // step 1.4 to 1.13 - use no_squish rotate to update p and q
  
    no_squish_rotate(3,conjqm[ibody],quat[ibody],inertia[ibody],dtq);
    no_squish_rotate(2,conjqm[ibody],quat[ibody],inertia[ibody],dtq);
    no_squish_rotate(1,conjqm[ibody],quat[ibody],inertia[ibody],dtv);
    no_squish_rotate(2,conjqm[ibody],quat[ibody],inertia[ibody],dtq);
    no_squish_rotate(3,conjqm[ibody],quat[ibody],inertia[ibody],dtq);
  
    // update exyz_space
    // transform p back to angmom
    // update angular velocity
    
    MathExtra::q_to_exyz(quat[ibody],ex_space[ibody],ey_space[ibody],
			 ez_space[ibody]);
    MathExtra::invquatvec(quat[ibody],conjqm[ibody],mbody);
    MathExtra::matvec(ex_space[ibody],ey_space[ibody],ez_space[ibody],
		      mbody,angmom[ibody]);
    
    angmom[ibody][0] *= 0.5;
    angmom[ibody][1] *= 0.5;
    angmom[ibody][2] *= 0.5;
    
    MathExtra::angmom_to_omega(angmom[ibody],ex_space[ibody],ey_space[ibody],
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

void FixRigidNVE::final_integrate()
{
  int i,ibody;
  double dtfm,xy,xz,yz;

  // sum over atoms to get force and torque on rigid body
  
  int *image = atom->image;
  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  if (triclinic) {
    xy = domain->xy;
    xz = domain->xz;
    yz = domain->yz;
  }

  int xbox,ybox,zbox;
  double xunwrap,yunwrap,zunwrap,dx,dy,dz;
  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;
  
  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    sum[ibody][0] += f[i][0];
    sum[ibody][1] += f[i][1];
    sum[ibody][2] += f[i][2];
      
    xbox = (image[i] & 1023) - 512;
    ybox = (image[i] >> 10 & 1023) - 512;
    zbox = (image[i] >> 20) - 512;

    if (triclinic == 0) {
      xunwrap = x[i][0] + xbox*xprd;
      yunwrap = x[i][1] + ybox*yprd;
      zunwrap = x[i][2] + zbox*zprd;
    } else {
      xunwrap = x[i][0] + xbox*xprd + ybox*xy + zbox*xz;
      yunwrap = x[i][1] + ybox*yprd + zbox*yz;
      zunwrap = x[i][2] + zbox*zprd;
    }

    dx = xunwrap - xcm[ibody][0];
    dy = yunwrap - xcm[ibody][1];
    dz = zunwrap - xcm[ibody][2];
    
    sum[ibody][3] += dy*f[i][2] - dz*f[i][1];
    sum[ibody][4] += dz*f[i][0] - dx*f[i][2];
    sum[ibody][5] += dx*f[i][1] - dy*f[i][0];
  }
  
  // extended particles add their torque to torque of body

  if (extended) {
    double **torque_one = atom->torque;

    for (i = 0; i < nlocal; i++) {
      if (body[i] < 0) continue;
      ibody = body[i];

      if (eflags[i] & TORQUE) {
	sum[ibody][3] += torque_one[i][0];
	sum[ibody][4] += torque_one[i][1];
	sum[ibody][5] += torque_one[i][2];
      }
    }
  }

  MPI_Allreduce(sum[0],all[0],6*nbody,MPI_DOUBLE,MPI_SUM,world);
  
  // update vcm and angmom
  // include Langevin thermostat forces
  // fflag,tflag = 0 for some dimensions in 2d

  double mbody[3],tbody[3],fquat[4];
  double dtf2 = dtf * 2.0;

  for (ibody = 0; ibody < nbody; ibody++) {
    fcm[ibody][0] = all[ibody][0] + langextra[ibody][0];
    fcm[ibody][1] = all[ibody][1] + langextra[ibody][1];
    fcm[ibody][2] = all[ibody][2] + langextra[ibody][2];
    torque[ibody][0] = all[ibody][3] + langextra[ibody][3];
    torque[ibody][1] = all[ibody][4] + langextra[ibody][4];
    torque[ibody][2] = all[ibody][5] + langextra[ibody][5];

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
    
    MathExtra::transpose_matvec(ex_space[ibody],ey_space[ibody],ez_space[ibody],
				torque[ibody],tbody);
    MathExtra::quatvec(quat[ibody],tbody,fquat);
    
    conjqm[ibody][0] += dtf2 * fquat[0];
    conjqm[ibody][1] += dtf2 * fquat[1];
    conjqm[ibody][2] += dtf2 * fquat[2];
    conjqm[ibody][3] += dtf2 * fquat[3];
    
    MathExtra::invquatvec(quat[ibody],conjqm[ibody],mbody);
    MathExtra::matvec(ex_space[ibody],ey_space[ibody],ez_space[ibody],
		      mbody,angmom[ibody]);
    
    angmom[ibody][0] *= 0.5;
    angmom[ibody][1] *= 0.5;
    angmom[ibody][2] *= 0.5;  
    
    MathExtra::angmom_to_omega(angmom[ibody],ex_space[ibody],ey_space[ibody],
			       ez_space[ibody],inertia[ibody],omega[ibody]);
  }
  
  // set velocity/rotation of atoms in rigid bodies
  // virial is already setup from initial_integrate

  set_v();
}
