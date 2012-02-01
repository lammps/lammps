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
#include "stdlib.h"
#include "string.h"
#include "fix_rigid_nvt.h"
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

FixRigidNVT::FixRigidNVT(LAMMPS *lmp, int narg, char **arg) :
  FixRigid(lmp, narg, arg)
{ 
  // other settings are made by FixRigid parent

  scalar_flag = 1;
  restart_global = 1;
  extscalar = 1;
  
  // error checking
  // convert input period to frequency

  if (tempflag == 0)
    error->all(FLERR,"Did not set temp for fix rigid/nvt");
  if (t_start < 0.0 || t_stop <= 0.0)
    error->all(FLERR,"Target temperature for fix rigid/nvt cannot be 0.0");
  if (t_period <= 0.0) error->all(FLERR,"Fix rigid/nvt period must be > 0.0");
  t_freq = 1.0 / t_period;

  if (t_chain < 1) error->all(FLERR,"Illegal fix_modify command");
  if (t_iter < 1) error->all(FLERR,"Illegal fix_modify command");
  if (t_order != 3 && t_order != 5) 
    error->all(FLERR,"Fix_modify order must be 3 or 5"); 
  
  allocate_chain();
  allocate_order();
  memory->create(conjqm,nbody,4,"nve_rigid:conjqm");
  
  // one-time initialize of thermostat variables
  
  eta_t[0] = eta_r[0] = 0.0;
  eta_dot_t[0] = eta_dot_r[0] = 0.0;
  f_eta_t[0] = f_eta_r[0] = 0.0;
  
  for (int i = 1; i < t_chain; i++) {
    eta_t[i] = eta_r[i] = 0.0;
    eta_dot_t[i] = eta_dot_r[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

FixRigidNVT::~FixRigidNVT()
{
  deallocate_chain();
  deallocate_order();
  memory->destroy(conjqm);
}

/* ---------------------------------------------------------------------- */

int FixRigidNVT::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= PRE_NEIGHBOR;
  mask |= THERMO_ENERGY;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRigidNVT::init()
{
  FixRigid::init();

  // initialize thermostats
  // set timesteps, constants 
  // store Yoshida-Suzuki integrator parameters
  
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dtq = 0.5 * update->dt;
  
  boltz = force->boltz;
  nf_t = nf_r = domain->dimension * nbody;
	
  for (int ibody = 0; ibody < nbody; ibody++)
    for (int k = 0; k < domain->dimension; k++)
      if (fabs(inertia[ibody][k]) < 1e-6) nf_r--;
  
  // see Table 1 in Kamberaj et al

  if (t_order == 3) {
    w[0] = 1.0 / (2.0 - pow(2.0, 1.0/3.0));
    w[1] = 1.0 - 2.0*w[0];
    w[2] = w[0];
  } else if (t_order == 5) {
    w[0] = 1.0 / (4.0 - pow(4.0, 1.0/3.0));
    w[1] = w[0];
    w[2] = 1.0 - 4.0 * w[0];
    w[3] = w[0];
    w[4] = w[0];
  }
  
  // initialize thermostat settings
  
  t_target = t_start;
  double kt = boltz * t_target;
  double t_mass = kt / (t_freq*t_freq);
  q_t[0] = nf_t * t_mass;
  q_r[0] = nf_r * t_mass;
  for (int i = 1; i < t_chain; i++)
    q_t[i] = q_r[i] = t_mass;
  
  // initialize thermostat chain positions, velocites, forces
  
  for (int i = 1; i < t_chain; i++) {
    f_eta_t[i] = q_t[i-1] * eta_dot_t[i-1] * eta_dot_t[i-1] - kt;
    f_eta_r[i] = q_r[i-1] * eta_dot_r[i-1] * eta_dot_r[i-1] - kt;
  }
}

/* ---------------------------------------------------------------------- */

void FixRigidNVT::setup(int vflag)
{
  FixRigid::setup(vflag);

  t_target = t_start;
  
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

void FixRigidNVT::initial_integrate(int vflag)
{
  double tmp,akin_t,akin_r,scale_t,scale_r;
  double dtfm,mbody[3],tbody[3],fquat[4];
  double dtf2 = dtf * 2.0;
  
  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  t_target = t_start + delta * (t_stop - t_start);
  
  // intialize velocity scale for translation and rotation
  
  akin_t = akin_r = 0.0;
  tmp = -1.0 * dtq * eta_dot_t[0];
  scale_t = exp(tmp);
  tmp = -1.0 * dtq * eta_dot_r[0];
  scale_r = exp(tmp);
  
  for (int ibody = 0; ibody < nbody; ibody++) {
    
    // step 1.1 - update vcm by 1/2 step
    
    dtfm = dtf / masstotal[ibody];
    vcm[ibody][0] += dtfm * fcm[ibody][0] * fflag[ibody][0];
    vcm[ibody][1] += dtfm * fcm[ibody][1] * fflag[ibody][1];
    vcm[ibody][2] += dtfm * fcm[ibody][2] * fflag[ibody][2];
    vcm[ibody][0] *= scale_t;
    vcm[ibody][1] *= scale_t;
    vcm[ibody][2] *= scale_t;
    
    tmp = vcm[ibody][0]*vcm[ibody][0] + vcm[ibody][1]*vcm[ibody][1] +
      vcm[ibody][2]*vcm[ibody][2];
    akin_t += masstotal[ibody]*tmp;
    
    // step 1.2 - update xcm by full step
    
    xcm[ibody][0] += dtv * vcm[ibody][0];
    xcm[ibody][1] += dtv * vcm[ibody][1];
    xcm[ibody][2] += dtv * vcm[ibody][2];
    
    torque[ibody][0] *= tflag[ibody][0];
    torque[ibody][1] *= tflag[ibody][1];
    torque[ibody][2] *= tflag[ibody][2];
    
    // step 1.3 - apply torque (body coords) to quaternion momentum

    MathExtra::transpose_matvec(ex_space[ibody],ey_space[ibody],ez_space[ibody],
				torque[ibody],tbody);
    MathExtra::quatvec(quat[ibody],tbody,fquat);
    
    conjqm[ibody][0] += dtf2 * fquat[0];
    conjqm[ibody][1] += dtf2 * fquat[1];
    conjqm[ibody][2] += dtf2 * fquat[2];
    conjqm[ibody][3] += dtf2 * fquat[3];
    conjqm[ibody][0] *= scale_r;
    conjqm[ibody][1] *= scale_r;
    conjqm[ibody][2] *= scale_r;
    conjqm[ibody][3] *= scale_r;
  
    // step 1.4 to 1.8 - use no_squish rotate to update p (i.e. conjqm) and q
  
    no_squish_rotate(3,conjqm[ibody],quat[ibody],inertia[ibody],dtq);
    no_squish_rotate(2,conjqm[ibody],quat[ibody],inertia[ibody],dtq);
    no_squish_rotate(1,conjqm[ibody],quat[ibody],inertia[ibody],dtv);
    no_squish_rotate(2,conjqm[ibody],quat[ibody],inertia[ibody],dtq);
    no_squish_rotate(3,conjqm[ibody],quat[ibody],inertia[ibody],dtq);
  
    // update the exyz_space from new quaternion
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
    
    akin_r += angmom[ibody][0]*omega[ibody][0] + 
      angmom[ibody][1]*omega[ibody][1] + angmom[ibody][2]*omega[ibody][2];
  }
  
  // update thermostat chains
  
  update_nhcp(akin_t,akin_r);
  
  // virial setup before call to set_xv
  
  if (vflag) v_setup(vflag);
  else evflag = 0;
  
  // set coords/orient and velocity/rotation of atoms in rigid bodies
  // from quarternion and omega
  
  set_xv();
}

/* ---------------------------------------------------------------------- */

void FixRigidNVT::final_integrate()
{
  int i,ibody;
  double tmp,scale_t,scale_r;
  double dtfm,xy,xz,yz;
  
  // compute velocity scales for translation and rotation
  
  tmp = -1.0 * dtq * eta_dot_t[0];
  scale_t = exp(tmp);
  tmp = -1.0 * dtq * eta_dot_r[0];
  scale_r = exp(tmp);

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
  
  double mbody[3],tbody[3],fquat[4];
  double dtf2 = dtf * 2.0;
  
  for (ibody = 0; ibody < nbody; ibody++) {
    fcm[ibody][0] = all[ibody][0];
    fcm[ibody][1] = all[ibody][1];
    fcm[ibody][2] = all[ibody][2];
    torque[ibody][0] = all[ibody][3];
    torque[ibody][1] = all[ibody][4];
    torque[ibody][2] = all[ibody][5];

    // 2.5-2.6 update vcm by 1/2 step
  
    dtfm = dtf / masstotal[ibody];
    vcm[ibody][0] *= scale_t;
    vcm[ibody][1] *= scale_t;
    vcm[ibody][2] *= scale_t;
    vcm[ibody][0] += dtfm * fcm[ibody][0] * fflag[ibody][0];
    vcm[ibody][1] += dtfm * fcm[ibody][1] * fflag[ibody][1];
    vcm[ibody][2] += dtfm * fcm[ibody][2] * fflag[ibody][2];
    
    // 2.1-2.4 update conjqm, angular momentum and angular velocity
    // apply body torque flags  
    
    torque[ibody][0] *= tflag[ibody][0];
    torque[ibody][1] *= tflag[ibody][1];
    torque[ibody][2] *= tflag[ibody][2];
    
    // convert torque to the body frame 
    
    MathExtra::transpose_matvec(ex_space[ibody],ey_space[ibody],ez_space[ibody],
				torque[ibody],tbody);
    
    // compute "force" for quaternion
    
    MathExtra::quatvec(quat[ibody],tbody,fquat);
    
    // update the conjugate quaternion momentum (conjqm)
    
    conjqm[ibody][0] = scale_r * conjqm[ibody][0] + dtf2 * fquat[0];
    conjqm[ibody][1] = scale_r * conjqm[ibody][1] + dtf2 * fquat[1];
    conjqm[ibody][2] = scale_r * conjqm[ibody][2] + dtf2 * fquat[2];
    conjqm[ibody][3] = scale_r * conjqm[ibody][3] + dtf2 * fquat[3];
    
    // compute angular momentum in the body frame
    // then convert to the space-fixed frame
    
    MathExtra::invquatvec(quat[ibody],conjqm[ibody],mbody);
    MathExtra::matvec(ex_space[ibody],ey_space[ibody],ez_space[ibody],
		      mbody,angmom[ibody]);
    
    angmom[ibody][0] *= 0.5;
    angmom[ibody][1] *= 0.5;
    angmom[ibody][2] *= 0.5;  
    
    // compute new angular velocity
    
    MathExtra::angmom_to_omega(angmom[ibody],ex_space[ibody],ey_space[ibody],
			       ez_space[ibody],inertia[ibody],omega[ibody]);
  }
  
  // set velocity/rotation of atoms in rigid bodies
  // virial is already setup from initial_integrate

  set_v();
}

/* ---------------------------------------------------------------------- */

void FixRigidNVT::update_nhcp(double akin_t, double akin_r)
{
  int i,j,k;
  double kt,gfkt_t,gfkt_r,tmp,ms,s,s2;
  
  kt = boltz * t_target;
  gfkt_t = nf_t * kt;
  gfkt_r = nf_r * kt;
  
  akin_t *= force->mvv2e;
  akin_r *= force->mvv2e;
  
  // update thermostat masses
  
  double t_mass = boltz * t_target / (t_freq * t_freq);
  q_t[0] = nf_t * t_mass;
  q_r[0] = nf_r * t_mass;
  for (i = 1; i < t_chain; i++)
    q_t[i] = q_r[i] = t_mass;
  
  // update order/timestep dependent coefficients
  
  for (i = 0; i < t_order; i++) {
    wdti1[i] = w[i] * dtv / t_iter;
    wdti2[i] = wdti1[i] / 2.0;
    wdti4[i] = wdti1[i] / 4.0;
  }
  
  // update force of thermostats coupled to particles
  
  f_eta_t[0] = (akin_t - gfkt_t) / q_t[0];
  f_eta_r[0] = (akin_r - gfkt_r) / q_r[0];
  
  // multiple timestep iteration
  
  for (i = 0; i < t_iter; i++) {
    for (j = 0; j < t_order; j++) {
  
      // update thermostat velocities half step
  
      eta_dot_t[t_chain-1] += wdti2[j] * f_eta_t[t_chain-1];
      eta_dot_r[t_chain-1] += wdti2[j] * f_eta_r[t_chain-1];
      
      for (k = 1; k < t_chain; k++) {
	tmp = wdti4[j] * eta_dot_t[t_chain-k];
	ms = maclaurin_series(tmp);
	s = exp(-1.0 * tmp);
	s2 = s * s;
	eta_dot_t[t_chain-k-1] = eta_dot_t[t_chain-k-1] * s2 + 
	  wdti2[j] * f_eta_t[t_chain-k-1] * s * ms;
	
	tmp = wdti4[j] * eta_dot_r[t_chain-k];
	ms = maclaurin_series(tmp);
	s = exp(-1.0 * tmp);
	s2 = s * s;
	eta_dot_r[t_chain-k-1] = eta_dot_r[t_chain-k-1] * s2 + 
	  wdti2[j] * f_eta_r[t_chain-k-1] * s * ms;
      }
      
      // update thermostat positions a full step
      
      for (k = 0; k < t_chain; k++) {
	eta_t[k] += wdti1[j] * eta_dot_t[k];
	eta_r[k] += wdti1[j] * eta_dot_r[k];
      }
      
      // update thermostat forces 
      
      for (k = 1; k < t_chain; k++) {
	f_eta_t[k] = q_t[k-1] * eta_dot_t[k-1] * eta_dot_t[k-1] - kt;
	f_eta_t[k] /= q_t[k];
	f_eta_r[k] = q_r[k-1] * eta_dot_r[k-1] * eta_dot_r[k-1] - kt;
	f_eta_r[k] /= q_r[k];
      }
      
      // update thermostat velocities a full step
      
      for (k = 0; k < t_chain-1; k++) {
	tmp = wdti4[j] * eta_dot_t[k+1];
	ms = maclaurin_series(tmp);
	s = exp(-1.0 * tmp);
	s2 = s * s;
	eta_dot_t[k] = eta_dot_t[k] * s2 + wdti2[j] * f_eta_t[k] * s * ms;
	tmp = q_t[k] * eta_dot_t[k] * eta_dot_t[k] - kt;
	f_eta_t[k+1] = tmp / q_t[k+1];
	
	tmp = wdti4[j] * eta_dot_r[k+1];
	ms = maclaurin_series(tmp);
	s = exp(-1.0 * tmp);
	s2 = s * s;
	eta_dot_r[k] = eta_dot_r[k] * s2 + wdti2[j] * f_eta_r[k] * s * ms;
	tmp = q_r[k] * eta_dot_r[k] * eta_dot_r[k] - kt;
	f_eta_r[k+1] = tmp / q_r[k+1];
      }
      
      eta_dot_t[t_chain-1] += wdti2[j] * f_eta_t[t_chain-1];
      eta_dot_r[t_chain-1] += wdti2[j] * f_eta_r[t_chain-1];
    }
  }
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write 
------------------------------------------------------------------------- */

void FixRigidNVT::write_restart(FILE *fp)
{
  int n = 0;
  double *list = new double[1+6*t_chain];
  list[n++] = t_chain;

  for (int i = 0; i < t_chain; i++) {
    list[n++] = eta_t[i];
    list[n++] = eta_r[i];
    list[n++] = eta_dot_t[i];
    list[n++] = eta_dot_r[i];
    list[n++] = f_eta_t[i];
    list[n++] = f_eta_r[i];
  }
  
  if (comm->me == 0) {
    int size = (1 + 6*t_chain)*sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),1+6*t_chain,fp);
  }
  
  delete list;
}

/* ---------------------------------------------------------------------- 
   compute kinetic energy in the extended Hamiltonian
   conserved quantity = sum of returned energy and potential energy
-----------------------------------------------------------------------*/

double FixRigidNVT::compute_scalar()
{
  int i,k,ibody;
  double kt = boltz * t_target;
  double energy,ke_t,ke_q,tmp,Pkq[4];
  
  // compute the kinetic parts of H_NVE in Kameraj et al (JCP 2005, pp 224114)
  
  // translational kinetic energy

  ke_t = 0.0;
  for (ibody = 0; ibody < nbody; ibody++)
    ke_t += 0.5 * masstotal[ibody] * (vcm[ibody][0]*vcm[ibody][0] +
				      vcm[ibody][1]*vcm[ibody][1] +
				      vcm[ibody][2]*vcm[ibody][2]);
  
  // rotational kinetic energy

  ke_q = 0.0;
  for (ibody = 0; ibody < nbody; ibody++) {
    for (k = 1; k < 4; k++) {
      if (k == 1) {
        Pkq[0] = -quat[ibody][1];
        Pkq[1] =  quat[ibody][0];
        Pkq[2] =  quat[ibody][3];
        Pkq[3] = -quat[ibody][2];
      } else if (k == 2) {
        Pkq[0] = -quat[ibody][2];
        Pkq[1] = -quat[ibody][3];
        Pkq[2] =  quat[ibody][0];
        Pkq[3] =  quat[ibody][1];
      } else if (k == 3) {
        Pkq[0] = -quat[ibody][3];
        Pkq[1] =  quat[ibody][2];
        Pkq[2] = -quat[ibody][1];
        Pkq[3] =  quat[ibody][0];      
      }
   
      tmp = conjqm[ibody][0]*Pkq[0] + conjqm[ibody][1]*Pkq[1] +
	conjqm[ibody][2]*Pkq[2] + conjqm[ibody][3]*Pkq[3];
      tmp *= tmp;
    
      if (fabs(inertia[ibody][k-1]) < 1e-6) tmp = 0.0;
      else tmp /= (8.0 * inertia[ibody][k-1]); 
      ke_q += tmp;
    }
  }
  
  energy = ke_t + ke_q;
  
  // thermostat chain energy: from equation 12 in Kameraj et al (JCP 2005)

  energy += kt * (nf_t * eta_t[0] + nf_r * eta_r[0]);
  
  for (i = 1; i < t_chain; i++) 
    energy += kt * (eta_t[i] + eta_r[i]);
  
  for (i = 0;  i < t_chain; i++) {
    energy += 0.5 * q_t[i] * (eta_dot_t[i] * eta_dot_t[i]);
    energy += 0.5 * q_r[i] * (eta_dot_r[i] * eta_dot_r[i]);
  }
    
  return energy;
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix 
------------------------------------------------------------------------- */

void FixRigidNVT::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;
  
  int t_chain_prev = static_cast<int> (list[n++]);
  if (t_chain_prev != t_chain)
    error->all(FLERR,"Cannot restart fix rigid/nvt with different # of chains");

  for (int i = 0; i < t_chain; i++) {
    eta_t[i] = list[n++];
    eta_r[i] = list[n++];
    eta_dot_t[i] = list[n++];
    eta_dot_r[i] = list[n++];
    f_eta_t[i] = list[n++];
    f_eta_r[i] = list[n++];
  }
}

/* ---------------------------------------------------------------------- */

void FixRigidNVT::reset_target(double t_new)
{
  t_start = t_stop = t_new;
}

/* ---------------------------------------------------------------------- */

void FixRigidNVT::allocate_chain()
{
  q_t = new double[t_chain];
  q_r = new double[t_chain];
  eta_t = new double[t_chain];
  eta_r = new double[t_chain];
  eta_dot_t = new double[t_chain];
  eta_dot_r = new double[t_chain];
  f_eta_t = new double[t_chain];
  f_eta_r = new double[t_chain];
}

/* ---------------------------------------------------------------------- */

void FixRigidNVT::allocate_order()
{
  w = new double[t_order];
  wdti1 = new double[t_order];
  wdti2 = new double[t_order];
  wdti4 = new double[t_order];
}

/* ---------------------------------------------------------------------- */

void FixRigidNVT::deallocate_chain()
{
  delete [] q_t;
  delete [] q_r;
  delete [] eta_t;
  delete [] eta_r;
  delete [] eta_dot_t;
  delete [] eta_dot_r;
  delete [] f_eta_t;
  delete [] f_eta_r;
}

/* ---------------------------------------------------------------------- */

void FixRigidNVT::deallocate_order()
{
  delete [] w;
  delete [] wdti1;
  delete [] wdti2;
  delete [] wdti4;
}
