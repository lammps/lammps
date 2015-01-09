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

#include "math.h"
#include "stdio.h"
#include "string.h"
#include "fix_rigid_nh.h"
#include "math_extra.h"
#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "update.h"
#include "modify.h"
#include "fix_deform.h"
#include "group.h"
#include "comm.h"
#include "force.h"
#include "kspace.h"
#include "output.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,XYZ,XY,YZ,XZ};     // same as in FixRigid
enum{ISO,ANISO,TRICLINIC};   // same as in FixRigid

#define EPSILON 1.0e-7

/* ---------------------------------------------------------------------- */

FixRigidNH::FixRigidNH(LAMMPS *lmp, int narg, char **arg) :
  FixRigid(lmp, narg, arg)
{
  // error checks: could be moved up to FixRigid
  
  if ((p_flag[0] == 1 && p_period[0] <= 0.0) || 
      (p_flag[1] == 1 && p_period[1] <= 0.0) || 
      (p_flag[2] == 1 && p_period[2] <= 0.0)) 
    error->all(FLERR,"Fix rigid npt/nph period must be > 0.0");
  
  if (dimension == 2 && p_flag[2])
    error->all(FLERR,"Invalid fix rigid npt/nph command for a 2d simulation");
  if (dimension == 2 && (pcouple == YZ || pcouple == XZ))
    error->all(FLERR,"Invalid fix rigid npt/nph command for a 2d simulation");

  if (pcouple == XYZ && (p_flag[0] == 0 || p_flag[1] == 0))
    error->all(FLERR,"Invalid fix rigid npt/nph command pressure settings");
  if (pcouple == XYZ && dimension == 3 && p_flag[2] == 0)
    error->all(FLERR,"Invalid fix rigid npt/nph command pressure settings");
  if (pcouple == XY && (p_flag[0] == 0 || p_flag[1] == 0))
    error->all(FLERR,"Invalid fix rigid npt/nph command pressure settings");
  if (pcouple == YZ && (p_flag[1] == 0 || p_flag[2] == 0))
    error->all(FLERR,"Invalid fix rigid npt/nph command pressure settings");
  if (pcouple == XZ && (p_flag[0] == 0 || p_flag[2] == 0))
    error->all(FLERR,"Invalid fix rigid npt/nph command pressure settings");

  // require periodicity in tensile dimension

  if (p_flag[0] && domain->xperiodic == 0)
    error->all(FLERR,
               "Cannot use fix rigid npt/nph on a non-periodic dimension");
  if (p_flag[1] && domain->yperiodic == 0)
    error->all(FLERR,
               "Cannot use fix rigid npt/nph on a non-periodic dimension");
  if (p_flag[2] && domain->zperiodic == 0)
    error->all(FLERR,
               "Cannot use fix rigid npt/nph on a non-periodic dimension");
  
    if (pcouple == XYZ && dimension == 3 &&
      (p_start[0] != p_start[1] || p_start[0] != p_start[2] ||
       p_stop[0] != p_stop[1] || p_stop[0] != p_stop[2] ||
       p_period[0] != p_period[1] || p_period[0] != p_period[2]))
    error->all(FLERR,"Invalid fix rigid npt/nph pressure settings");
  if (pcouple == XYZ && dimension == 2 &&
      (p_start[0] != p_start[1] || p_stop[0] != p_stop[1] ||
       p_period[0] != p_period[1]))
    error->all(FLERR,"Invalid fix rigid npt/nph pressure settings");
  if (pcouple == XY &&
      (p_start[0] != p_start[1] || p_stop[0] != p_stop[1] ||
       p_period[0] != p_period[1]))
    error->all(FLERR,"Invalid fix rigid npt/nph pressure settings");
  if (pcouple == YZ &&
      (p_start[1] != p_start[2] || p_stop[1] != p_stop[2] ||
       p_period[1] != p_period[2]))
    error->all(FLERR,"Invalid fix rigid npt/nph pressure settings");
  if (pcouple == XZ &&
      (p_start[0] != p_start[2] || p_stop[0] != p_stop[2] ||
       p_period[0] != p_period[2]))
    error->all(FLERR,"Invalid fix rigid npt/nph pressure settings");

  if ((tstat_flag && t_period <= 0.0) ||
      (p_flag[0] && p_period[0] <= 0.0) ||
      (p_flag[1] && p_period[1] <= 0.0) ||
      (p_flag[2] && p_period[2] <= 0.0))
    error->all(FLERR,"Fix rigid nvt/npt/nph damping parameters must be > 0.0");

  // memory allocation and initialization
  
  memory->create(conjqm,nbody,4,"rigid_nh:conjqm");
  if (tstat_flag || pstat_flag) {
    allocate_chain();
    allocate_order();
  }
  
  if (tstat_flag) {
    eta_t[0] = eta_r[0] = 0.0;
    eta_dot_t[0] = eta_dot_r[0] = 0.0;
    f_eta_t[0] = f_eta_r[0] = 0.0;
  
    for (int i = 1; i < t_chain; i++) {
      eta_t[i] = eta_r[i] = 0.0;
      eta_dot_t[i] = eta_dot_r[i] = 0.0;
    }
  }
  
  if (pstat_flag) {
    epsilon_dot[0] = epsilon_dot[1] = epsilon_dot[2] = 0.0;
    eta_b[0] = eta_dot_b[0] = f_eta_b[0] = 0.0;
    for (int i = 1; i < p_chain; i++) 
      eta_b[i] = eta_dot_b[i] = 0.0;
  }

  // rigid body pointers
  
  nrigidfix = 0;
  rfix = NULL;

  vol0 = 0.0;
  t0 = 1.0;
  
  tcomputeflag = 0;
  pcomputeflag = 0;
}

/* ---------------------------------------------------------------------- */

FixRigidNH::~FixRigidNH()
{
  memory->destroy(conjqm);
  if (tstat_flag || pstat_flag) {
    deallocate_chain();
    deallocate_order();
  }
  
  if (rfix) delete [] rfix;

  if (tcomputeflag) {
    modify->delete_compute(id_temp);
    delete [] id_temp;
  }
  
  // delete pressure if fix created it
  
  if (pstat_flag) {
    if (pcomputeflag) modify->delete_compute(id_press);
    delete [] id_press;
  }
}

/* ---------------------------------------------------------------------- */

int FixRigidNH::setmask()
{
  int mask = 0;
  mask = FixRigid::setmask();
  if (tstat_flag || pstat_flag) mask |= THERMO_ENERGY;
  
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRigidNH::init()
{
  FixRigid::init();

  // recheck that dilate group has not been deleted

  if (allremap == 0) {
    int idilate = group->find(id_dilate);
    if (idilate == -1)
      error->all(FLERR,"Fix rigid npt/nph dilate group ID does not exist");
    dilate_group_bit = group->bitmask[idilate];
  }

  // initialize thermostats
  // set timesteps, constants 
  // store Yoshida-Suzuki integrator parameters
  
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dtq = 0.5 * update->dt;
  
  boltz = force->boltz;
  nktv2p = force->nktv2p;
  mvv2e = force->mvv2e;
  
  if (force->kspace) kspace_flag = 1;
  else kspace_flag = 0;

  nf_t = dimension * nbody;  
  if (dimension == 3) {
    nf_r = dimension * nbody;
    for (int ibody = 0; ibody < nbody; ibody++)
      for (int k = 0; k < domain->dimension; k++) 
        if (fabs(inertia[ibody][k]) < EPSILON) nf_r--;
  } else if (dimension == 2) {
    nf_r = nbody;
    for (int ibody = 0; ibody < nbody; ibody++)
      if (fabs(inertia[ibody][2]) < EPSILON) nf_r--;
  }

  g_f = nf_t + nf_r;  
  onednft = 1.0 + (double)(dimension) / (double)g_f;
  onednfr = (double) (dimension) / (double)g_f;

  // see Table 1 in Kamberaj et al
  
  if (tstat_flag || pstat_flag) {
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
  }

  int icompute;
  if (tcomputeflag) {  
    icompute = modify->find_compute(id_temp);
    if (icompute < 0) 
      error->all(FLERR,"Temperature ID for fix rigid nvt/npt/nph does not exist");
    temperature = modify->compute[icompute];
  }

  if (pstat_flag) {
    if (domain->triclinic) 
      error->all(FLERR,"Fix rigid npt/nph does not yet allow triclinic box");
  
    // ensure no conflict with fix deform

    for (int i = 0; i < modify->nfix; i++)
      if (strcmp(modify->fix[i]->style,"deform") == 0) {
      	int *dimflag = ((FixDeform *) modify->fix[i])->dimflag;
      	if ((p_flag[0] && dimflag[0]) || (p_flag[1] && dimflag[1]) || 
      	    (p_flag[2] && dimflag[2]))
          error->all(FLERR,"Cannot use fix rigid npt/nph and fix deform on "
                     "same component of stress tensor");
      }

    // set frequency
  
    p_freq_max = 0.0;
    p_freq_max = MAX(p_freq[0],p_freq[1]);
    p_freq_max = MAX(p_freq_max,p_freq[2]);

    // tally the number of dimensions that are barostatted
    // set initial volume and reference cell, if not already done

    pdim = p_flag[0] + p_flag[1] + p_flag[2];
    if (vol0 == 0.0) {
      if (dimension == 2) vol0 = domain->xprd * domain->yprd;
      else vol0 = domain->xprd * domain->yprd * domain->zprd;
    }

    // set pressure compute ptr

    icompute = modify->find_compute(id_press);
    if (icompute < 0) 
      error->all(FLERR,"Pressure ID for fix rigid npt/nph does not exist");
    pressure = modify->compute[icompute];
    
    // detect if any rigid fixes exist so rigid bodies move on remap
    // rfix[] = indices to each fix rigid
    // this will include self

    if (rfix) delete [] rfix;
    nrigidfix = 0;
    rfix = NULL;

    for (int i = 0; i < modify->nfix; i++)
      if (modify->fix[i]->rigid_flag) nrigidfix++;
    if (nrigidfix) {
      rfix = new int[nrigidfix];
      nrigidfix = 0;
      for (int i = 0; i < modify->nfix; i++)
        if (modify->fix[i]->rigid_flag) rfix[nrigidfix++] = i;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixRigidNH::setup(int vflag)
{
  FixRigid::setup(vflag);
  
  double mbody[3];
  akin_t = akin_r = 0.0;
  for (int ibody = 0; ibody < nbody; ibody++) {
    MathExtra::transpose_matvec(ex_space[ibody],ey_space[ibody],ez_space[ibody],
                                angmom[ibody],mbody);
    MathExtra::quatvec(quat[ibody],mbody,conjqm[ibody]);
    conjqm[ibody][0] *= 2.0;
    conjqm[ibody][1] *= 2.0;
    conjqm[ibody][2] *= 2.0;
    conjqm[ibody][3] *= 2.0;
    
    if (tstat_flag || pstat_flag) {
      akin_t += masstotal[ibody]*(vcm[ibody][0]*vcm[ibody][0] + 
        vcm[ibody][1]*vcm[ibody][1] + vcm[ibody][2]*vcm[ibody][2]);
      akin_r += angmom[ibody][0]*omega[ibody][0] + 
        angmom[ibody][1]*omega[ibody][1] + angmom[ibody][2]*omega[ibody][2];
    }
  }

  // compute target temperature
  
  if (tstat_flag) compute_temp_target();
  else if (pstat_flag) {
    t0 = temperature->compute_scalar();
    if (t0 == 0.0) {
      if (strcmp(update->unit_style,"lj") == 0) t0 = 1.0;
      else t0 = 300.0;
    }
    t_target = t0;
  }

  // compute target pressure
  // compute current pressure
  // trigger virial computation on next timestep
    
  if (pstat_flag) { 
    if (pstyle == ISO) {
      temperature->compute_scalar();
      pressure->compute_scalar();
    } else {
      temperature->compute_vector();
      pressure->compute_vector();
    }
    couple();
    pressure->addstep(update->ntimestep+1);
  }
  
  // initialize thermostat/barostat settings
  
  double kt, t_mass, tb_mass;
  kt = boltz * t_target;

  if (tstat_flag) {
    t_mass = kt / (t_freq*t_freq);
    q_t[0] = nf_t * t_mass;
    q_r[0] = nf_r * t_mass;
    for (int i = 1; i < t_chain; i++) 
      q_t[i] = q_r[i] = t_mass;

    for (int i = 1; i < t_chain; i++) {
      f_eta_t[i] = (q_t[i-1] * eta_dot_t[i-1] * eta_dot_t[i-1] - kt)/q_t[i];
      f_eta_r[i] = (q_r[i-1] * eta_dot_r[i-1] * eta_dot_r[i-1] - kt)/q_r[i];
    }
  }
  
  // initial forces on barostat thermostat variables
  
  if (pstat_flag) {
    for (int i = 0; i < 3; i++)
      if (p_flag[i]) {
        epsilon_mass[i] = (g_f + dimension) * kt / (p_freq[i]*p_freq[i]);
        epsilon[i] = log(vol0)/dimension;
      } 
    
    tb_mass = kt / (p_freq_max * p_freq_max);
    q_b[0] = dimension * dimension * tb_mass;
    for (int i = 1; i < p_chain; i++) {
      q_b[i] = tb_mass;
      f_eta_b[i] = (q_b[i] * eta_dot_b[i-1] * eta_dot_b[i-1] - kt)/q_b[i];
    }
  }
  
  // update order/timestep dependent coefficients
  
  if (tstat_flag || pstat_flag) {
    for (int i = 0; i < t_order; i++) {
      wdti1[i] = w[i] * dtv / t_iter;
      wdti2[i] = wdti1[i] / 2.0;
      wdti4[i] = wdti1[i] / 4.0;
    }
  }

  if (pstat_flag) { 
    compute_press_target();
    nh_epsilon_dot();
  }
}

/* ----------------------------------------------------------------------
   perform preforce velocity Verlet integration
   see Kamberaj paper for step references
------------------------------------------------------------------------- */

void FixRigidNH::initial_integrate(int vflag)
{
  double tmp,scale_r,scale_t[3],scale_v[3];
  double dtfm,mbody[3],tbody[3],fquat[4];
  double dtf2 = dtf * 2.0;
  
  // compute scale variables

  scale_t[0] = scale_t[1] = scale_t[2] = 1.0;
  scale_v[0] = scale_v[1] = scale_v[2] = 1.0;
  scale_r = 1.0;

  if (tstat_flag) {
    akin_t = akin_r = 0.0;
    tmp = exp(-dtq * eta_dot_t[0]);
    scale_t[0] = scale_t[1] = scale_t[2] = tmp;
    tmp = exp(-dtq * eta_dot_r[0]);
    scale_r = tmp;
  } 

  if (pstat_flag) {
    akin_t = akin_r = 0.0;
    scale_t[0] *= exp(-dtq * (epsilon_dot[0] + mtk_term2));
    scale_t[1] *= exp(-dtq * (epsilon_dot[1] + mtk_term2));
    scale_t[2] *= exp(-dtq * (epsilon_dot[2] + mtk_term2));
    scale_r *= exp(-dtq * (pdim * mtk_term2));

    tmp = dtq * epsilon_dot[0];
    scale_v[0] = dtv * exp(tmp) * maclaurin_series(tmp);
    tmp = dtq * epsilon_dot[1];
    scale_v[1] = dtv * exp(tmp) * maclaurin_series(tmp);
    tmp = dtq * epsilon_dot[2];
    scale_v[2] = dtv * exp(tmp) * maclaurin_series(tmp);
  }
    
  // update xcm, vcm, quat, conjqm and angmom

  for (int ibody = 0; ibody < nbody; ibody++) {
    
    // step 1.1 - update vcm by 1/2 step
    
    dtfm = dtf / masstotal[ibody];
    vcm[ibody][0] += dtfm * fcm[ibody][0] * fflag[ibody][0];
    vcm[ibody][1] += dtfm * fcm[ibody][1] * fflag[ibody][1];
    vcm[ibody][2] += dtfm * fcm[ibody][2] * fflag[ibody][2];
    
    if (tstat_flag || pstat_flag) {
      vcm[ibody][0] *= scale_t[0];
      vcm[ibody][1] *= scale_t[1];
      vcm[ibody][2] *= scale_t[2];
      
      tmp = vcm[ibody][0]*vcm[ibody][0] + vcm[ibody][1]*vcm[ibody][1] +
        vcm[ibody][2]*vcm[ibody][2];
      akin_t += masstotal[ibody]*tmp;
    }
    
    // step 1.2 - update xcm by full step

    if (!pstat_flag) {
      xcm[ibody][0] += dtv * vcm[ibody][0];
      xcm[ibody][1] += dtv * vcm[ibody][1];
      xcm[ibody][2] += dtv * vcm[ibody][2];
    } else {
      xcm[ibody][0] += scale_v[0] * vcm[ibody][0];
      xcm[ibody][1] += scale_v[1] * vcm[ibody][1];
      xcm[ibody][2] += scale_v[2] * vcm[ibody][2];
    }
    
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
    
    if (tstat_flag || pstat_flag) {
      conjqm[ibody][0] *= scale_r;
      conjqm[ibody][1] *= scale_r;
      conjqm[ibody][2] *= scale_r;
      conjqm[ibody][3] *= scale_r;
    }
    
    // step 1.4 to 1.13 - use no_squish rotate to update p and q
  
    MathExtra::no_squish_rotate(3,conjqm[ibody],quat[ibody],inertia[ibody],dtq);
    MathExtra::no_squish_rotate(2,conjqm[ibody],quat[ibody],inertia[ibody],dtq);
    MathExtra::no_squish_rotate(1,conjqm[ibody],quat[ibody],inertia[ibody],dtv);
    MathExtra::no_squish_rotate(2,conjqm[ibody],quat[ibody],inertia[ibody],dtq);
    MathExtra::no_squish_rotate(3,conjqm[ibody],quat[ibody],inertia[ibody],dtq);
  
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
    
    if (tstat_flag || pstat_flag) {
      akin_r += angmom[ibody][0]*omega[ibody][0] + 
        angmom[ibody][1]*omega[ibody][1] + angmom[ibody][2]*omega[ibody][2];
    }
  }

  // compute target temperature
  // update thermostat chains using akin_t and akin_r
  // refer to update_nhcp() in Kamberaj et al.

  if (tstat_flag) {
    compute_temp_target();
    nhc_temp_integrate();
  }
    
  // update thermostat chains coupled with barostat
  // refer to update_nhcb() in Kamberaj et al.

  if (pstat_flag) {
    nhc_press_integrate();
  }

  // virial setup before call to set_xv

  if (vflag) v_setup(vflag);
  else evflag = 0;

  // remap simulation box by 1/2 step

  if (pstat_flag) remap();

  // set coords/orient and velocity/rotation of atoms in rigid bodies
  // from quarternion and omega

  set_xv();
  
  // remap simulation box by full step
  // redo KSpace coeffs since volume has changed

  if (pstat_flag) {
    remap();
    if (kspace_flag) force->kspace->setup();
  }  
}

/* ---------------------------------------------------------------------- */

void FixRigidNH::final_integrate()
{
  int i,ibody;
  double tmp,scale_t[3],scale_r;
  double dtfm,xy,xz,yz;
  double mbody[3],tbody[3],fquat[4];
  double dtf2 = dtf * 2.0;

  // compute scale variables
  
  scale_t[0] = scale_t[1] = scale_t[2] = 1.0;
  scale_r = 1.0;

  if (tstat_flag) {
    tmp = exp(-1.0 * dtq * eta_dot_t[0]);
    scale_t[0] = scale_t[1] = scale_t[2] = tmp;
    scale_r = exp(-1.0 * dtq * eta_dot_r[0]);
  } 
  
  if (pstat_flag) {
    scale_t[0] *= exp(-dtq * (epsilon_dot[0] + mtk_term2));
    scale_t[1] *= exp(-dtq * (epsilon_dot[1] + mtk_term2));
    scale_t[2] *= exp(-dtq * (epsilon_dot[2] + mtk_term2));
    scale_r *= exp(-dtq * (pdim * mtk_term2));
    
    // reset akin_t and akin_r, to be accumulated for use in nh_epsilon_dot()

    akin_t = akin_r = 0.0;
  }
  
  // sum over atoms to get force and torque on rigid body
  
  imageint *image = atom->image;
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
      
    xbox = (image[i] & IMGMASK) - IMGMAX;
    ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (image[i] >> IMG2BITS) - IMGMAX;

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

  for (ibody = 0; ibody < nbody; ibody++) {
    fcm[ibody][0] = all[ibody][0] + langextra[ibody][0];
    fcm[ibody][1] = all[ibody][1] + langextra[ibody][1];
    fcm[ibody][2] = all[ibody][2] + langextra[ibody][2];
    torque[ibody][0] = all[ibody][3] + langextra[ibody][3];
    torque[ibody][1] = all[ibody][4] + langextra[ibody][4];
    torque[ibody][2] = all[ibody][5] + langextra[ibody][5];

    // update vcm by 1/2 step
  
    dtfm = dtf / masstotal[ibody];
    if (tstat_flag || pstat_flag) {
      vcm[ibody][0] *= scale_t[0];
      vcm[ibody][1] *= scale_t[1];
      vcm[ibody][2] *= scale_t[2];
    }
    
    vcm[ibody][0] += dtfm * fcm[ibody][0] * fflag[ibody][0];
    vcm[ibody][1] += dtfm * fcm[ibody][1] * fflag[ibody][1];
    vcm[ibody][2] += dtfm * fcm[ibody][2] * fflag[ibody][2];
    
    if (pstat_flag) {
      tmp = vcm[ibody][0]*vcm[ibody][0] + vcm[ibody][1]*vcm[ibody][1] +
        vcm[ibody][2]*vcm[ibody][2];
      akin_t += masstotal[ibody]*tmp;
    }
    
    // update conjqm, then transform to angmom, set velocity again
    // virial is already setup from initial_integrate
    
    torque[ibody][0] *= tflag[ibody][0];
    torque[ibody][1] *= tflag[ibody][1];
    torque[ibody][2] *= tflag[ibody][2];
    
    MathExtra::transpose_matvec(ex_space[ibody],ey_space[ibody],
                                ez_space[ibody],torque[ibody],tbody);
    MathExtra::quatvec(quat[ibody],tbody,fquat);
    
    if (tstat_flag || pstat_flag) {
      conjqm[ibody][0] = scale_r * conjqm[ibody][0] + dtf2 * fquat[0];
      conjqm[ibody][1] = scale_r * conjqm[ibody][1] + dtf2 * fquat[1];
      conjqm[ibody][2] = scale_r * conjqm[ibody][2] + dtf2 * fquat[2];
      conjqm[ibody][3] = scale_r * conjqm[ibody][3] + dtf2 * fquat[3];
    } else {
      conjqm[ibody][0] += dtf2 * fquat[0];
      conjqm[ibody][1] += dtf2 * fquat[1];
      conjqm[ibody][2] += dtf2 * fquat[2];
      conjqm[ibody][3] += dtf2 * fquat[3];
    }

    MathExtra::invquatvec(quat[ibody],conjqm[ibody],mbody);
    MathExtra::matvec(ex_space[ibody],ey_space[ibody],ez_space[ibody],
                      mbody,angmom[ibody]);
    
    angmom[ibody][0] *= 0.5;
    angmom[ibody][1] *= 0.5;
    angmom[ibody][2] *= 0.5;  
    
    MathExtra::angmom_to_omega(angmom[ibody],ex_space[ibody],ey_space[ibody],
                               ez_space[ibody],inertia[ibody],omega[ibody]);
 
    if (pstat_flag) {
      akin_r += angmom[ibody][0]*omega[ibody][0] + 
        angmom[ibody][1]*omega[ibody][1] + 
        angmom[ibody][2]*omega[ibody][2];
    }
  }

  // set velocity/rotation of atoms in rigid bodies
  // virial is already setup from initial_integrate

  set_v();

  // compute current temperature
  if (tcomputeflag) t_current = temperature->compute_scalar();

  // compute current and target pressures
  // update epsilon dot using akin_t and akin_r
  
  if (pstat_flag) {
    if (pstyle == ISO) {
      temperature->compute_scalar();
      pressure->compute_scalar();
    } else {
      temperature->compute_vector();
      pressure->compute_vector();
    }
    couple();
    pressure->addstep(update->ntimestep+1);

    compute_press_target();

    nh_epsilon_dot();
  }
}

/* ---------------------------------------------------------------------- */

void FixRigidNH::nhc_temp_integrate()
{
  int i,j,k;
  double kt,gfkt_t,gfkt_r,tmp,ms,s,s2;
  
  kt = boltz * t_target;
  gfkt_t = nf_t * kt;
  gfkt_r = nf_r * kt;

  // update thermostat masses
  
  double t_mass = boltz * t_target / (t_freq * t_freq);
  q_t[0] = nf_t * t_mass;
  q_r[0] = nf_r * t_mass;
  for (i = 1; i < t_chain; i++)
    q_t[i] = q_r[i] = t_mass;
  
  // update force of thermostats coupled to particles
  
  f_eta_t[0] = (akin_t * mvv2e - gfkt_t) / q_t[0];
  f_eta_r[0] = (akin_r * mvv2e - gfkt_r) / q_r[0];
  
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

/* ---------------------------------------------------------------------- */

void FixRigidNH::nhc_press_integrate()
{
  int i,j,k;
  double tmp,s,s2,ms,kecurrent;
  double kt = boltz * t_target;
  double lkt_press = kt;
  
  // update thermostat masses
  
  double tb_mass = kt / (p_freq_max * p_freq_max);
  q_b[0] = dimension * dimension * tb_mass;
  for (int i = 1; i < p_chain; i++) {
    q_b[i] = tb_mass;
    f_eta_b[i] = q_b[i-1] * eta_dot_b[i-1] * eta_dot_b[i-1] - kt;
    f_eta_b[i] /= q_b[i];
  }
     
  // update forces acting on thermostat
  
  kecurrent = 0.0;
  for (i = 0; i < 3; i++) 
    if (p_flag[i]) {
      epsilon_mass[i] = (g_f + dimension) * kt / (p_freq[i] * p_freq[i]);
      kecurrent += epsilon_mass[i] * epsilon_dot[i] * epsilon_dot[i];
    }
  kecurrent /= pdim;

  f_eta_b[0] = (kecurrent - lkt_press) / q_b[0];

  // multiple timestep iteration
  
  for (i = 0; i < t_iter; i++) {
    for (j = 0; j < t_order; j++) {

      // update thermostat velocities a half step
  
      eta_dot_b[p_chain-1] += wdti2[j] * f_eta_b[p_chain-1];
  
      for (k = 1; k < p_chain; k++) {
        tmp = wdti4[j] * eta_dot_b[p_chain-k];
        ms = maclaurin_series(tmp);
        s = exp(-0.5 * tmp);
        s2 = s * s;
        eta_dot_b[p_chain-k-1] = eta_dot_b[p_chain-k-1] * s2 + 
          wdti2[j] * f_eta_b[p_chain-k-1] * s * ms;
      }
  
      // update thermostat positions
  
      for (k = 0; k < p_chain; k++)
        eta_b[k] += wdti1[j] * eta_dot_b[k];
  
      // update thermostat forces
  
      for (k = 1; k < p_chain; k++) {
        f_eta_b[k] = q_b[k-1] * eta_dot_b[k-1] * eta_dot_b[k-1] - kt;
        f_eta_b[k] /= q_b[k];
      }

      // update thermostat velocites a full step
  
      for (k = 0; k < p_chain-1; k++) {
        tmp = wdti4[j] * eta_dot_b[k+1];
        ms = maclaurin_series(tmp);
        s = exp(-0.5 * tmp);
        s2 = s * s;
        eta_dot_b[k] = eta_dot_b[k] * s2 + wdti2[j] * f_eta_b[k] * s * ms;
        tmp = q_b[k] * eta_dot_b[k] * eta_dot_b[k] - kt;
        f_eta_b[k+1] = tmp / q_b[k+1];
      }
  
      eta_dot_b[p_chain-1] += wdti2[j] * f_eta_b[p_chain-1];
    }
  }
}

/* ---------------------------------------------------------------------- 
   compute kinetic energy in the extended Hamiltonian
   conserved quantity = sum of returned energy and potential energy
-----------------------------------------------------------------------*/

double FixRigidNH::compute_scalar()
{
  const double kt = boltz * t_target;
  double energy;
  int i;
  
  energy = FixRigid::compute_scalar();
  
  if (tstat_flag) {
  
    // thermostat chain energy: from equation 12 in Kameraj et al (JCP 2005)

    energy += kt * (nf_t * eta_t[0] + nf_r * eta_r[0]);
  
    for (i = 1; i < t_chain; i++) 
      energy += kt * (eta_t[i] + eta_r[i]);
  
    for (i = 0;  i < t_chain; i++) {
      energy += 0.5 * q_t[i] * (eta_dot_t[i] * eta_dot_t[i]);
      energy += 0.5 * q_r[i] * (eta_dot_r[i] * eta_dot_r[i]);
    }
  }
  
  if (pstat_flag) {

    // using equation 22 in Kameraj et al for H_NPT

    double e = 0.0;
    for (i = 0; i < 3; i++)
      if (p_flag[i])
        e += epsilon_mass[i] * epsilon_dot[i] * epsilon_dot[i];
    energy += e*(0.5/pdim);
  
    double vol;
    if (dimension == 2) vol = domain->xprd * domain->yprd;
    else vol = domain->xprd * domain->yprd * domain->zprd;

    double p0 = (p_target[0] + p_target[1] + p_target[2]) / 3.0;
    energy += p0 * vol / nktv2p;

    for (i = 0;  i < p_chain; i++) {
      energy += kt * eta_b[i];
      energy += 0.5 * q_b[i] * (eta_dot_b[i] * eta_dot_b[i]);
    }
  }
  
  return energy;
}

/* ---------------------------------------------------------------------- */

void FixRigidNH::couple()
{
  double *tensor = pressure->vector;
 
  if (pstyle == ISO) {
    p_current[0] = p_current[1] = p_current[2] = pressure->scalar;
  } else if (pcouple == XYZ) {
    double ave = 1.0/3.0 * (tensor[0] + tensor[1] + tensor[2]);
    p_current[0] = p_current[1] = p_current[2] = ave;
  } else if (pcouple == XY) {
    double ave = 0.5 * (tensor[0] + tensor[1]);
    p_current[0] = p_current[1] = ave;
    p_current[2] = tensor[2];
  } else if (pcouple == YZ) {
    double ave = 0.5 * (tensor[1] + tensor[2]);
    p_current[1] = p_current[2] = ave;
    p_current[0] = tensor[0];
  } else if (pcouple == XZ) {
    double ave = 0.5 * (tensor[0] + tensor[2]);
    p_current[0] = p_current[2] = ave;
    p_current[1] = tensor[1];
  } else {
    p_current[0] = tensor[0];
    p_current[1] = tensor[1];
    p_current[2] = tensor[2];
  }
}

/* ---------------------------------------------------------------------- */

void FixRigidNH::remap()
{
  int i;
  double oldlo,oldhi,ctr,expfac;
  
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  // epsilon is not used, except for book-keeping
  
  for (i = 0; i < 3; i++) epsilon[i] += dtq * epsilon_dot[i];
  
  // convert pertinent atoms and rigid bodies to lamda coords
  
  if (allremap) domain->x2lamda(nlocal);
  else {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & dilate_group_bit)
        domain->x2lamda(x[i],x[i]);
  }
  
  if (nrigidfix)
    for (i = 0; i < nrigidfix; i++)
      modify->fix[rfix[i]]->deform(0);
  
  // reset global and local box to new size/shape
  
  for (i = 0; i < 3; i++) {
    if (p_flag[i]) {
      oldlo = domain->boxlo[i];
      oldhi = domain->boxhi[i];
      ctr = 0.5 * (oldlo + oldhi);
      expfac = exp(dtq * epsilon_dot[i]);
      domain->boxlo[i] = (oldlo-ctr)*expfac + ctr;
      domain->boxhi[i] = (oldhi-ctr)*expfac + ctr;
    }
  }

  domain->set_global_box();
  domain->set_local_box();
  
  // convert pertinent atoms and rigid bodies back to box coords
  
  if (allremap) domain->lamda2x(nlocal);
  else {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & dilate_group_bit)
        domain->lamda2x(x[i],x[i]);
  }
  
  if (nrigidfix)
    for (i = 0; i< nrigidfix; i++)
      modify->fix[rfix[i]]->deform(1);
}

/* ----------------------------------------------------------------------
   compute target temperature and kinetic energy
-----------------------------------------------------------------------*/

void FixRigidNH::compute_temp_target()
{
  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;
      
  t_target = t_start + delta * (t_stop-t_start);
}

/* ----------------------------------------------------------------------
   compute hydrostatic target pressure
-----------------------------------------------------------------------*/

void FixRigidNH::compute_press_target()
{
  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;
      
  p_hydro = 0.0;
  for (int i = 0; i < 3; i++)
    if (p_flag[i]) {
      p_target[i] = p_start[i] + delta * (p_stop[i]-p_start[i]);
      p_hydro += p_target[i];
    }
  p_hydro /= pdim;
}

/* ----------------------------------------------------------------------
   update epsilon_dot
-----------------------------------------------------------------------*/

void FixRigidNH::nh_epsilon_dot()
{
  int i;
  double volume,scale,f_epsilon;

  if (dimension == 2) volume = domain->xprd*domain->yprd;
  else volume = domain->xprd*domain->yprd*domain->zprd;

  // MTK terms
  
  mtk_term1 = (akin_t + akin_r) * mvv2e / g_f;
  
  scale = exp(-1.0 * dtq * eta_dot_b[0]);

  for (i = 0; i < 3; i++)
    if (p_flag[i]) {
      f_epsilon = (p_current[i]-p_hydro)*volume / nktv2p + mtk_term1;
      f_epsilon /= epsilon_mass[i];
      epsilon_dot[i] += dtq * f_epsilon;
      epsilon_dot[i] *= scale;
    }
  
  mtk_term2 = 0.0;
  for (i = 0; i < 3; i++)
    if (p_flag[i]) mtk_term2 += epsilon_dot[i];
  mtk_term2 /= g_f;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write 
------------------------------------------------------------------------- */

void FixRigidNH::write_restart(FILE *fp)
{
  if (tstat_flag == 0 && pstat_flag == 0) return;
  
  int nsize = 2; // tstat_flag and pstat_flag
  
  if (tstat_flag) {
    nsize += 1;         // t_chain
    nsize += 4*t_chain; // eta_t, eta_r, eta_dot_t, eta_dot_r
  }

  if (pstat_flag) {
    nsize += 7;         // p_chain, epsilon(3) and epsilon_dot(3)
    nsize += 2*p_chain;
  }
  
  double *list;
  memory->create(list,nsize,"rigid_nh:list");
  
  int n = 0;
  
  list[n++] = tstat_flag;
  if (tstat_flag) {
    list[n++] = t_chain;
    for (int i = 0; i < t_chain; i++) {
      list[n++] = eta_t[i];
      list[n++] = eta_r[i];
      list[n++] = eta_dot_t[i];
      list[n++] = eta_dot_r[i];
    }
  }
  
  list[n++] = pstat_flag;
  if (pstat_flag) {
    list[n++] = epsilon[0];
    list[n++] = epsilon[1];
    list[n++] = epsilon[2];
    list[n++] = epsilon_dot[0];
    list[n++] = epsilon_dot[1];
    list[n++] = epsilon_dot[2];
        
    list[n++] = p_chain;
    for (int i = 0; i < p_chain; i++) {
      list[n++] = eta_b[i];
      list[n++] = eta_dot_b[i];
    }
  }
  
  if (comm->me == 0) {
    int size = (nsize)*sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),nsize,fp);
  }
    
  memory->destroy(list);
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix 
------------------------------------------------------------------------- */

void FixRigidNH::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;
  int flag = static_cast<int> (list[n++]);
  
  if (flag) {
    int m = static_cast<int> (list[n++]);
    if (tstat_flag && m == t_chain) {
      for (int i = 0; i < t_chain; i++) {
        eta_t[i] = list[n++];
        eta_r[i] = list[n++];
        eta_dot_t[i] = list[n++];
        eta_dot_r[i] = list[n++];
      }
    } else n += 4*m;
  }

  flag = static_cast<int> (list[n++]);
  if (flag) {
    epsilon[0] = list[n++];
    epsilon[1] = list[n++];
    epsilon[2] = list[n++];
    epsilon_dot[0] = list[n++];
    epsilon_dot[1] = list[n++];
    epsilon_dot[2] = list[n++];
    
    int m = static_cast<int> (list[n++]);
    if (pstat_flag && m == p_chain) {
      for (int i = 0; i < p_chain; i++) {
        eta_b[i] = list[n++];
        eta_dot_b[i] = list[n++];
      }
    } else n += 2*m;
  }
}

/* ---------------------------------------------------------------------- */

int FixRigidNH::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (!pstat_flag) error->all(FLERR,"Illegal fix_modify command");
    if (tcomputeflag) {
      modify->delete_compute(id_temp);
      tcomputeflag = 0;
    }
    delete [] id_temp;
    int n = strlen(arg[1]) + 1;
    id_temp = new char[n];
    strcpy(id_temp,arg[1]);

    int icompute = modify->find_compute(arg[1]);
    if (icompute < 0)
      error->all(FLERR,"Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR,
                 "Fix_modify temperature ID does not compute temperature");
    if (temperature->igroup != 0 && comm->me == 0)
      error->warning(FLERR,"Temperature for fix modify is not for group all");

    // reset id_temp of pressure to new temperature ID

    if (pstat_flag) {
      icompute = modify->find_compute(id_press);
      if (icompute < 0)
        error->all(FLERR,"Pressure ID for fix modify does not exist");
      modify->compute[icompute]->reset_extra_compute_fix(id_temp);
    }

    return 2;

  } else if (strcmp(arg[0],"press") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (!pstat_flag) error->all(FLERR,"Illegal fix_modify command");
    if (pcomputeflag) {
      modify->delete_compute(id_press);
      pcomputeflag = 0;
    }
    delete [] id_press;
    int n = strlen(arg[1]) + 1;
    id_press = new char[n];
    strcpy(id_press,arg[1]);

    int icompute = modify->find_compute(arg[1]);
    if (icompute < 0) error->all(FLERR,"Could not find fix_modify pressure ID");
    pressure = modify->compute[icompute];

    if (pressure->pressflag == 0)
      error->all(FLERR,"Fix_modify pressure ID does not compute pressure");
    return 2;
  }

  return 0;
}

/* ---------------------------------------------------------------------- */

void FixRigidNH::allocate_chain()
{
  if (tstat_flag) {
    q_t = new double[t_chain];
    q_r = new double[t_chain];
    eta_t = new double[t_chain];
    eta_r = new double[t_chain];
    eta_dot_t = new double[t_chain];
    eta_dot_r = new double[t_chain];
    f_eta_t = new double[t_chain];
    f_eta_r = new double[t_chain];
  }
  
  if (pstat_flag) {
    q_b = new double[p_chain];
    eta_b = new double[p_chain];
    eta_dot_b = new double[p_chain];
    f_eta_b = new double[p_chain];
  }
}

/* ---------------------------------------------------------------------- */

void FixRigidNH::reset_target(double t_new)
{
  t_start = t_stop = t_new;
}

/* ---------------------------------------------------------------------- */

void FixRigidNH::allocate_order()
{
  w = new double[t_order];
  wdti1 = new double[t_order];
  wdti2 = new double[t_order];
  wdti4 = new double[t_order];
}

/* ---------------------------------------------------------------------- */

void FixRigidNH::deallocate_chain()
{
  if (tstat_flag) {
    delete [] q_t;
    delete [] q_r;
    delete [] eta_t;
    delete [] eta_r;
    delete [] eta_dot_t;
    delete [] eta_dot_r;
    delete [] f_eta_t;
    delete [] f_eta_r;
  }
  
  if (pstat_flag) {
    delete [] q_b;
    delete [] eta_b;
    delete [] eta_dot_b;
    delete [] f_eta_b;
  }
}

/* ---------------------------------------------------------------------- */

void FixRigidNH::deallocate_order()
{
  delete [] w;
  delete [] wdti1;
  delete [] wdti2;
  delete [] wdti4;
}

