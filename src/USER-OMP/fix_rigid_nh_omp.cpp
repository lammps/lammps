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

#include "fix_rigid_nh_omp.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "force.h"
#include "kspace.h"
#include "modify.h"
#include "update.h"

#include <string.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "math_extra.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace MathExtra;

enum{SINGLE,MOLECULE,GROUP};	// same as in FixRigid
enum{ISO,ANISO,TRICLINIC};	// same as in FixRigid

#define EINERTIA 0.4            // moment of inertia prefactor for ellipsoid

typedef struct { double x,y,z; } dbl3_t;

/* ----------------------------------------------------------------------
   perform preforce velocity Verlet integration
   see Kamberaj paper for step references
------------------------------------------------------------------------- */

void FixRigidNHOMP::initial_integrate(int vflag)
{
  double scale_r,scale_t[3],scale_v[3];
  
  // compute target temperature
  // update thermostat chains coupled to particles
  
  if (tstat_flag) {
    compute_temp_target();
    nhc_temp_integrate();
  }

  // compute target pressure
  // update epsilon dot
  // update thermostat coupled to barostat
  
  if (pstat_flag) {
    nhc_press_integrate();
    
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
  
  // compute scale variables

  scale_t[0] = scale_t[1] = scale_t[2] = 1.0;
  scale_v[0] = scale_v[1] = scale_v[2] = 1.0;
  scale_r = 1.0;

  if (tstat_flag) {
    akin_t = akin_r = 0.0;
    double tmp = exp(-dtq * eta_dot_t[0]);
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

    double tmp = dtq * epsilon_dot[0];
    scale_v[0] = dtv * exp(tmp) * maclaurin_series(tmp);
    tmp = dtq * epsilon_dot[1];
    scale_v[1] = dtv * exp(tmp) * maclaurin_series(tmp);
    tmp = dtq * epsilon_dot[2];
    scale_v[2] = dtv * exp(tmp) * maclaurin_series(tmp);
  }
    
  // update xcm, vcm, quat, conjqm and angmom
  double akt=0.0, akr=0.0;
  int ibody;

#if defined(_OPENMP)
#pragma omp parallel for default(none) private(ibody) shared(scale_r,scale_t,scale_v) schedule(static) reduction(+:akt,akr)
#endif
  for (ibody = 0; ibody < nbody; ibody++) {
    double mbody[3],tbody[3],fquat[4];
    const double dtf2 = dtf * 2.0;
    
    // step 1.1 - update vcm by 1/2 step
    
    const double dtfm = dtf / masstotal[ibody];
    vcm[ibody][0] += dtfm * fcm[ibody][0] * fflag[ibody][0];
    vcm[ibody][1] += dtfm * fcm[ibody][1] * fflag[ibody][1];
    vcm[ibody][2] += dtfm * fcm[ibody][2] * fflag[ibody][2];
    
    if (tstat_flag || pstat_flag) {
      vcm[ibody][0] *= scale_t[0];
      vcm[ibody][1] *= scale_t[1];
      vcm[ibody][2] *= scale_t[2];
      
      double tmp = vcm[ibody][0]*vcm[ibody][0] + vcm[ibody][1]*vcm[ibody][1] +
        vcm[ibody][2]*vcm[ibody][2];
      akt += masstotal[ibody]*tmp;
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
    
    if (tstat_flag || pstat_flag) {
      akr += angmom[ibody][0]*omega[ibody][0] + 
        angmom[ibody][1]*omega[ibody][1] + angmom[ibody][2]*omega[ibody][2];
    }
  } // end of parallel for

  if (pstat_flag || tstat_flag) {
    akin_t = akt;
    akin_r = akr;
  }

  // virial setup before call to set_xv

  if (vflag) v_setup(vflag);
  else evflag = 0;
  
  // remap simulation box by 1/2 step

  if (pstat_flag) remap();
  
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

  // remap simulation box by full step
  // redo KSpace coeffs since volume has changed

  if (pstat_flag) {
    remap();
    if (kspace_flag) force->kspace->setup();
  }  
}

/* ---------------------------------------------------------------------- */

void FixRigidNHOMP::final_integrate()
{
  double scale_t[3],scale_r;

  // compute scale variables
  
  scale_t[0] = scale_t[1] = scale_t[2] = 1.0;
  scale_r = 1.0;

  if (tstat_flag) {
    double tmp = exp(-1.0 * dtq * eta_dot_t[0]);
    scale_t[0] = scale_t[1] = scale_t[2] = tmp;
    scale_r = exp(-1.0 * dtq * eta_dot_r[0]);
  } 
  
  if (pstat_flag) {
    scale_t[0] *= exp(-dtq * (epsilon_dot[0] + mtk_term2));
    scale_t[1] *= exp(-dtq * (epsilon_dot[1] + mtk_term2));
    scale_t[2] *= exp(-dtq * (epsilon_dot[2] + mtk_term2));
    scale_r *= exp(-dtq * (pdim * mtk_term2));
    
    akin_t = akin_r = 0.0;
  }
  
  double * const * _noalias const x = atom->x;
  const dbl3_t * _noalias const f = (dbl3_t *) atom->f[0];
  const double * const * const torque_one = atom->torque;
  const imageint * _noalias const image = atom->image;
  const int nlocal = atom->nlocal;

  // sum over atoms to get force and torque on rigid body
  // we have 3 different strategies for multi-threading this.

   if (rstyle == SINGLE) {
     // we have just one rigid body. use OpenMP reduction to get sum[]
     double s0=0.0,s1=0.0,s2=0.0,s3=0.0,s4=0.0,s5=0.0;
     int i;

#if defined(_OPENMP)
#pragma omp parallel for default(none) private(i) reduction(+:s0,s1,s2,s3,s4,s5)
#endif
     for (i = 0; i < nlocal; i++) {
       const int ibody = body[i];
       if (ibody < 0) continue;

       double unwrap[3];
       domain->unmap(x[i],image[i],unwrap);
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
       int i;

#if defined(_OPENMP)
#pragma omp parallel for default(none) private(i) shared(ib) reduction(+:s0,s1,s2,s3,s4,s5)
#endif
       for (i = 0; i < nlocal; i++) {
	 const int ibody = body[i];
	 if (ibody != ib) continue;

	 s0 += f[i].x;
	 s1 += f[i].y;
	 s2 += f[i].z;

	 double unwrap[3];
	 domain->unmap(x[i],image[i],unwrap);
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
	 domain->unmap(x[i],image[i],unwrap);
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
  double akt=0.0,akr=0.0;
  const double dtf2 = dtf * 2.0;
  int ibody;

#if defined(_OPENMP)
#pragma omp parallel for default(none) private(ibody) shared(scale_t,scale_r) schedule(static) reduction(+:akt,akr)
#endif
  for (ibody = 0; ibody < nbody; ibody++) {
    double mbody[3],tbody[3],fquat[4];

    fcm[ibody][0] = all[ibody][0] + langextra[ibody][0];
    fcm[ibody][1] = all[ibody][1] + langextra[ibody][1];
    fcm[ibody][2] = all[ibody][2] + langextra[ibody][2];
    torque[ibody][0] = all[ibody][3] + langextra[ibody][3];
    torque[ibody][1] = all[ibody][4] + langextra[ibody][4];
    torque[ibody][2] = all[ibody][5] + langextra[ibody][5];

    // update vcm by 1/2 step

    const double dtfm = dtf / masstotal[ibody];
    if (tstat_flag || pstat_flag) {
      vcm[ibody][0] *= scale_t[0];
      vcm[ibody][1] *= scale_t[1];
      vcm[ibody][2] *= scale_t[2];
    }
    
    vcm[ibody][0] += dtfm * fcm[ibody][0] * fflag[ibody][0];
    vcm[ibody][1] += dtfm * fcm[ibody][1] * fflag[ibody][1];
    vcm[ibody][2] += dtfm * fcm[ibody][2] * fflag[ibody][2];
    
    if (pstat_flag) {
      double tmp = vcm[ibody][0]*vcm[ibody][0] + vcm[ibody][1]*vcm[ibody][1] +
        vcm[ibody][2]*vcm[ibody][2];
      akt += masstotal[ibody]*tmp;
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
      akr += angmom[ibody][0]*omega[ibody][0] + 
        angmom[ibody][1]*omega[ibody][1] + 
        angmom[ibody][2]*omega[ibody][2];
    }
  } // end of parallel for
  if (pstat_flag) {
    akin_t += akt;
    akin_r += akr;
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
  
  // compute temperature and pressure tensor
  // couple to compute current pressure components
  // trigger virial computation on next timestep
  
  if (tcomputeflag) t_current = temperature->compute_scalar();
  if (pstat_flag) {
    if (pstyle == ISO) pressure->compute_scalar();
    else pressure->compute_vector();
    couple();
    pressure->addstep(update->ntimestep+1);
  }

  if (pstat_flag) nh_epsilon_dot();  
  
  // update eta_dot_t and eta_dot_r
  // update eta_dot_b
      
  if (tstat_flag) nhc_temp_integrate();
  if (pstat_flag) nhc_press_integrate();  
}

/* ---------------------------------------------------------------------- */

void FixRigidNHOMP::remap()
{
  double * const * _noalias const x = atom->x;
  const int * _noalias const mask = atom->mask;
  const int nlocal = atom->nlocal;

  // epsilon is not used, except for book-keeping
  
  for (int i = 0; i < 3; i++) epsilon[i] += dtq * epsilon_dot[i];
  
  // convert pertinent atoms and rigid bodies to lamda coords
  
  if (allremap) domain->x2lamda(nlocal);
  else {
    int i;
#if defined (_OPENMP)
#pragma omp parallel for private(i) default(none) schedule(static)
#endif
    for (i = 0; i < nlocal; i++)
      if (mask[i] & dilate_group_bit)
        domain->x2lamda(x[i],x[i]);
  }
  
  if (nrigid)
    for (int i = 0; i < nrigidfix; i++)
      modify->fix[rfix[i]]->deform(0);
  
  // reset global and local box to new size/shape
  
  for (int i = 0; i < 3; i++) {
    if (p_flag[i]) {
      const double oldlo = domain->boxlo[i];
      const double oldhi = domain->boxhi[i];
      const double ctr = 0.5 * (oldlo + oldhi);
      const double expfac = exp(dtq * epsilon_dot[i]);
      domain->boxlo[i] = (oldlo-ctr)*expfac + ctr;
      domain->boxhi[i] = (oldhi-ctr)*expfac + ctr;
    }
  }

  domain->set_global_box();
  domain->set_local_box();
  
  // convert pertinent atoms and rigid bodies back to box coords
  
  if (allremap) domain->lamda2x(nlocal);
  else {
    int i;
#if defined (_OPENMP)
#pragma omp parallel for private(i) default(none) schedule(static)
#endif
    for (i = 0; i < nlocal; i++)
      if (mask[i] & dilate_group_bit)
        domain->lamda2x(x[i],x[i]);
  }

  if (nrigid)
    for (int i = 0; i< nrigidfix; i++)
      modify->fix[rfix[i]]->deform(1);
}

/* ----------------------------------------------------------------------
   set space-frame coords and velocity of each atom in each rigid body
   set orientation and rotation of extended particles
   x = Q displace + Xcm, mapped back to periodic box
   v = Vcm + (W cross (x - Xcm))

   NOTE: this needs to be kept in sync with FixRigidOMP
------------------------------------------------------------------------- */
template <int TRICLINIC, int EVFLAG>
void FixRigidNHOMP::set_xv_thr()
{
  dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const v = (dbl3_t *) atom->v[0];
  const dbl3_t * _noalias const f = (dbl3_t *) atom->f[0];
  const double * _noalias const rmass = atom->rmass;
  const double * _noalias const mass = atom->mass;
  const int * _noalias const type = atom->type;
  const imageint * _noalias const image = atom->image;

  double v0=0.0,v1=0.0,v2=0.0,v3=0.0,v4=0.0,v5=0.0;

  const double xprd = domain->xprd;
  const double yprd = domain->yprd;
  const double zprd = domain->zprd;
  const double xy = domain->xy;
  const double xz = domain->xz;
  const double yz = domain->yz;

  // set x and v of each atom

  const int nlocal = atom->nlocal;
  int i;

#if defined(_OPENMP)
#pragma omp parallel for default(none) private(i) reduction(+:v0,v1,v2,v3,v4,v5)
#endif
  for (i = 0; i < nlocal; i++) {
    const int ibody = body[i];
    if (ibody < 0) continue;

    const dbl3_t &xcmi = * ((dbl3_t *) xcm[ibody]);
    const dbl3_t &vcmi = * ((dbl3_t *) vcm[ibody]);
    const dbl3_t &omegai = * ((dbl3_t *) omega[ibody]);

    const int xbox = (image[i] & IMGMASK) - IMGMAX;
    const int ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    const int zbox = (image[i] >> IMG2BITS) - IMGMAX;
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
        while (theta <= MINUSPI) theta += TWOPI;
        while (theta > MY_PI) theta -= TWOPI;
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

   NOTE: this needs to be kept in sync with FixRigidOMP
------------------------------------------------------------------------- */
template <int TRICLINIC, int EVFLAG>
void FixRigidNHOMP::set_v_thr()
{
  dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const v = (dbl3_t *) atom->v[0];
  const dbl3_t * _noalias const f = (dbl3_t *) atom->f[0];
  const double * _noalias const rmass = atom->rmass;
  const double * _noalias const mass = atom->mass;
  const int * _noalias const type = atom->type;
  const imageint * _noalias const image = atom->image;

  const double xprd = domain->xprd;
  const double yprd = domain->yprd;
  const double zprd = domain->zprd;
  const double xy = domain->xy;
  const double xz = domain->xz;
  const double yz = domain->yz;

  double v0=0.0,v1=0.0,v2=0.0,v3=0.0,v4=0.0,v5=0.0;

  // set v of each atom

  const int nlocal = atom->nlocal;
  int i;

#if defined(_OPENMP)
#pragma omp parallel for default(none) private(i) reduction(+:v0,v1,v2,v3,v4,v5)
#endif
  for (i = 0; i < nlocal; i++) {
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

      const int xbox = (image[i] & IMGMASK) - IMGMAX;
      const int ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      const int zbox = (image[i] >> IMG2BITS) - IMGMAX;
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
