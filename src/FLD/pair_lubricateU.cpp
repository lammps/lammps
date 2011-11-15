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
   Contributing authors: Amit Kumar and Michael Bybee (UIUC)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_lubricateU.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "domain.h"
#include "update.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define TOL 1E-4   // tolerance for conjugate gradient

/* ---------------------------------------------------------------------- */

PairLubricateU::PairLubricateU(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;

  // pair lubricateU cannot compute virial as F dot r
  // due to how drag forces are applied to atoms
  // correct method is how per-atom virial does it

  no_virial_fdotr_compute = 1;

  nmax = 0;
  fl = Tl = xl = NULL;

  cgmax = 0;
  bcg = xcg = rcg = rcg1 = pcg = RU =  NULL;

  // set comm size needed by this Pair

  comm_forward = 6;
}

/* ---------------------------------------------------------------------- */

PairLubricateU::~PairLubricateU()
{
  memory->destroy(fl);
  memory->destroy(Tl);
  memory->destroy(xl);

  memory->destroy(bcg);
  memory->destroy(xcg);
  memory->destroy(rcg);
  memory->destroy(rcg1);
  memory->destroy(pcg);
  memory->destroy(RU);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(cut_inner);
  }
}

/* ----------------------------------------------------------------------
   It first has to solve for the velocity of the particles such that
   the net force on the particles is zero. NOTE: it has to be the last
   type of pair interaction specified in the input file. Also, it
   assumes that no other types of interactions, like k-space, is
   present. As already mentioned, the net force on the particles after
   this pair interaction would be identically zero.
   ---------------------------------------------------------------------- */

void PairLubricateU::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype; 

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int nall = nlocal + nghost;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;
  
  // skip compute() if called from integrate::setup()
  // this is b/c do not want compute() to update velocities twice on a restart
  // when restarting, call compute on step N (last step of prev run),
  // again on step N (setup of restart run),
  // then on step N+1 (first step of restart)
  // so this is one extra time which leads to bad dynamics

  if (update->setupflag) return;

  // grow per-atom arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(fl);
    memory->destroy(Tl);
    memory->destroy(xl);
    nmax = atom->nmax;
    memory->create(fl,nmax,3,"pair:fl");
    memory->create(Tl,nmax,3,"pair:Tl");
    memory->create(xl,nmax,3,"pair:xl");
  }

  // Added to implement midpoint integration scheme
  // Save force, torque found so far. Also save the positions

  for (i=0;i<nlocal+nghost;i++) {
    for (j=0;j<3;j++) {
      fl[i][j] = f[i][j];
      Tl[i][j] = torque[i][j];
      xl[i][j] = x[i][j];
    }
  }
  
  // Stage one of Midpoint method
  // Solve for velocities based on intial positions

  stage_one();
  
  // find positions at half the timestep and store in xl

  intermediates(nall,xl);
  
  // store back the saved forces and torques in original arrays

  for(i=0;i<nlocal+nghost;i++) {
    for(j=0;j<3;j++) {
      f[i][j] = fl[i][j];
      torque[i][j] = Tl[i][j];
    }
  }
  
  // stage two: this will give the final velocities

  stage_two(xl);  
}

/* ------------------------------------------------------------------------
   Stage one of midpoint method
------------------------------------------------------------------------- */

void PairLubricateU::stage_one()
{
  int i,j,ii,jj,inum,jnum,itype,jtype; 
  
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  int newton_pair = force->newton_pair;
  double vxmu2f = force->vxmu2f;
  double inv_inertia,mo_inertia;
  int *ilist;
  
  double radi;
  int nprocs = comm->nprocs;
  
  inum = list->inum;
  ilist = list->ilist;
  
  if (6*inum > cgmax) {
    memory->destroy(bcg);
    memory->destroy(xcg);
    memory->destroy(rcg);
    memory->destroy(rcg1);
    memory->destroy(pcg);
    memory->destroy(RU);
    cgmax = 6*inum;
    memory->create(bcg,cgmax,"pair:bcg");
    memory->create(xcg,cgmax,"pair:bcg");
    memory->create(rcg,cgmax,"pair:bcg");
    memory->create(rcg1,cgmax,"pair:bcg");
    memory->create(pcg,cgmax,"pair:bcg");
    memory->create(RU,cgmax,"pair:bcg");
  }

  double alpha,beta;
  double normi,error,normig;
  double send[2],recv[2],rcg_dot_rcg;
  
  // First compute R_FE*E  

  compute_RE();
  
  // Reverse communication of forces and torques to 
  // accumulate the net force on each of the particles
  
  if (newton_pair) comm->reverse_comm();
  
  // CONJUGATE GRADIENT
  // Find the right hand side= -ve of all forces/torques
  // b = 6*Npart in overall size
  
  for(ii = 0; ii < inum; ii++) {
    i = ilist[ii];  
    for (j = 0; j < 3; j++) {
      bcg[6*ii+j] = -f[i][j];
      bcg[6*ii+j+3] = -torque[i][j];  
    }
  } 
  
  // Start solving the equation : F^H = -F^P -F^B - F^H_{Ef}   
  // Store initial guess for velocity and angular-velocities/angular momentum
  // NOTE velocities and angular velocities are assumed relative to the fluid
  
  for (ii=0;ii<inum;ii++)
    for (j=0;j<3;j++) {      
      xcg[6*ii+j] = 0.0;
      xcg[6*ii+j+3] = 0.0;       
    }
  
  // Copy initial guess to the global arrays to be acted upon by R_{FU}
  // and returned by f and torque arrays
  
  copy_vec_uo(inum,xcg,v,omega);
  
  // set velocities for ghost particles
  
  comm->forward_comm_pair(this);
  
  // Find initial residual
  
  compute_RU();
  
  // reverse communication of forces and torques
  
  if (newton_pair) comm->reverse_comm();
  
  copy_uo_vec(inum,f,torque,RU);
  
  for (i=0;i<6*inum;i++)
    rcg[i] = bcg[i] - RU[i];
  
  // Set initial conjugate direction
  
  for (i=0;i<6*inum;i++)
    pcg[i] = rcg[i];
  
  // Find initial norm of the residual or norm of the RHS (either is fine)
  
  normi = dot_vec_vec(6*inum,bcg,bcg);
  
  MPI_Allreduce(&normi,&normig,1,MPI_DOUBLE,MPI_SUM,world);
  
  // Loop until convergence
  
  do {    
    // find R*p
    
    copy_vec_uo(inum,pcg,v,omega);
    
    // set velocities for ghost particles   
    
    comm->forward_comm_pair(this);
    
    compute_RU();
    
    // reverse communication of forces and torques
    
    if (newton_pair) comm->reverse_comm();
    
    
    copy_uo_vec(inum,f,torque,RU);  
    
    // Find alpha
    
    send[0] = dot_vec_vec(6*inum,rcg,rcg);
    send[1] = dot_vec_vec(6*inum,RU,pcg);
    
    MPI_Allreduce(send,recv,2,MPI_DOUBLE,MPI_SUM,world);
    
    alpha = recv[0]/recv[1];
    rcg_dot_rcg = recv[0];
    
    // Find new x
    
    for (i=0;i<6*inum;i++)
      xcg[i] = xcg[i] + alpha*pcg[i];
    
    // find new residual
    
    for (i=0;i<6*inum;i++)
      rcg1[i] = rcg[i] - alpha*RU[i];
    
    // find beta
    
    send[0] = dot_vec_vec(6*inum,rcg1,rcg1);
    
    MPI_Allreduce(send,recv,1,MPI_DOUBLE,MPI_SUM,world);
    
    beta = recv[0]/rcg_dot_rcg;
    
    // Find new conjugate direction
    
    for (i=0;i<6*inum;i++)
      pcg[i] = rcg1[i] + beta*pcg[i];
    
    for (i=0;i<6*inum;i++)
      rcg[i] = rcg1[i];
    
    // Find relative error
    
    error = sqrt(recv[0]/normig);    
    
  } while (error > TOL);  
  
  // update the final converged velocities in respective arrays
  
  copy_vec_uo(inum,xcg,v,omega);

  // set velocities for ghost particles  

  comm->forward_comm_pair(this);
  
  // Find actual particle's velocities from relative velocities 
  // Only non-zero component of fluid's vel : vx=gdot*y and wz=-gdot/2  

  for (ii=0;ii<inum;ii++) {
    i = ilist[ii];
    itype = type[i];
    radi = radius[i];
    
    v[i][0] = v[i][0] + gdot*x[i][1];
    omega[i][2] = omega[i][2] - gdot/2.0;
  }
}

/*---------------------------------------------------------------
  Finds the position of the particles at half the time step
----------------------------------------------------------------*/

void PairLubricateU::intermediates(int nall, double **xl)
{
  int i;
  double **x = atom->x;
  double **v = atom->v;
  double dtv = update->dt;
  
  for (i=0;i<nall;i++) {
    xl[i][0] = x[i][0] + 0.5*dtv*v[i][0];
    xl[i][1] = x[i][1] + 0.5*dtv*v[i][1];
    xl[i][2] = x[i][2] + 0.5*dtv*v[i][2];
  }
}

/* ------------------------------------------------------------------------
   Stage one of midpoint method
------------------------------------------------------------------------- */

void PairLubricateU::stage_two(double **x)
{
  int i,j,ii,jj,inum,jnum,itype,jtype; 
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  int newton_pair = force->newton_pair;
  double vxmu2f = force->vxmu2f;
  double inv_inertia,mo_inertia;
  int *ilist;
  
  double radi;
  int nprocs = comm->nprocs;
  
  inum = list->inum;
  ilist = list->ilist;
  
  double alpha,beta;
  double normi,error,normig;
  double send[2],recv[2],rcg_dot_rcg;
  
  // First compute R_FE*E  

  compute_RE(x);
  
  // Reverse communication of forces and torques to 
  // accumulate the net force on each of the particles
  
  if (newton_pair) comm->reverse_comm();
  
  // CONJUGATE GRADIENT
  // Find the right hand side= -ve of all forces/torques
  // b = 6*Npart in overall size
  
  for(ii = 0; ii < inum; ii++) {
    i = ilist[ii];  
    for (j = 0; j < 3; j++) {
      bcg[6*ii+j] = -f[i][j];
      bcg[6*ii+j+3] = -torque[i][j];  
    }
  } 
  
  // Start solving the equation : F^H = -F^P -F^B - F^H_{Ef}   
  // Store initial guess for velocity and angular-velocities/angular momentum
  // NOTE velocities and angular velocities are assumed relative to the fluid
  
  for (ii=0;ii<inum;ii++)
    for (j=0;j<3;j++) {      
      xcg[6*ii+j] = 0.0;
      xcg[6*ii+j+3] = 0.0;       
    }
  
  // Copy initial guess to the global arrays to be acted upon by R_{FU}
  // and returned by f and torque arrays
  
  copy_vec_uo(inum,xcg,v,omega);
  
  // set velocities for ghost particles
  
  comm->forward_comm_pair(this);
  
  // Find initial residual
  
  compute_RU(x);
  
  // reverse communication of forces and torques
  
  if (newton_pair) comm->reverse_comm();
  
  copy_uo_vec(inum,f,torque,RU);
  
  for (i=0;i<6*inum;i++)
    rcg[i] = bcg[i] - RU[i];
  
  // Set initial conjugate direction
  
  for (i=0;i<6*inum;i++)
    pcg[i] = rcg[i];
  
  // Find initial norm of the residual or norm of the RHS (either is fine)
  
  normi = dot_vec_vec(6*inum,bcg,bcg);
  
  MPI_Allreduce(&normi,&normig,1,MPI_DOUBLE,MPI_SUM,world);
  
  // Loop until convergence
  
  do {    
    // find R*p
    
    copy_vec_uo(inum,pcg,v,omega);
    
    // set velocities for ghost particles   
    
    comm->forward_comm_pair(this);
    
    compute_RU(x);
    
    // reverse communication of forces and torques
    
    if (newton_pair) comm->reverse_comm();
    
    copy_uo_vec(inum,f,torque,RU);  
    
    // Find alpha
    
    send[0] = dot_vec_vec(6*inum,rcg,rcg);
    send[1] = dot_vec_vec(6*inum,RU,pcg);
    
    MPI_Allreduce(send,recv,2,MPI_DOUBLE,MPI_SUM,world);
    
    alpha = recv[0]/recv[1];
    rcg_dot_rcg = recv[0];
    
    // Find new x
    
    for (i=0;i<6*inum;i++)
      xcg[i] = xcg[i] + alpha*pcg[i];
    
    // find new residual
    
    for (i=0;i<6*inum;i++)
      rcg1[i] = rcg[i] - alpha*RU[i];
    
    // find beta
    
    send[0] = dot_vec_vec(6*inum,rcg1,rcg1);
    
    MPI_Allreduce(send,recv,1,MPI_DOUBLE,MPI_SUM,world);
    
    beta = recv[0]/rcg_dot_rcg;
    
    // Find new conjugate direction
    
    for (i=0;i<6*inum;i++)
      pcg[i] = rcg1[i] + beta*pcg[i];
    
    for (i=0;i<6*inum;i++)
      rcg[i] = rcg1[i];
    
    // Find relative error
    
    error = sqrt(recv[0]/normig);    
    
  } while (error > TOL);  
  
  
  // update the final converged velocities in respective arrays
  
  copy_vec_uo(inum,xcg,v,omega);
  
  // set velocities for ghost particles
  
  comm->forward_comm_pair(this);
  
  // Compute the viscosity/pressure
  
  if (evflag) compute_Fh(x);
  
  // Find actual particle's velocities from relative velocities 
  // Only non-zero component of fluid's vel : vx=gdot*y and wz=-gdot/2
  
  for (ii=0;ii<inum;ii++) {
    i = ilist[ii];
    itype = type[i];
    radi = radius[i];
    
    v[i][0] = v[i][0] + gdot*x[i][1];
    omega[i][2] = omega[i][2] - gdot/2.0;
  }
}

/* ------------------------------------------------------------------------
   This function computes the final hydrodynamic force once the
   velocities have converged.
   ------------------------------------------------------------------------- */

void PairLubricateU::compute_Fh(double **x)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fpair,fx,fy,fz,tx,ty,tz;
  double rsq,r,h_sep,radi,tfmag;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3;
  double vt1,vt2,vt3,wdotn,wt1,wt2,wt3;   
  double inv_inertia;
  int *ilist,*jlist,*numneigh,**firstneigh;
  
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int newton_pair = force->newton_pair;

  double vxmu2f = force->vxmu2f;
  int overlaps = 0;
  double vi[3],vj[3],wi[3],wj[3],xl[3],a_sq,a_sh,a_pu,Fbmag,del,delmin,eta;
  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // Set force to zero which is the final value after this pair interaction

  for (i=0;i<nlocal+nghost;i++)
    for (j=0;j<3;j++) {
      f[i][j] = 0.0;
      torque[i][j] = 0.0;
    }
  
  // reverse communication of forces and torques

  if (newton_pair) comm->reverse_comm(); // not really needed
  
  // Find additional contribution from the stresslets

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];   
    
    // Find the contribution to stress from isotropic RS0
    // Set psuedo force to obtain the required contribution
    // need to set delx  and fy only

    fx = 0.0; delx = radi;
    fy = vxmu2f*RS0*gdot/2.0/radi; dely = 0.0;
    fz = 0.0; delz = 0.0;
    if (evflag)
      ev_tally_xyz(i,i,nlocal,newton_pair,0.0,0.0,-fx,-fy,-fz,delx,dely,delz);
    
    // Find angular velocity

    wi[0] = omega[i][0];
    wi[1] = omega[i][1];
    wi[2] = omega[i][2];          
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      
      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);  
        
        // Use omega directly if it exists, else angmom
        // angular momentum = I*omega = 2/5 * M*R^2 * omega
        
	wj[0] = omega[j][0];
	wj[1] = omega[j][1];
	wj[2] = omega[j][2];              
        
        // loc of the point of closest approach on particle i from its cente

        xl[0] = -delx/r*radi;
        xl[1] = -dely/r*radi;
        xl[2] = -delz/r*radi;
        
        // velocity at the point of closest approach on both particles
        // v = v + omega_cross_xl

        // particle i

        vi[0] = v[i][0] + (wi[1]*xl[2] - wi[2]*xl[1]);
        vi[1] = v[i][1] + (wi[2]*xl[0] - wi[0]*xl[2]);
        vi[2] = v[i][2] + (wi[0]*xl[1] - wi[1]*xl[0]);
        
        // particle j

        vj[0] = v[j][0] - (wj[1]*xl[2] - wj[2]*xl[1]);
        vj[1] = v[j][1] - (wj[2]*xl[0] - wj[0]*xl[2]);
        vj[2] = v[j][2] - (wj[0]*xl[1] - wj[1]*xl[0]);
        
        
        // Relative  velocity at the point of closest approach
	// include contribution from Einf of the fluid

        vr1 = vi[0] - vj[0] - 
	  2.0*(Ef[0][0]*xl[0] + Ef[0][1]*xl[1] + Ef[0][2]*xl[2]);
        vr2 = vi[1] - vj[1] - 
	  2.0*(Ef[1][0]*xl[0] + Ef[1][1]*xl[1] + Ef[1][2]*xl[2]);
        vr3 = vi[2] - vj[2] - 
	  2.0*(Ef[2][0]*xl[0] + Ef[2][1]*xl[1] + Ef[2][2]*xl[2]);
        
        // Normal component (vr.n)n

        vnnr = (vr1*delx + vr2*dely + vr3*delz)/r;
        vn1 = vnnr*delx/r;
        vn2 = vnnr*dely/r;
        vn3 = vnnr*delz/r;
        
        // Tangential component vr - (vr.n)n

        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;
        
        // Find the scalar resistances a_sq, a_sh and a_pu

        h_sep = r - 2.0*radi;
        
        // check for overlaps

        if (h_sep < 0.0) overlaps++;
        
        // If less than the minimum gap use the minimum gap instead

        if (r < cut_inner[itype][jtype])
          h_sep = cut_inner[itype][jtype] - 2.0*radi;          
        
        // Scale h_sep by radi

        h_sep = h_sep/radi;
        
        // Scalar resistances

        if (flaglog) {
          a_sq = 6.0*MY_PI*mu*radi*(1.0/4.0/h_sep + 9.0/40.0*log(1.0/h_sep));
          a_sh = 6.0*MY_PI*mu*radi*(1.0/6.0*log(1.0/h_sep));
        } else
          a_sq = 6.0*MY_PI*mu*radi*(1.0/4.0/h_sep);
                
        // Find force due to squeeze type motion

        fx  = a_sq*vn1;
        fy  = a_sq*vn2;
        fz  = a_sq*vn3;
        
        // Find force due to all shear kind of motions

        if (flaglog) {
          fx = fx + a_sh*vt1;
          fy = fy + a_sh*vt2;
          fz = fz + a_sh*vt3;                  
        }
        
        // Scale forces to obtain in appropriate units

        fx = vxmu2f*fx;
        fy = vxmu2f*fy;
        fz = vxmu2f*fz;
        
        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
				 0.0,0.0,-fx,-fy,-fz,delx,dely,delz);
      }
    }
  }
}

/* ---------------------------------------------------------------------- 
  computes R_FU * U 
---------------------------------------------------------------------- */
  
void PairLubricateU::compute_RU()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fpair,fx,fy,fz,tx,ty,tz;
  double rsq,r,h_sep,radi,tfmag;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3;
  double vt1,vt2,vt3,wdotn,wt1,wt2,wt3;   
  double inv_inertia;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int newton_pair = force->newton_pair;

  double vxmu2f = force->vxmu2f;
  int overlaps = 0;
  double vi[3],vj[3],wi[3],wj[3],xl[3],a_sq,a_sh,a_pu,Fbmag,del,delmin,eta;
  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // Initialize f to zero

  for (i=0;i<nlocal+nghost;i++)
    for (j=0;j<3;j++) {
      f[i][j] = 0.0;
      torque[i][j] = 0.0;
    }
  
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];   
    
    // Find angular velocity

    wi[0] = omega[i][0];
    wi[1] = omega[i][1];
    wi[2] = omega[i][2];          
     
    // Contribution due to the isotropic terms

    f[i][0] += -vxmu2f*R0*v[i][0];
    f[i][1] += -vxmu2f*R0*v[i][1];
    f[i][2] += -vxmu2f*R0*v[i][2];    
    
    torque[i][0] += -vxmu2f*RT0*wi[0];
    torque[i][1] += -vxmu2f*RT0*wi[1];
    torque[i][2] += -vxmu2f*RT0*wi[2];   
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);  
  
        // Use omega directly if it exists, else angmom
        // angular momentum = I*omega = 2/5 * M*R^2 * omega

	wj[0] = omega[j][0];
	wj[1] = omega[j][1];
	wj[2] = omega[j][2];              

        // loc of the point of closest approach on particle i from its center

        xl[0] = -delx/r*radi;
        xl[1] = -dely/r*radi;
        xl[2] = -delz/r*radi;
  
        // velocity at the point of closest approach on both particles
        // v = v + omega_cross_xl

        // particle i

        vi[0] = v[i][0] + (wi[1]*xl[2] - wi[2]*xl[1]);
        vi[1] = v[i][1] + (wi[2]*xl[0] - wi[0]*xl[2]);
        vi[2] = v[i][2] + (wi[0]*xl[1] - wi[1]*xl[0]);
  
        // particle j

        vj[0] = v[j][0] - (wj[1]*xl[2] - wj[2]*xl[1]);
        vj[1] = v[j][1] - (wj[2]*xl[0] - wj[0]*xl[2]);
        vj[2] = v[j][2] - (wj[0]*xl[1] - wj[1]*xl[0]);
  
        // Find the scalar resistances a_sq and a_sh

        h_sep = r - 2.0*radi;
        
        // check for overlaps

        if(h_sep < 0.0) overlaps++;
        
        // If less than the minimum gap use the minimum gap instead

        if (r < cut_inner[itype][jtype])
          h_sep = cut_inner[itype][jtype] - 2.0*radi;          
        
        // Scale h_sep by radi

        h_sep = h_sep/radi;
  
        // Scalar resistances

        if (flaglog) {
          a_sq = 6.0*MY_PI*mu*radi*(1.0/4.0/h_sep + 9.0/40.0*log(1.0/h_sep));
          a_sh = 6.0*MY_PI*mu*radi*(1.0/6.0*log(1.0/h_sep));
          a_pu = 8.0*MY_PI*mu*pow(radi,3)*(3.0/160.0*log(1.0/h_sep));
        } else
          a_sq = 6.0*MY_PI*mu*radi*(1.0/4.0/h_sep);
  
        // Relative  velocity at the point of closest approach

        vr1 = vi[0] - vj[0];
        vr2 = vi[1] - vj[1];
        vr3 = vi[2] - vj[2];

        // Normal component (vr.n)n

        vnnr = (vr1*delx + vr2*dely + vr3*delz)/r;
        vn1 = vnnr*delx/r;
        vn2 = vnnr*dely/r;
        vn3 = vnnr*delz/r;

        // Tangential component vr - (vr.n)n

        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;

        // Find force due to squeeze type motion

        fx  = a_sq*vn1;
        fy  = a_sq*vn2;
        fz  = a_sq*vn3;
        
        // Find force due to all shear kind of motions

        if (flaglog) {
          fx = fx + a_sh*vt1;
          fy = fy + a_sh*vt2;
          fz = fz + a_sh*vt3;                  
        }
        
        // Scale forces to obtain in appropriate units

        fx = vxmu2f*fx;
        fy = vxmu2f*fy;
        fz = vxmu2f*fz;
        
        // Add to the total forc

        f[i][0] -= fx;
        f[i][1] -= fy;
        f[i][2] -= fz;    
        
        if (newton_pair || j < nlocal) {
          f[j][0] += fx;
          f[j][1] += fy;
          f[j][2] += fz;    
        }
  
        // Find torque due to this force

        if (flaglog) {
          tx = xl[1]*fz - xl[2]*fy;
          ty = xl[2]*fx - xl[0]*fz;
          tz = xl[0]*fy - xl[1]*fx;                  
  
          // Why a scale factor ?

          torque[i][0] -= vxmu2f*tx;
          torque[i][1] -= vxmu2f*ty;
          torque[i][2] -= vxmu2f*tz;        
            
          if(newton_pair || j < nlocal) {
            torque[j][0] -= vxmu2f*tx;
            torque[j][1] -= vxmu2f*ty;
            torque[j][2] -= vxmu2f*tz;
          }
          
          // Torque due to a_pu

          wdotn = ((wi[0]-wj[0])*delx + 
		   (wi[1]-wj[1])*dely + (wi[2]-wj[2])*delz)/r;
          wt1 = (wi[0]-wj[0]) - wdotn*delx/r;
          wt2 = (wi[1]-wj[1]) - wdotn*dely/r;
          wt3 = (wi[2]-wj[2]) - wdotn*delz/r;
          
          tx = a_pu*wt1;
          ty = a_pu*wt2;
          tz = a_pu*wt3;
          
          // add to total

          torque[i][0] -= vxmu2f*tx;
          torque[i][1] -= vxmu2f*ty;
          torque[i][2] -= vxmu2f*tz;        
            
          if (newton_pair || j < nlocal) {
            torque[j][0] += vxmu2f*tx;
            torque[j][1] += vxmu2f*ty;
            torque[j][2] += vxmu2f*tz;
          }
        }      
      }
    }
  }
}

/* ---------------------------------------------------------------------- 
  computes R_FU * U 
---------------------------------------------------------------------- */

void PairLubricateU::compute_RU(double **x)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fpair,fx,fy,fz,tx,ty,tz;
  double rsq,r,h_sep,radi,tfmag;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3;
  double vt1,vt2,vt3,wdotn,wt1,wt2,wt3;   
  double inv_inertia;
  int *ilist,*jlist,*numneigh,**firstneigh;
  
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int newton_pair = force->newton_pair;

  double vxmu2f = force->vxmu2f;
  int overlaps = 0;
  double vi[3],vj[3],wi[3],wj[3],xl[3],a_sq,a_sh,a_pu,Fbmag,del,delmin,eta;
  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // Initialize f to zero
  
  for (i=0;i<nlocal+nghost;i++)
    for (j=0;j<3;j++) {
      f[i][j] = 0.0;
      torque[i][j] = 0.0;
    }
  
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];   
    
    // Find angular velocity
    
    wi[0] = omega[i][0];
    wi[1] = omega[i][1];
    wi[2] = omega[i][2];
    
    // Contribution due to the isotropic terms
    
    f[i][0] += -vxmu2f*R0*v[i][0];
    f[i][1] += -vxmu2f*R0*v[i][1];
    f[i][2] += -vxmu2f*R0*v[i][2];    
    
    torque[i][0] += -vxmu2f*RT0*wi[0];
    torque[i][1] += -vxmu2f*RT0*wi[1];
    torque[i][2] += -vxmu2f*RT0*wi[2];   
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      
      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);  
        
        // Use omega directly if it exists, else angmom
        // angular momentum = I*omega = 2/5 * M*R^2 * omega
        
	wj[0] = omega[j][0];
	wj[1] = omega[j][1];
	wj[2] = omega[j][2];              
        
        // loc of the point of closest approach on particle i from its center
        
        xl[0] = -delx/r*radi;
        xl[1] = -dely/r*radi;
        xl[2] = -delz/r*radi;
        
        // velocity at the point of closest approach on both particles
        // v = v + omega_cross_xl
        
        // particle i
        
        vi[0] = v[i][0] + (wi[1]*xl[2] - wi[2]*xl[1]);
        vi[1] = v[i][1] + (wi[2]*xl[0] - wi[0]*xl[2]);
        vi[2] = v[i][2] + (wi[0]*xl[1] - wi[1]*xl[0]);
        
        // particle j
        
        vj[0] = v[j][0] - (wj[1]*xl[2] - wj[2]*xl[1]);
        vj[1] = v[j][1] - (wj[2]*xl[0] - wj[0]*xl[2]);
        vj[2] = v[j][2] - (wj[0]*xl[1] - wj[1]*xl[0]);
        
        // Find the scalar resistances a_sq and a_sh
        
        h_sep = r - 2.0*radi;
        
        // check for overlaps
        
        if(h_sep < 0.0) overlaps++;
        
        // If less than the minimum gap use the minimum gap instead
        
        if (r < cut_inner[itype][jtype])
          h_sep = cut_inner[itype][jtype] - 2.0*radi;          
        
        // Scale h_sep by radi
        
        h_sep = h_sep/radi;
        
        // Scalar resistances
        
        if (flaglog) {
          a_sq = 6.0*MY_PI*mu*radi*(1.0/4.0/h_sep + 9.0/40.0*log(1.0/h_sep));
          a_sh = 6.0*MY_PI*mu*radi*(1.0/6.0*log(1.0/h_sep));
          a_pu = 8.0*MY_PI*mu*pow(radi,3)*(3.0/160.0*log(1.0/h_sep));
        } else
          a_sq = 6.0*MY_PI*mu*radi*(1.0/4.0/h_sep);
        
        // Relative  velocity at the point of closest approach
        
        vr1 = vi[0] - vj[0];
        vr2 = vi[1] - vj[1];
        vr3 = vi[2] - vj[2];
        
        // Normal component (vr.n)n
        
        vnnr = (vr1*delx + vr2*dely + vr3*delz)/r;
        vn1 = vnnr*delx/r;
        vn2 = vnnr*dely/r;
        vn3 = vnnr*delz/r;
        
        // Tangential component vr - (vr.n)n
        
        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;
        
        // Find force due to squeeze type motion
        
        fx  = a_sq*vn1;
        fy  = a_sq*vn2;
        fz  = a_sq*vn3;
        
        // Find force due to all shear kind of motions
        
        if (flaglog) {
          fx = fx + a_sh*vt1;
          fy = fy + a_sh*vt2;
          fz = fz + a_sh*vt3;                  
        }
        
        // Scale forces to obtain in appropriate units
        
        fx = vxmu2f*fx;
        fy = vxmu2f*fy;
        fz = vxmu2f*fz;
        
        // Add to the total force
        
        f[i][0] -= fx;
        f[i][1] -= fy;
        f[i][2] -= fz;    
        
        if (newton_pair || j < nlocal) {
          f[j][0] += fx;
          f[j][1] += fy;
          f[j][2] += fz;    
        }
        
        // Find torque due to this force
        
        if (flaglog) {
          tx = xl[1]*fz - xl[2]*fy;
          ty = xl[2]*fx - xl[0]*fz;
          tz = xl[0]*fy - xl[1]*fx;                  
          
          // Why a scale factor ?
          
          torque[i][0] -= vxmu2f*tx;
          torque[i][1] -= vxmu2f*ty;
          torque[i][2] -= vxmu2f*tz;        
          
          if(newton_pair || j < nlocal) {
            torque[j][0] -= vxmu2f*tx;
            torque[j][1] -= vxmu2f*ty;
            torque[j][2] -= vxmu2f*tz;
          }
          
          // Torque due to a_pu
          
          wdotn = ((wi[0]-wj[0])*delx + 
                   (wi[1]-wj[1])*dely + (wi[2]-wj[2])*delz)/r;
          wt1 = (wi[0]-wj[0]) - wdotn*delx/r;
          wt2 = (wi[1]-wj[1]) - wdotn*dely/r;
          wt3 = (wi[2]-wj[2]) - wdotn*delz/r;
          
          tx = a_pu*wt1;
          ty = a_pu*wt2;
          tz = a_pu*wt3;
          
          // add to total
          
          torque[i][0] -= vxmu2f*tx;
          torque[i][1] -= vxmu2f*ty;
          torque[i][2] -= vxmu2f*tz;        
          
          if (newton_pair || j < nlocal) {
            torque[j][0] += vxmu2f*tx;
            torque[j][1] += vxmu2f*ty;
            torque[j][2] += vxmu2f*tz;
          }
        }      
      }
    }
  }
}

/* ----------------------------------------------------------------------
   This computes R_{FE}*E , where E is the rate of strain of tensor which is 
   known apriori, as it depends only on the known fluid velocity.
   So, this part of the hydrodynamic interaction can be pre computed and
   transferred to the RHS
   ---------------------------------------------------------------------- */
  
void PairLubricateU::compute_RE()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fpair,fx,fy,fz,tx,ty,tz;
  double rsq,r,h_sep,radi,tfmag;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3;
  double vt1,vt2,vt3;  
  double inv_inertia;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int newton_pair = force->newton_pair;

  double vxmu2f = force->vxmu2f;
  int overlaps = 0;
  double vi[3],vj[3],wi[3],wj[3],xl[3],a_sq,a_sh,a_pu,Fbmag,del,delmin,eta;
  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];   
       
    // No contribution from isotropic terms due to E    
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);  
  
        // loc of the point of closest approach on particle i from its center

        xl[0] = -delx/r*radi;
        xl[1] = -dely/r*radi;
        xl[2] = -delz/r*radi;
  
        // Find the scalar resistances a_sq and a_sh

        h_sep = r - 2.0*radi;
        
        // check for overlaps

        if(h_sep < 0.0) overlaps++;
        
        // If less than the minimum gap use the minimum gap instead

        if (r < cut_inner[itype][jtype])
          h_sep = cut_inner[itype][jtype] - 2.0*radi;          
        
        // Scale h_sep by radi

        h_sep = h_sep/radi;
  
        // Scalar resistance for Squeeze type motions

        if (flaglog)
          a_sq = 6*MY_PI*mu*radi*(1.0/4.0/h_sep + 9.0/40.0*log(1/h_sep));
        else
          a_sq = 6*MY_PI*mu*radi*(1.0/4.0/h_sep);
        
        // Scalar resistance for Shear type motions

        if (flaglog) {
          a_sh = 6*MY_PI*mu*radi*(1.0/6.0*log(1/h_sep));
          a_pu = 8.0*MY_PI*mu*pow(radi,3)*(3.0/160.0*log(1.0/h_sep));          
        }
         
        // Relative velocity at the point of closest approach due to Ef only

        vr1 = -2.0*(Ef[0][0]*xl[0] + Ef[0][1]*xl[1] + Ef[0][2]*xl[2]);
        vr2 = -2.0*(Ef[1][0]*xl[0] + Ef[1][1]*xl[1] + Ef[1][2]*xl[2]);
        vr3 = -2.0*(Ef[2][0]*xl[0] + Ef[2][1]*xl[1] + Ef[2][2]*xl[2]);        
            
        // Normal component (vr.n)n

        vnnr = (vr1*delx + vr2*dely + vr3*delz)/r;
        vn1 = vnnr*delx/r;
        vn2 = vnnr*dely/r;
        vn3 = vnnr*delz/r;

        // Tangential component vr - (vr.n)n

        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;

        // Find force due to squeeze type motion

        fx  = a_sq*vn1;
        fy  = a_sq*vn2;
        fz  = a_sq*vn3;
        
        // Find force due to all shear kind of motions

        if (flaglog) {
          fx = fx + a_sh*vt1;
          fy = fy + a_sh*vt2;
          fz = fz + a_sh*vt3;                  
        }
        
        // Scale forces to obtain in appropriate units

        fx = vxmu2f*fx;
        fy = vxmu2f*fy;
        fz = vxmu2f*fz;
        
        // Add to the total forc

        f[i][0] -= fx;
        f[i][1] -= fy;
        f[i][2] -= fz;    
        
        if (newton_pair || j < nlocal) {
          f[j][0] += fx;
          f[j][1] += fy;
          f[j][2] += fz;    
        }
        
        // Find torque due to this force

        if (flaglog) {
          tx = xl[1]*fz - xl[2]*fy;
          ty = xl[2]*fx - xl[0]*fz;
          tz = xl[0]*fy - xl[1]*fx;                  
  
          // Why a scale factor ?

          torque[i][0] -= vxmu2f*tx;
          torque[i][1] -= vxmu2f*ty;
          torque[i][2] -= vxmu2f*tz;        
            
          if (newton_pair || j < nlocal) {
            torque[j][0] -= vxmu2f*tx;
            torque[j][1] -= vxmu2f*ty;
            torque[j][2] -= vxmu2f*tz;
          }
          
          // NOTE No a_pu term needed as they add up to zero
        }      
      }
    }
  }

  int print_overlaps = 0;
  if (print_overlaps && overlaps)
    printf("Number of overlaps=%d\n",overlaps);
}

/* ----------------------------------------------------------------------
   This computes R_{FE}*E , where E is the rate of strain of tensor which is 
   known apriori, as it depends only on the known fluid velocity.
   So, this part of the hydrodynamic interaction can be pre computed and
   transferred to the RHS
   ---------------------------------------------------------------------- */

void PairLubricateU::compute_RE(double **x)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fpair,fx,fy,fz,tx,ty,tz;
  double rsq,r,h_sep,radi,tfmag;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3;
  double vt1,vt2,vt3;  
  double inv_inertia;
  int *ilist,*jlist,*numneigh,**firstneigh;
  
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int newton_pair = force->newton_pair;

  double vxmu2f = force->vxmu2f;
  int overlaps = 0;
  double vi[3],vj[3],wi[3],wj[3],xl[3],a_sq,a_sh,a_pu,Fbmag,del,delmin,eta;
  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];   
    
    // No contribution from isotropic terms due to E    
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      
      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);  
        
        // loc of the point of closest approach on particle i from its center
        
        xl[0] = -delx/r*radi;
        xl[1] = -dely/r*radi;
        xl[2] = -delz/r*radi;
        
        // Find the scalar resistances a_sq and a_sh
        
        h_sep = r - 2.0*radi;
        
        // check for overlaps
        
        if(h_sep < 0.0) overlaps++;
        
        // If less than the minimum gap use the minimum gap instead
        
        if (r < cut_inner[itype][jtype])
          h_sep = cut_inner[itype][jtype] - 2.0*radi;          
        
        // Scale h_sep by radi
        
        h_sep = h_sep/radi;
        
        // Scalar resistance for Squeeze type motions
        
        if (flaglog)
          a_sq = 6*MY_PI*mu*radi*(1.0/4.0/h_sep + 9.0/40.0*log(1/h_sep));
        else
          a_sq = 6*MY_PI*mu*radi*(1.0/4.0/h_sep);
        
        // Scalar resistance for Shear type motions
        
        if (flaglog) {
          a_sh = 6*MY_PI*mu*radi*(1.0/6.0*log(1/h_sep));
          a_pu = 8.0*MY_PI*mu*pow(radi,3)*(3.0/160.0*log(1.0/h_sep));          
        }
        
        // Relative velocity at the point of closest approach due to Ef only
        
        vr1 = -2.0*(Ef[0][0]*xl[0] + Ef[0][1]*xl[1] + Ef[0][2]*xl[2]);
        vr2 = -2.0*(Ef[1][0]*xl[0] + Ef[1][1]*xl[1] + Ef[1][2]*xl[2]);
        vr3 = -2.0*(Ef[2][0]*xl[0] + Ef[2][1]*xl[1] + Ef[2][2]*xl[2]);        
        
        // Normal component (vr.n)n
        
        vnnr = (vr1*delx + vr2*dely + vr3*delz)/r;
        vn1 = vnnr*delx/r;
        vn2 = vnnr*dely/r;
        vn3 = vnnr*delz/r;
        
        // Tangential component vr - (vr.n)n
        
        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;
        
        // Find force due to squeeze type motion
        
        fx  = a_sq*vn1;
        fy  = a_sq*vn2;
        fz  = a_sq*vn3;
        
        // Find force due to all shear kind of motions
        
        if (flaglog) {
          fx = fx + a_sh*vt1;
          fy = fy + a_sh*vt2;
          fz = fz + a_sh*vt3;                  
        }
        
        // Scale forces to obtain in appropriate units
        
        fx = vxmu2f*fx;
        fy = vxmu2f*fy;
        fz = vxmu2f*fz;
        
        // Add to the total forc
        
        f[i][0] -= fx;
        f[i][1] -= fy;
        f[i][2] -= fz;    
        
        if (newton_pair || j < nlocal) {
          f[j][0] += fx;
          f[j][1] += fy;
          f[j][2] += fz;    
        }
        
        // Find torque due to this force
        
        if (flaglog) {
          tx = xl[1]*fz - xl[2]*fy;
          ty = xl[2]*fx - xl[0]*fz;
          tz = xl[0]*fy - xl[1]*fx;                  
          
          // Why a scale factor ?
          
          torque[i][0] -= vxmu2f*tx;
          torque[i][1] -= vxmu2f*ty;
          torque[i][2] -= vxmu2f*tz;        
          
          if (newton_pair || j < nlocal) {
            torque[j][0] -= vxmu2f*tx;
            torque[j][1] -= vxmu2f*ty;
            torque[j][2] -= vxmu2f*tz;
          }
          
          // NOTE No a_pu term needed as they add up to zero
        }      
      }
    }
  }
  
  int print_overlaps = 0;
  if (print_overlaps && overlaps)
    printf("Number of overlaps=%d\n",overlaps);
}


/* ----------------------------------------------------------------------
   allocate all arrays 
------------------------------------------------------------------------- */

void PairLubricateU::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  setflag = memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  cutsq = memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(cut_inner,n+1,n+1,"pair:cut_inner");
}

/*-----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void PairLubricateU::settings(int narg, char **arg)
{
  if (narg != 5) error->all(FLERR,"Illegal pair_style command");

  mu = atof(arg[0]);
  flaglog = atoi(arg[1]);
  cut_inner_global = atof(arg[2]);
  cut_global = atof(arg[3]);
  gdot =  atof(arg[4]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_inner[i][j] = cut_inner_global;
          cut[i][j] = cut_global;
        }
  }
  
  // store the rate of strain tensor

  Ef[0][0] = 0.0;
  Ef[0][1] = 0.5*gdot;
  Ef[0][2] = 0.0;
  Ef[1][0] = 0.5*gdot;
  Ef[1][1] = 0.0;
  Ef[1][2] = 0.0;
  Ef[2][0] = 0.0;
  Ef[2][1] = 0.0;
  Ef[2][2] = 0.0; 
}

/*-----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLubricateU::coeff(int narg, char **arg)
{
  if (narg != 2 && narg != 4)
    error->all(FLERR,"Incorrect args for pair coefficients");

  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double cut_inner_one = cut_inner_global;
  double cut_one = cut_global;
  if (narg == 4) {
    cut_inner_one = atof(arg[2]);
    cut_one = atof(arg[3]);
  }

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut_inner[i][j] = cut_inner_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLubricateU::init_style()
{
  if (!atom->sphere_flag)
    error->all(FLERR,"Pair lubricateU requires atom style sphere");
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair lubricateU requires ghost atoms store velocity");

  int irequest = neighbor->request(this);

  // require that atom radii are identical within each type
  // require monodisperse system with same radii for all types

  double rad,radtype;
  for (int i = 1; i <= atom->ntypes; i++) {
    if (!atom->radius_consistency(i,radtype))
      error->all(FLERR,"Pair lubricateU requires monodisperse particles");
    if (i > 1 && radtype != rad)
      error->all(FLERR,"Pair lubricateU requires monodisperse particles");
    rad = radtype;
  }
  
  // set the isotropic constants depending on the volume fraction
  // vol_T = total volume

  double vol_T = domain->xprd*domain->yprd*domain->zprd; 
  
  // assuming monodisperse spheres, vol_P = volume of the particles

  double tmp = 0.0;
  if (atom->radius) tmp = atom->radius[0];
  double radi;
  MPI_Allreduce(&tmp,&radi,1,MPI_DOUBLE,MPI_MAX,world);

  double vol_P = atom->natoms * (4.0/3.0)*MY_PI*pow(radi,3);
  
  // vol_f = volume fraction

  double vol_f = vol_P/vol_T;
  
  // set the isotropic constant

  if (flaglog == 0) {
    R0  = 6*MY_PI*mu*radi*(1.0 + 2.16*vol_f);
    RT0 = 8*MY_PI*mu*pow(radi,3);  // not actually needed
    RS0 = 20.0/3.0*MY_PI*mu*pow(radi,3)*(1.0 + 3.33*vol_f + 2.80*vol_f*vol_f);
  } else {
    R0  = 6*MY_PI*mu*radi*(1.0 + 2.725*vol_f - 6.583*vol_f*vol_f);
    RT0 = 8*MY_PI*mu*pow(radi,3)*(1.0 + 0.749*vol_f - 2.469*vol_f*vol_f);
    RS0 = 20.0/3.0*MY_PI*mu*pow(radi,3)*(1.0 + 3.64*vol_f - 6.95*vol_f*vol_f);
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLubricateU::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    cut_inner[i][j] = mix_distance(cut_inner[i][i],cut_inner[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  cut_inner[j][i] = cut_inner[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file 
------------------------------------------------------------------------- */

void PairLubricateU::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
	fwrite(&cut_inner[i][j],sizeof(double),1,fp);
	fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLubricateU::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
	if (me == 0) {
	  fread(&cut_inner[i][j],sizeof(double),1,fp);
	  fread(&cut[i][j],sizeof(double),1,fp);
	}
	MPI_Bcast(&cut_inner[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLubricateU::write_restart_settings(FILE *fp)
{
  fwrite(&mu,sizeof(double),1,fp);
  fwrite(&flaglog,sizeof(int),1,fp);
  fwrite(&cut_inner_global,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLubricateU::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&mu,sizeof(double),1,fp);
    fread(&flaglog,sizeof(int),1,fp);
    fread(&cut_inner_global,sizeof(double),1,fp);
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&mu,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&flaglog,1,MPI_INT,0,world);
  MPI_Bcast(&cut_inner_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/*---------------------------------------------------------------------------*/

void PairLubricateU::copy_vec_uo(int inum, double *xcg, 
				 double **v, double **omega)
{
  int i,j,ii;
  int *ilist;
  int itype;
  double radi;
  double inertia;

  double *rmass = atom->rmass;
  int *type = atom->type;
  
  ilist = list->ilist;
  
  for (ii=0;ii<inum;ii++) {
    i = ilist[ii];
    itype = type[i];
    radi = atom->radius[i];
    inertia = 0.4*rmass[i]*radi*radi;    
    
    for (j=0;j<3;j++) {       
      v[i][j] = xcg[6*ii+j];
      omega[i][j] = xcg[6*ii+j+3];
    }
  }
}

/*---------------------------------------------------------------------------*/

void PairLubricateU::copy_uo_vec(int inum, double **f, double **torque, 
				 double *RU)
{
  int i,j,ii;
  int *ilist;
  
  ilist = list->ilist;
  
  for (ii=0;ii<inum;ii++) {
    i = ilist[ii];
    for (j=0;j<3;j++) {       
      RU[6*ii+j] = f[i][j];
      RU[6*ii+j+3] = torque[i][j];
    }
  }
}

/* ---------------------------------------------------------------------- */

int PairLubricateU::pack_comm(int n, int *list, double *buf,
			      int pbc_flag, int *pbc)
{
  int i,j,m;

  double **v = atom->v;
  double **omega = atom->omega;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = v[j][0];
    buf[m++] = v[j][1];
    buf[m++] = v[j][2];
    buf[m++] = omega[j][0];
    buf[m++] = omega[j][1];
    buf[m++] = omega[j][2];
  }
  return 6;
}

/* ---------------------------------------------------------------------- */

void PairLubricateU::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  double **v = atom->v;
  double **omega = atom->omega;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    omega[i][0] = buf[m++];
    omega[i][1] = buf[m++];
    omega[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

double PairLubricateU::dot_vec_vec(int N, double *x, double *y)
{
  double dotp=0.0;
  for (int i = 0; i < N; i++) dotp += x[i]*y[i];
  return dotp;
}
