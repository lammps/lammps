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
                         Pieter in 't Veld (BASF), code restructuring
                         Dave Heine (Corning), polydispersity
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_lubricateU_poly.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "domain.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "fix_deform.h"
#include "fix_wall.h"
#include "input.h"
#include "variable.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define TOL 1E-3   // tolerance for conjugate gradient

// same as fix_wall.cpp

enum{EDGE,CONSTANT,VARIABLE};


/* ---------------------------------------------------------------------- */

PairLubricateUPoly::PairLubricateUPoly(LAMMPS *lmp) :
  PairLubricateU(lmp) {}

/* ----------------------------------------------------------------------
   It first has to solve for the velocity of the particles such that
   the net force on the particles is zero. NOTE: it has to be the last
   type of pair interaction specified in the input file. Also, it
   assumes that no other types of interactions, like k-space, is
   present. As already mentioned, the net force on the particles after
   this pair interaction would be identically zero.
   ---------------------------------------------------------------------- */

void PairLubricateUPoly::compute(int eflag, int vflag)
{
  int i,j;

  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int nall = nlocal + nghost;

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

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

  if (6*list->inum > cgmax) {
    memory->sfree(bcg);
    memory->sfree(xcg);
    memory->sfree(rcg);
    memory->sfree(rcg1);
    memory->sfree(pcg);
    memory->sfree(RU);
    cgmax = 6*list->inum;
    memory->create(bcg,cgmax,"pair:bcg");
    memory->create(xcg,cgmax,"pair:bcg");
    memory->create(rcg,cgmax,"pair:bcg");
    memory->create(rcg1,cgmax,"pair:bcg");
    memory->create(pcg,cgmax,"pair:bcg");
    memory->create(RU,cgmax,"pair:bcg");
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

  iterate(atom->x,1);

  // Find positions at half the timestep and store in xl

  intermediates(nall,xl);

  // Store back the saved forces and torques in original arrays

  for(i=0;i<nlocal+nghost;i++) {
    for(j=0;j<3;j++) {
      f[i][j] = fl[i][j];
      torque[i][j] = Tl[i][j];
    }
  }

  // Stage two: This will give the final velocities

  iterate(xl,2);
}

/* ------------------------------------------------------------------------
   Stage one of midpoint method
------------------------------------------------------------------------- */

void PairLubricateUPoly::iterate(double **x, int stage)
{
  int i,j,ii,itype;

  int inum = list->inum;
  int *ilist = list->ilist;
  int *type = atom->type;
  int newton_pair = force->newton_pair;

  double alpha,beta;
  double normi,error,normig;
  double send[2],recv[2],rcg_dot_rcg;
  double mo_inertia,radi;

  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *radius = atom->radius;

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

  // compute the viscosity/pressure

  if (evflag && stage == 2) compute_Fh(x);

  // find actual particle's velocities from relative velocities
  // only non-zero component of fluid's vel : vx=gdot*y and wz=-gdot/2

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

void PairLubricateUPoly::compute_Fh(double **x)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  int *ilist,*jlist,*numneigh,**firstneigh;

  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int newton_pair = force->newton_pair;
  int overlaps = 0;

  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz;
  double rsq,r,h_sep,radi,radj;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3;
  double vt1,vt2,vt3;
  double inv_inertia;
  double vi[3],vj[3],wi[3],wj[3],xl[3],jl[3],pre[2];

  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *radius = atom->radius;

  double beta[2][5];
  double vxmu2f = force->vxmu2f;
  double a_sq = 0.0;
  double a_sh = 0.0;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  beta[0][0] = beta[1][0] = beta[1][4] = 0.0;

  // This section of code adjusts R0/RT0/RS0 if necessary due to changes
  // in the volume fraction as a result of fix deform or moving walls

  double dims[3], wallcoord;
  if (flagVF) // Flag for volume fraction corrections
    if (flagdeform || flagwall == 2){ // Possible changes in volume fraction
      if (flagdeform && !flagwall)
        for (j = 0; j < 3; j++)
          dims[j] = domain->prd[j];
      else if (flagwall == 2 || (flagdeform && flagwall == 1)){
         double wallhi[3], walllo[3];
         for (int j = 0; j < 3; j++){
           wallhi[j] = domain->prd[j];
           walllo[j] = 0;
         }
         for (int m = 0; m < wallfix->nwall; m++){
           int dim = wallfix->wallwhich[m] / 2;
           int side = wallfix->wallwhich[m] % 2;
           if (wallfix->xstyle[m] == VARIABLE){
             wallcoord = input->variable->compute_equal(wallfix->xindex[m]);
           }
           else wallcoord = wallfix->coord0[m];
           if (side == 0) walllo[dim] = wallcoord;
           else wallhi[dim] = wallcoord;
         }
         for (int j = 0; j < 3; j++)
           dims[j] = wallhi[j] - walllo[j];
      }
      double vol_T = dims[0]*dims[1]*dims[2];
      double vol_f = vol_P/vol_T;
      if (flaglog == 0) {
        //R0  = 6*MY_PI*mu*(1.0 + 2.16*vol_f);
        //RT0 = 8*MY_PI*mu;
        RS0 = 20.0/3.0*MY_PI*mu*(1.0 + 3.33*vol_f + 2.80*vol_f*vol_f);
      } else {
        //R0  = 6*MY_PI*mu*(1.0 + 2.725*vol_f - 6.583*vol_f*vol_f);
        //RT0 = 8*MY_PI*mu*(1.0 + 0.749*vol_f - 2.469*vol_f*vol_f);
        RS0 = 20.0/3.0*MY_PI*mu*(1.0 + 3.64*vol_f - 6.95*vol_f*vol_f);
      }
    }

  // end of R0 adjustment code
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
    pre[1] = 8.0*(pre[0] = MY_PI*mu*radi)*radi*radi; // BROKEN?? Should be "+"??
    pre[0] *= 6.0;

    // Find the contribution to stress from isotropic RS0
    // Set psuedo force to obtain the required contribution
    // need to set delx  and fy only

    fx = 0.0; delx = radi;
    fy = vxmu2f*RS0*pow(radi,3.0)*gdot/2.0/radi; dely = 0.0;
    fz = 0.0; delz = 0.0;
    if (evflag)
      ev_tally_xyz(i,i,nlocal,newton_pair,0.0,0.0,-fx,-fy,-fz,delx,dely,delz);

    // Find angular velocity

    wi[0] = omega[i][0];
    wi[1] = omega[i][1];
    wi[2] = omega[i][2];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      radj = radius[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);

        // Use omega directly if it exists, else angmom
        // angular momentum = I*omega = 2/5 * M*R^2 * omega

        wj[0] = omega[j][0];
        wj[1] = omega[j][1];
        wj[2] = omega[j][2];

        // loc of the point of closest approach on particle i from its cente
        // POC for j is in opposite direction as for i

        xl[0] = -delx/r*radi;
        xl[1] = -dely/r*radi;
        xl[2] = -delz/r*radi;
        jl[0] = delx/r*radj;
        jl[1] = dely/r*radj;
        jl[2] = delz/r*radj;

        h_sep = r - radi-radj;

        // velocity at the point of closest approach on both particles
        // v = v + omega_cross_xl

        // particle i

        vi[0] = v[i][0] + (wi[1]*xl[2] - wi[2]*xl[1]);
        vi[1] = v[i][1] + (wi[2]*xl[0] - wi[0]*xl[2]);
        vi[2] = v[i][2] + (wi[0]*xl[1] - wi[1]*xl[0]);

        // particle j

        vj[0] = v[j][0] + (wj[1]*jl[2] - wj[2]*jl[1]);
        vj[1] = v[j][1] + (wj[2]*jl[0] - wj[0]*jl[2]);
        vj[2] = v[j][2] + (wj[0]*jl[1] - wj[1]*jl[0]);


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

        // check for overlaps

        if (h_sep < 0.0) overlaps++;

        // If less than the minimum gap use the minimum gap instead

        if (r < cut_inner[itype][jtype])
          h_sep = cut_inner[itype][jtype] - radi-radj;

        // Scale h_sep by radi

        h_sep = h_sep/radi;
        beta[0][1] = radj/radi;
        beta[1][1] = 1.0 + beta[0][1];

        /*beta0 = radj/radi;
        beta1 = 1.0 + beta0;*/

        // Scalar resistances

        if (flaglog) {
          beta[0][2] = beta[0][1]*beta[0][1];
          beta[0][3] = beta[0][2]*beta[0][1];
          beta[0][4] = beta[0][3]*beta[0][1];
          beta[1][2] = beta[1][1]*beta[1][1];
          beta[1][3] = beta[1][2]*beta[1][1];
          double log_h_sep_beta13 = log(1.0/h_sep)/beta[1][3];
          double h_sep_beta11 = h_sep/beta[1][1];

          a_sq = pre[0]*(beta[0][2]/beta[1][2]/h_sep
                +((0.2+1.4*beta[0][1]+0.2*beta[0][2])
                  +(1.0+18.0*(beta[0][1]+beta[0][3])-29.0*beta[0][2]
                    +beta[0][4])*h_sep_beta11/21.0)*log_h_sep_beta13);

          a_sh = pre[0]*((8.0*(beta[0][1]+beta[0][3])+4.0*beta[0][2])/15.0
                +(64.0-180.0*(beta[0][1]+beta[0][3])+232.0*beta[0][2]
                  +64.0*beta[0][4])*h_sep_beta11/375.0)*log_h_sep_beta13;

          /*a_sq = beta0*beta0/beta1/beta1/h_sep
                  +(1.0+7.0*beta0+beta0*beta0)/5.0/pow(beta1,3)*log(1.0/h_sep);
          a_sq += (1.0+18.0*beta0-29.0*beta0*beta0+18.0*pow(beta0,3)
                  +pow(beta0,4))/21.0/pow(beta1,4)*h_sep*log(1.0/h_sep);
          a_sq *= 6.0*MY_PI*mu*radi;

          a_sh = 4.0*beta0*(2.0+beta0
                  +2.0*beta0*beta0)/15.0/pow(beta1,3)*log(1.0/h_sep);
          a_sh += 4.0*(16.0-45.0*beta0+58.0*beta0*beta0-45.0*pow(beta0,3)
                  +16.0*pow(beta0,4))/375.0/pow(beta1,4)*h_sep*log(1.0/h_sep);
          a_sh *= 6.0*MY_PI*mu*radi;*/
        } else {
          //a_sq = 6.0*MY_PI*mu*radi*(beta0*beta0/beta1/beta1/h_sep);
          a_sq = pre[0]*(beta[0][1]*beta[0][1]/(beta[1][1]*beta[1][1]*h_sep));
        }

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

        // set j = nlocal so that only I gets tallied

        if (evflag) ev_tally_xyz(i,nlocal,nlocal,0,
                                 0.0,0.0,-fx,-fy,-fz,delx,dely,delz);
      }
    }
  }
}

/* ----------------------------------------------------------------------
  computes R_FU * U
---------------------------------------------------------------------- */

void PairLubricateUPoly::compute_RU(double **x)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  int *ilist,*jlist,*numneigh,**firstneigh;

  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int overlaps = 0;

  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz,tx,ty,tz;
  double rsq,r,radi,radj,h_sep;
  //double beta0,beta1;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3;
  double vt1,vt2,vt3,wdotn,wt1,wt2,wt3;
  double inv_inertia;
  double vi[3],vj[3],wi[3],wj[3],xl[3],jl[3],pre[2];

  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *radius = atom->radius;

  double beta[2][5];
  double vxmu2f = force->vxmu2f;
  double a_sq = 0.0;
  double a_sh = 0.0;
  double a_pu = 0.0;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  beta[0][0] = beta[1][0] = beta[1][4] = 0.0;

 // This section of code adjusts R0/RT0/RS0 if necessary due to changes
  // in the volume fraction as a result of fix deform or moving walls

  double dims[3], wallcoord;
  if (flagVF) // Flag for volume fraction corrections
    if (flagdeform || flagwall == 2){ // Possible changes in volume fraction
      if (flagdeform && !flagwall)
        for (j = 0; j < 3; j++)
          dims[j] = domain->prd[j];
      else if (flagwall == 2 || (flagdeform && flagwall == 1)){
         double wallhi[3], walllo[3];
         for (j = 0; j < 3; j++){
           wallhi[j] = domain->prd[j];
           walllo[j] = 0;
         }
         for (int m = 0; m < wallfix->nwall; m++){
           int dim = wallfix->wallwhich[m] / 2;
           int side = wallfix->wallwhich[m] % 2;
           if (wallfix->xstyle[m] == VARIABLE){
             wallcoord = input->variable->compute_equal(wallfix->xindex[m]);
           }
           else wallcoord = wallfix->coord0[m];
           if (side == 0) walllo[dim] = wallcoord;
           else wallhi[dim] = wallcoord;
         }
         for (j = 0; j < 3; j++)
           dims[j] = wallhi[j] - walllo[j];
      }
      double vol_T = dims[0]*dims[1]*dims[2];
      double vol_f = vol_P/vol_T;
      if (flaglog == 0) {
        R0  = 6*MY_PI*mu*(1.0 + 2.16*vol_f);
        RT0 = 8*MY_PI*mu;
        //        RS0 = 20.0/3.0*MY_PI*mu*(1.0 + 3.33*vol_f + 2.80*vol_f*vol_f);
      } else {
        R0  = 6*MY_PI*mu*(1.0 + 2.725*vol_f - 6.583*vol_f*vol_f);
        RT0 = 8*MY_PI*mu*(1.0 + 0.749*vol_f - 2.469*vol_f*vol_f);
        //        RS0 = 20.0/3.0*MY_PI*mu*(1.0 + 3.64*vol_f - 6.95*vol_f*vol_f);
      }
    }

  // end of R0 adjustment code

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
    pre[1] = 8.0*(pre[0] = MY_PI*mu*radi)*radi*radi;
    pre[0] *= 6.0;

    // Find angular velocity

    wi[0] = omega[i][0];
    wi[1] = omega[i][1];
    wi[2] = omega[i][2];

    // Contribution due to the isotropic terms

    f[i][0] += -vxmu2f*R0*radi*v[i][0];
    f[i][1] += -vxmu2f*R0*radi*v[i][1];
    f[i][2] += -vxmu2f*R0*radi*v[i][2];

    const double radi3 = radi*radi*radi;
    torque[i][0] += -vxmu2f*RT0*radi3*wi[0];
    torque[i][1] += -vxmu2f*RT0*radi3*wi[1];
    torque[i][2] += -vxmu2f*RT0*radi3*wi[2];

    if (!flagHI) continue;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      radj = radius[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);

        wj[0] = omega[j][0];
        wj[1] = omega[j][1];
        wj[2] = omega[j][2];

        // loc of the point of closest approach on particle i from its center

        xl[0] = -delx/r*radi;
        xl[1] = -dely/r*radi;
        xl[2] = -delz/r*radi;
        jl[0] = delx/r*radj;
        jl[1] = dely/r*radj;
        jl[2] = delz/r*radj;

        // velocity at the point of closest approach on both particles
        // v = v + omega_cross_xl

        // particle i

        vi[0] = v[i][0] + (wi[1]*xl[2] - wi[2]*xl[1]);
        vi[1] = v[i][1] + (wi[2]*xl[0] - wi[0]*xl[2]);
        vi[2] = v[i][2] + (wi[0]*xl[1] - wi[1]*xl[0]);

        // particle j

        vj[0] = v[j][0] + (wj[1]*jl[2] - wj[2]*jl[1]);
        vj[1] = v[j][1] + (wj[2]*jl[0] - wj[0]*jl[2]);
        vj[2] = v[j][2] + (wj[0]*jl[1] - wj[1]*jl[0]);

        // Find the scalar resistances a_sq and a_sh

        h_sep = r - radi-radj;

        // check for overlaps

        if(h_sep < 0.0) overlaps++;

        // If less than the minimum gap use the minimum gap instead

        if (r < cut_inner[itype][jtype])
          h_sep = cut_inner[itype][jtype] - radi-radj;

        // Scale h_sep by radi

        h_sep = h_sep/radi;
        beta[0][1] = radj/radi;
        beta[1][1] = 1.0 + beta[0][1];

        // Scalar resistances

        if (flaglog) {
          beta[0][2] = beta[0][1]*beta[0][1];
          beta[0][3] = beta[0][2]*beta[0][1];
          beta[0][4] = beta[0][3]*beta[0][1];
          beta[1][2] = beta[1][1]*beta[1][1];
          beta[1][3] = beta[1][2]*beta[1][1];
          double log_h_sep_beta13 = log(1.0/h_sep)/beta[1][3];
          double h_sep_beta11 = h_sep/beta[1][1];

          a_sq = pre[0]*(beta[0][2]/beta[1][2]/h_sep
                +((0.2+1.4*beta[0][1]+0.2*beta[0][2])
                  +(1.0+18.0*(beta[0][1]+beta[0][3])-29.0*beta[0][2]
                    +beta[0][4])*h_sep_beta11/21.0)*log_h_sep_beta13);

          a_sh = pre[0]*((8.0*(beta[0][1]+beta[0][3])+4.0*beta[0][2])/15.0
                +(64.0-180.0*(beta[0][1]+beta[0][3])+232.0*beta[0][2]
                  +64.0*beta[0][4])*h_sep_beta11/375.0)*log_h_sep_beta13;

          a_pu = pre[1]*((0.4*beta[0][1]+0.1*beta[0][2])*beta[1][1]
                +(0.128-0.132*beta[0][1]+0.332*beta[0][2]
                  +0.172*beta[0][3])*h_sep)*log_h_sep_beta13;

          /*//a_sq = 6*MY_PI*mu*radi*(1.0/4.0/h_sep + 9.0/40.0*log(1/h_sep));
          a_sq = beta0*beta0/beta1/beta1/h_sep
                  +(1.0+7.0*beta0+beta0*beta0)/5.0/pow(beta1,3)*log(1.0/h_sep);
          a_sq += (1.0+18.0*beta0-29.0*beta0*beta0+18.0*pow(beta0,3)
                  +pow(beta0,4))/21.0/pow(beta1,4)*h_sep*log(1.0/h_sep);
          a_sq *= 6.0*MY_PI*mu*radi;

          a_sh = 4.0*beta0*(2.0+beta0
                  +2.0*beta0*beta0)/15.0/pow(beta1,3)*log(1.0/h_sep);
          a_sh += 4.0*(16.0-45.0*beta0+58.0*beta0*beta0-45.0*pow(beta0,3)
                  +16.0*pow(beta0,4))/375.0/pow(beta1,4)*h_sep*log(1.0/h_sep);
          a_sh *= 6.0*MY_PI*mu*radi;

          a_pu = beta0*(4.0+beta0)/10.0/beta1/beta1*log(1.0/h_sep);
          a_pu += (32.0-33.0*beta0+83.0*beta0*beta0
                  +43.0*pow(beta0,3))/250.0/pow(beta1,3)*h_sep*log(1.0/h_sep);
          a_pu *= 8.0*MY_PI*mu*pow(radi,3);*/
        } else
          a_sq = pre[0]*(beta[0][1]*beta[0][1]/(beta[1][1]*beta[1][1]*h_sep));

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

        // Find torque due to this force

        if (flaglog) {
          tx = xl[1]*fz - xl[2]*fy;
          ty = xl[2]*fx - xl[0]*fz;
          tz = xl[0]*fy - xl[1]*fx;

          // Why a scale factor ?

          torque[i][0] -= vxmu2f*tx;
          torque[i][1] -= vxmu2f*ty;
          torque[i][2] -= vxmu2f*tz;

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

void PairLubricateUPoly::compute_RE(double **x)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  int *ilist,*jlist,*numneigh,**firstneigh;

  int *type = atom->type;
  int overlaps = 0;

  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz,tx,ty,tz;
  double rsq,r,h_sep,radi,radj;
  //double beta0,beta1,lhsep;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3;
  double vt1,vt2,vt3;
  double xl[3],pre[2];

  double **f = atom->f;
  double **torque = atom->torque;
  double *radius = atom->radius;

  double beta[2][5];
  double vxmu2f = force->vxmu2f;
  double a_sq = 0.0;
  double a_sh = 0.0;
  double a_pu = 0.0;

  if (!flagHI) return;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  beta[0][0] = beta[1][0] = beta[1][4] = 0.0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    pre[1] = 8.0*(pre[0] = MY_PI*mu*radi)*radi*radi;
    pre[0] *= 6.0;

    // No contribution from isotropic terms due to E
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      radj = radius[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);

        // loc of the point of closest approach on particle i from its center

        xl[0] = -delx/r*radi;
        xl[1] = -dely/r*radi;
        xl[2] = -delz/r*radi;

        // Find the scalar resistances a_sq and a_sh

        h_sep = r - radi-radj;

        // check for overlaps

        if (h_sep < 0.0) overlaps++;

        // If less than the minimum gap use the minimum gap instead

        if (r < cut_inner[itype][jtype])
          h_sep = cut_inner[itype][jtype] - radi-radj;

        // Scale h_sep by radi

        h_sep = h_sep/radi;
        beta[0][1] = radj/radi;
        beta[1][1] = 1.0 + beta[0][1];

        /*beta0 = radj/radi;
        beta1 = 1.0 + beta0;
        lhsep = log(1.0/h_sep);*/

        // Scalar resistance for Squeeze type motions


        if (flaglog) {
          beta[0][2] = beta[0][1]*beta[0][1];
          beta[0][3] = beta[0][2]*beta[0][1];
          beta[0][4] = beta[0][3]*beta[0][1];
          beta[1][2] = beta[1][1]*beta[1][1];
          beta[1][3] = beta[1][2]*beta[1][1];
          double log_h_sep_beta13 = log(1.0/h_sep)/beta[1][3];
          double h_sep_beta11 = h_sep/beta[1][1];

          a_sq = pre[0]*(beta[0][2]/beta[1][2]/h_sep
                +((0.2+1.4*beta[0][1]+0.2*beta[0][2])
                  +(1.0+18.0*(beta[0][1]+beta[0][3])-29.0*beta[0][2]
                    +beta[0][4])*h_sep_beta11/21.0)*log_h_sep_beta13);

          a_sh = pre[0]*((8.0*(beta[0][1]+beta[0][3])+4.0*beta[0][2])/15.0
                +(64.0-180.0*(beta[0][1]+beta[0][3])+232.0*beta[0][2]
                  +64.0*beta[0][4])*h_sep_beta11/375.0)*log_h_sep_beta13;

          a_pu = pre[1]*((0.4*beta[0][1]+0.1*beta[0][2])*beta[1][1]
                +(0.128-0.132*beta[0][1]+0.332*beta[0][2]
                  +0.172*beta[0][3])*h_sep)*log_h_sep_beta13;

          /*//a_sq = 6*MY_PI*mu*radi*(1.0/4.0/h_sep + 9.0/40.0*log(1/h_sep));
          a_sq = beta0*beta0/beta1/beta1/h_sep
                  +(1.0+7.0*beta0+beta0*beta0)/5.0/pow(beta1,3)*lhsep;
          a_sq += (1.0+18.0*beta0-29.0*beta0*beta0+18.0*pow(beta0,3)
                  +pow(beta0,4))/21.0/pow(beta1,4)*h_sep*lhsep;
          a_sq *= 6.0*MY_PI*mu*radi;

          a_sh = 4.0*beta0*(2.0+beta0
                  +2.0*beta0*beta0)/15.0/pow(beta1,3)*log(1.0/h_sep);
          a_sh += 4.0*(16.0-45.0*beta0+58.0*beta0*beta0-45.0*pow(beta0,3)
                  +16.0*pow(beta0,4))/375.0/pow(beta1,4)*h_sep*log(1.0/h_sep);
          a_sh *= 6.0*MY_PI*mu*radi;

          a_pu = beta0*(4.0+beta0)/10.0/beta1/beta1*log(1.0/h_sep);
          a_pu += (32.0-33.0*beta0+83.0*beta0*beta0
                  +43.0*pow(beta0,3))/250.0/pow(beta1,3)*h_sep*log(1.0/h_sep);
          a_pu *= 8.0*MY_PI*mu*pow(radi,3);*/
        } else
          a_sq = pre[0]*(beta[0][1]*beta[0][1]/(beta[1][1]*beta[1][1]*h_sep));

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

        // Find torque due to this force

        if (flaglog) {
          tx = xl[1]*fz - xl[2]*fy;
          ty = xl[2]*fx - xl[0]*fz;
          tz = xl[0]*fy - xl[1]*fx;

          // Why a scale factor ?

          torque[i][0] -= vxmu2f*tx;
          torque[i][1] -= vxmu2f*ty;
          torque[i][2] -= vxmu2f*tz;

          // NOTE No a_pu term needed as they add up to zero
        }
      }
    }
  }

  int print_overlaps = 0;
  if (print_overlaps && overlaps)
    printf("Number of overlaps=%d\n",overlaps);
}

/*-----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLubricateUPoly::settings(int narg, char **arg)
{
  int itype;

  if (narg < 5 || narg > 7) error->all(FLERR,"Illegal pair_style command");

  mu = atof(arg[0]);
  flaglog = atoi(arg[1]);
  cut_inner_global = atof(arg[2]);
  cut_global = atof(arg[3]);
  gdot =  atof(arg[4]);

  flagHI = flagVF = 1;
  if (narg >= 6) flagHI = atoi(arg[5]);
  if (narg == 7) flagVF = atoi(arg[6]);

  if (flaglog == 1 && flagHI == 0) {
    error->warning(FLERR,"Cannot include log terms without 1/r terms; "
                   "setting flagHI to 1");
    flagHI = 1;
  }

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

  // Store the rate of strain tensor

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

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLubricateUPoly::init_style()
{
  if (force->newton_pair == 1)
    error->all(FLERR,"Pair lubricateU/poly requires newton pair off");
  if (comm->ghost_velocity == 0)
    error->all(FLERR,
               "Pair lubricateU/poly requires ghost atoms store velocity");
  if (!atom->sphere_flag)
    error->all(FLERR,"Pair lubricate/poly requires atom style sphere");

  // insure all particles are finite-size
  // for pair hybrid, should limit test to types using the pair style

  double *radius = atom->radius;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (radius[i] == 0.0)
      error->one(FLERR,"Pair lubricate/poly requires extended particles");

  // Set the isotropic constants depending on the volume fraction

  // Find the total volume
  // check for fix deform, if exists it must use "remap v"
  // If box will change volume, set appropriate flag so that volume
  // and v.f. corrections are re-calculated at every step.
  //
  // If available volume is different from box volume
  // due to walls, set volume appropriately; if walls will
  // move, set appropriate flag so that volume and v.f. corrections
  // are re-calculated at every step.

  flagdeform = flagwall = 0;
  for (int i = 0; i < modify->nfix; i++){
    if (strcmp(modify->fix[i]->style,"deform") == 0)
      flagdeform = 1;
    else if (strstr(modify->fix[i]->style,"wall") != NULL){
      if (flagwall) 
        error->all(FLERR,
                   "Cannot use multiple fix wall commands with "
                   "pair lubricateU");
      flagwall = 1; // Walls exist
      wallfix = (FixWall *) modify->fix[i];
      if (wallfix->xflag) flagwall = 2; // Moving walls exist
    }
  }

  // set the isotropic constants depending on the volume fraction
  // vol_T = total volumeshearing = flagdeform = flagwall = 0;

  double vol_T, wallcoord;
    if (!flagwall) vol_T = domain->xprd*domain->yprd*domain->zprd;
  else {
    double wallhi[3], walllo[3];
    for (int j = 0; j < 3; j++){
      wallhi[j] = domain->prd[j];
      walllo[j] = 0;
    }
    for (int m = 0; m < wallfix->nwall; m++){
      int dim = wallfix->wallwhich[m] / 2;
      int side = wallfix->wallwhich[m] % 2;
      if (wallfix->xstyle[m] == VARIABLE){
        wallfix->xindex[m] = input->variable->find(wallfix->xstr[m]);
        //Since fix->wall->init happens after pair->init_style
        wallcoord = input->variable->compute_equal(wallfix->xindex[m]);
      }

      else wallcoord = wallfix->coord0[m];

      if (side == 0) walllo[dim] = wallcoord;
      else wallhi[dim] = wallcoord;
    }
    vol_T = (wallhi[0] - walllo[0]) * (wallhi[1] - walllo[1]) *
      (wallhi[2] - walllo[2]);
  }

  // Assuming monodisperse spheres, find the volume of the particles

  double volP = 0.0;
  for (int i = 0; i < nlocal; i++)
    volP += (4.0/3.0)*MY_PI*pow(atom->radius[i],3.0);
  MPI_Allreduce(&volP,&vol_P,1,MPI_DOUBLE,MPI_SUM,world);

  double vol_f = vol_P/vol_T;

  //DRH volume fraction needs to be defined manually
  // if excluded volume regions are present
  //  vol_f=0.5;

  if (!flagVF) vol_f = 0;

  if (!comm->me) {
    if(logfile)
      fprintf(logfile, "lubricateU: vol_f = %g, vol_p = %g, vol_T = %g\n",
          vol_f,vol_P,vol_T);
    if (screen)
      fprintf(screen, "lubricateU: vol_f = %g, vol_p = %g, vol_T = %g\n",
          vol_f,vol_P,vol_T);
  }

  // Set the isotropic constant

  if (flaglog == 0) {
    R0  = 6*MY_PI*mu*(1.0 + 2.16*vol_f);
    RT0 = 8*MY_PI*mu;  // Not needed actually
    RS0 = 20.0/3.0*MY_PI*mu*(1.0 + 3.33*vol_f + 2.80*vol_f*vol_f);
  } else {
    R0  = 6*MY_PI*mu*(1.0 + 2.725*vol_f - 6.583*vol_f*vol_f);
    RT0 = 8*MY_PI*mu*(1.0 + 0.749*vol_f - 2.469*vol_f*vol_f);
    RS0 = 20.0/3.0*MY_PI*mu*(1.0 + 3.64*vol_f - 6.95*vol_f*vol_f);
  }

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}
