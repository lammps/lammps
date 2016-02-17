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
   Contributing authors: 
   James Larentzos (U.S. Army Research Laboratory)
   and Timothy I. Mattox (Engility Corporation)

   Martin Lisal (Institute of Chemical Process Fundamentals 
   of the Czech Academy of Sciences and J. E. Purkinje University)

   John Brennan, Joshua Moore and William Mattson (Army Research Lab)

   Please cite the related publications:
   J. P. Larentzos, J. K. Brennan, J. D. Moore, M. Lisal, W. D. Mattson,
   "Parallel implementation of isothermal and isoenergetic Dissipative
   Particle Dynamics using Shardlow-like splitting algorithms", 
   Computer Physics Communications, 2014, 185, pp 1987--1998.

   M. Lisal, J. K. Brennan, J. Bonet Avalos, "Dissipative particle dynamics
   at isothermal, isobaric, isoenergetic, and isoenthalpic conditions using
   Shardlow-like splitting algorithms", Journal of Chemical Physics, 2011,
   135, 204105.
------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix_shardlow.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include <math.h>
#include "atom_vec.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "random_mars.h"
#include "memory.h"
#include "domain.h"
#include "modify.h"
#include "pair_dpd_fdt.h"
#include "pair_dpd_fdt_energy.h"
#include "pair.h"
#include "citeme.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON 1.0e-10

static const char cite_fix_shardlow[] =
  "fix shardlow command:\n\n"
  "@Article{Larentzos14,\n"
  " author = {J. P. Larentzos, J. K. Brennan, J. D. Moore, M. Lisal, W. D. Mattson},\n"
  " title = {Parallel implementation of isothermal and isoenergetic Dissipative Particle Dynamics using Shardlow-like splitting algorithms},\n"
  " journal = {Computer Physics Communications},\n"
  " year =    2014,\n"
  " volume =  185\n"
  " pages =   {1987--1998}\n"
  "}\n\n"
  "@Article{Lisal11,\n"
  " author = {M. Lisal, J. K. Brennan, J. Bonet Avalos},\n"
  " title = {Dissipative particle dynamics at isothermal, isobaric, isoenergetic, and isoenthalpic conditions using Shardlow-like splitting algorithms},\n"
  " journal = {Journal of Chemical Physics},\n"
  " year =    2011,\n"
  " volume =  135\n"
  " pages =   {204105}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

FixShardlow::FixShardlow(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_shardlow);

  if (narg != 3) error->all(FLERR,"Illegal fix shardlow command");

  pairDPD = NULL;
  pairDPDE = NULL;
  pairDPD = (PairDPDfdt *) force->pair_match("dpd/fdt",1);
  pairDPDE = (PairDPDfdtEnergy *) force->pair_match("dpd/fdt/energy",1);

  if(pairDPDE){
    comm_forward = 10;
    comm_reverse = 5;
  } else {
    comm_forward = 6;
    comm_reverse = 3;
  }

  if(pairDPD == NULL && pairDPDE == NULL)
    error->all(FLERR,"Must use pair_style dpd/fdt or dpd/fdt/energy with fix shardlow");
  
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"nve") == 0 || strcmp(modify->fix[i]->style,"nph") == 0)
      error->all(FLERR,"A deterministic integrator must be specified after fix shardlow in input file (e.g. fix nve or fix nph).");
}

/* ---------------------------------------------------------------------- */

int FixShardlow::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= PRE_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixShardlow::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixShardlow::setup(int vflag)
{
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"nvt") == 0 || strcmp(modify->fix[i]->style,"npt") == 0)
      error->all(FLERR,"Cannot use constant temperature integration routines with DPD.");
}

/* ---------------------------------------------------------------------- */

void FixShardlow::setup_pre_force(int vflag)
{
  neighbor->build_one(list);
}

/* ----------------------------------------------------------------------
   Perform the stochastic integration and Shardlow update
   Allow for both per-type and per-atom mass

   NOTE: only implemented for orthogonal boxes, not triclinic
------------------------------------------------------------------------- */

void FixShardlow::initial_integrate(int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double delvx,delvy,delvz;
  double rsq,r,rinv;
  double dot,wd,wr,randnum,factor_dpd,factor_dpd1;
  double dpx,dpy,dpz;
  double denom, mu_ij;

  double **x = atom->x;
  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int nall = nlocal + nghost;
  int newton_pair = force->newton_pair;
  double randPair;

  double *uCond = atom->uCond;
  double *uMech = atom->uMech;
  double *duCond = atom->duCond;
  double *duMech = atom->duMech;
  double *dpdTheta = atom->dpdTheta;
  double kappa_ij, alpha_ij, theta_ij, gamma_ij, sigma_ij;
  double vxi, vyi, vzi, vxj, vyj, vzj;
  double vx0i, vy0i, vz0i, vx0j, vy0j, vz0j;
  double dot1, dot2, dot3, dot4;
  double mass_i, mass_j;
  double massinv_i, massinv_j;
  double cut, cut2;

  const double dt     = update->dt;
  const double dtsqrt = sqrt(dt);

  // NOTE: this logic is specific to orthogonal boxes, not triclinic

  // Enforce the constraint that ghosts must be contained in the nearest sub-domains
  double bbx = domain->subhi[0] - domain->sublo[0];
  double bby = domain->subhi[1] - domain->sublo[1];
  double bbz = domain->subhi[2] - domain->sublo[2];

  double rcut = double(2.0)*neighbor->cutneighmax;

  if (domain->triclinic)
    error->all(FLERR,"Fix shardlow does not yet support triclinic geometries");

  if(rcut >= bbx || rcut >= bby || rcut>= bbz )
    error->all(FLERR,"Shardlow algorithm requires sub-domain length > 2*(rcut+skin). Either reduce the number of processors requested, or change the cutoff/skin\n");

  // Allocate memory for the dvSSA arrays
  dvSSA = new double*[nall];
  for (ii = 0; ii < nall; ii++) {
    dvSSA[ii] = new double[3];
  }

  // Zero the momenta
  for (ii = 0; ii < nlocal; ii++) {
    dvSSA[ii][0] = double(0.0);
    dvSSA[ii][1] = double(0.0);
    dvSSA[ii][2] = double(0.0);
    if(pairDPDE){
      duCond[ii] = double(0.0);
      duMech[ii] = double(0.0);
    }
  }

  // Communicate the updated momenta and velocities to all nodes
  comm->forward_comm_fix(this);

  // Define pointers to access the neighbor list
  if(pairDPDE){
    inum = pairDPDE->list->inum;
    ilist = pairDPDE->list->ilist;
    numneigh = pairDPDE->list->numneigh;
    firstneigh = pairDPDE->list->firstneigh;
  } else {
    inum = pairDPD->list->inum;
    ilist = pairDPD->list->ilist;
    numneigh = pairDPD->list->numneigh;
    firstneigh = pairDPD->list->firstneigh;
  }

  //Loop over all 14 directions (8 stages)  
  for (int idir = 1; idir <=8; idir++){
    
    // Loop over neighbors of my atoms
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      
      itype = type[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];
      
      // Loop over Directional Neighbors only
      for (jj = 0; jj < jnum; jj++) {
	j = jlist[jj];
	j &= NEIGHMASK;
	if (neighbor->ssa_airnum[j] != idir) continue;
	jtype = type[j];
	
	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;

	if(pairDPDE){
	  cut2 = pairDPDE->cutsq[itype][jtype];
	  cut  = pairDPDE->cut[itype][jtype];
	} else {
	  cut2 = pairDPD->cutsq[itype][jtype];
	  cut  = pairDPD->cut[itype][jtype];
	}
	  
	// if (rsq < pairDPD->cutsq[itype][jtype]) {
	if (rsq < cut2) {
	  r = sqrt(rsq);
	  if (r < EPSILON) continue;     // r can be 0.0 in DPD systems
	  rinv = double(1.0)/r;
	  
	  // Store the velocities from previous Shardlow step
	  vx0i = v[i][0] + dvSSA[i][0];
	  vy0i = v[i][1] + dvSSA[i][1];
	  vz0i = v[i][2] + dvSSA[i][2];
	  
	  vx0j = v[j][0] + dvSSA[j][0];
	  vy0j = v[j][1] + dvSSA[j][1];
	  vz0j = v[j][2] + dvSSA[j][2];
	  
	  // Compute the velocity difference between atom i and atom j
	  delvx = vx0i - vx0j;
	  delvy = vy0i - vy0j;
	  delvz = vz0i - vz0j;
	  
	  dot = (delx*delvx + dely*delvy + delz*delvz);
	  // wr = double(1.0) - r/pairDPD->cut[itype][jtype];
	  wr = double(1.0) - r/cut;
	  wd = wr*wr;

	  if(pairDPDE){
	    // Compute the current temperature
	    theta_ij = double(0.5)*(double(1.0)/dpdTheta[i] + double(1.0)/dpdTheta[j]);
	    theta_ij = double(1.0)/theta_ij;
	    sigma_ij = pairDPDE->sigma[itype][jtype];
	    randnum = pairDPDE->random->gaussian();
	  } else {
	    theta_ij = pairDPD->temperature;
	    sigma_ij = pairDPD->sigma[itype][jtype];
	    randnum = pairDPD->random->gaussian();
	  }

	  gamma_ij = sigma_ij*sigma_ij / (2.0*force->boltz*theta_ij);
	  randPair = sigma_ij*wr*randnum*dtsqrt;
	  
	  factor_dpd = -dt*gamma_ij*wd*dot*rinv; 
	  factor_dpd += randPair;
	  factor_dpd *= double(0.5);
	  
	  // Compute momentum change between t and t+dt 
	  dpx  = factor_dpd*delx*rinv;
	  dpy  = factor_dpd*dely*rinv;
	  dpz  = factor_dpd*delz*rinv;
	  
	  if (rmass) {
	    mass_i = rmass[i];
	    mass_j = rmass[j];
	  } else {
	    mass_i = mass[itype];
	    mass_j = mass[jtype];
	  }
	  massinv_i = double(1.0) / mass_i;
	  massinv_j = double(1.0) / mass_j;
	  
	  // Update the delta velocity on i
	  dvSSA[i][0] += dpx*force->ftm2v*massinv_i;
	  dvSSA[i][1] += dpy*force->ftm2v*massinv_i;
	  dvSSA[i][2] += dpz*force->ftm2v*massinv_i;
	  
	  if (newton_pair || j < nlocal) {
	    // Update the delta velocity on j
	    dvSSA[j][0] -= dpx*force->ftm2v*massinv_j;
	    dvSSA[j][1] -= dpy*force->ftm2v*massinv_j;
	    dvSSA[j][2] -= dpz*force->ftm2v*massinv_j;
	  }
	  
	  //ii.   Compute the velocity diff
	  delvx = v[i][0] + dvSSA[i][0] - v[j][0] - dvSSA[j][0];
	  delvy = v[i][1] + dvSSA[i][1] - v[j][1] - dvSSA[j][1];
	  delvz = v[i][2] + dvSSA[i][2] - v[j][2] - dvSSA[j][2];
	  
	  dot = delx*delvx + dely*delvy + delz*delvz;
	  
	  //iii.    Compute dpi again
	  mu_ij = massinv_i + massinv_j;
	  denom = double(1.0) + double(0.5)*mu_ij*gamma_ij*wd*dt*force->ftm2v;
	  factor_dpd = -double(0.5)*dt*gamma_ij*wd*force->ftm2v/denom;
	  factor_dpd1 = factor_dpd*(mu_ij*randPair);
	  factor_dpd1 += randPair;
	  factor_dpd1 *= double(0.5);
	  
	  // Compute the momentum change between t and t+dt  
	  dpx  = (factor_dpd*dot*rinv/force->ftm2v + factor_dpd1)*delx*rinv;
	  dpy  = (factor_dpd*dot*rinv/force->ftm2v + factor_dpd1)*dely*rinv;
	  dpz  = (factor_dpd*dot*rinv/force->ftm2v + factor_dpd1)*delz*rinv;
	  
	  //Update the velocity change on i
	  dvSSA[i][0] += dpx*force->ftm2v*massinv_i;
	  dvSSA[i][1] += dpy*force->ftm2v*massinv_i;
	  dvSSA[i][2] += dpz*force->ftm2v*massinv_i;
	  
	  if (newton_pair || j < nlocal) {
	    //Update the velocity change on j
	    dvSSA[j][0] -= dpx*force->ftm2v*massinv_j;
	    dvSSA[j][1] -= dpy*force->ftm2v*massinv_j;
	    dvSSA[j][2] -= dpz*force->ftm2v*massinv_j;
	  }

	  if(pairDPDE){
	    // Compute uCond
	    randnum = pairDPDE->random->gaussian();
	    kappa_ij = pairDPDE->kappa[itype][jtype];
	    alpha_ij = sqrt(2.0*force->boltz*kappa_ij);
	    randPair = alpha_ij*wr*randnum*dtsqrt;
	    
	    factor_dpd = kappa_ij*(double(1.0)/dpdTheta[i] - double(1.0)/dpdTheta[j])*wd*dt;
	    factor_dpd += randPair;
	    
	    duCond[i] += factor_dpd;
	    if (newton_pair || j < nlocal) {
	      duCond[j] -= factor_dpd;
	    }
	    
	    // Compute uMech
	    vxi = v[i][0] + dvSSA[i][0];
	    vyi = v[i][1] + dvSSA[i][1];
	    vzi = v[i][2] + dvSSA[i][2];
	    
	    vxj = v[j][0] + dvSSA[j][0];
	    vyj = v[j][1] + dvSSA[j][1];
	    vzj = v[j][2] + dvSSA[j][2];
	    
	    dot1 = vxi*vxi + vyi*vyi + vzi*vzi;
	    dot2 = vxj*vxj + vyj*vyj + vzj*vzj;
	    dot3 = vx0i*vx0i + vy0i*vy0i + vz0i*vz0i;
	    dot4 = vx0j*vx0j + vy0j*vy0j + vz0j*vz0j;
	    
	    dot1 = dot1*mass_i;
	    dot2 = dot2*mass_j;
	    dot3 = dot3*mass_i;
	    dot4 = dot4*mass_j;
	    
	    factor_dpd = double(0.25)*(dot1+dot2-dot3-dot4)/force->ftm2v;
	    duMech[i] -= factor_dpd;
	    if (newton_pair || j < nlocal) {
	      duMech[j] -= factor_dpd;
	    }
	  }
	}
      }
    }
    
    // Communicate the ghost delta velocities to the locally owned atoms
    comm->reverse_comm_fix(this);
    
    // Zero the dv
    for (ii = 0; ii < nlocal; ii++) {
      // Shardlow update
      v[ii][0] += dvSSA[ii][0];
      v[ii][1] += dvSSA[ii][1];
      v[ii][2] += dvSSA[ii][2];
      dvSSA[ii][0] = double(0.0);
      dvSSA[ii][1] = double(0.0);
      dvSSA[ii][2] = double(0.0);
      if(pairDPDE){
	uCond[ii] += duCond[ii];
	uMech[ii] += duMech[ii];
	duCond[ii] = double(0.0);
	duMech[ii] = double(0.0);
      }
    }
    
    // Communicate the updated momenta and velocities to all nodes
    comm->forward_comm_fix(this);
    
  }  //End Loop over all directions For idir = Top, Top-Right, Right, Bottom-Right, Back
  
  for (ii = 0; ii < nall; ii++) {
    delete dvSSA[ii];
  }
  delete [] dvSSA;
}

/* ----------------------------------------------------------------------
 *    assign owned and ghost atoms their ssa active interaction region numbers
------------------------------------------------------------------------- */

void FixShardlow::setup_pre_neighbor()
{
  neighbor->assign_ssa_airnums();
}

void FixShardlow::pre_neighbor()
{
  neighbor->assign_ssa_airnums();
}

/* ---------------------------------------------------------------------- */

int FixShardlow::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int ii,jj,m;
  double **v  = atom->v;
  double *duCond = atom->duCond;
  double *duMech = atom->duMech;
  double *uCond = atom->uCond;
  double *uMech = atom->uMech;

  m = 0;
  for (ii = 0; ii < n; ii++) {
    jj = list[ii];
    buf[m++] = dvSSA[jj][0];
    buf[m++] = dvSSA[jj][1];
    buf[m++] = dvSSA[jj][2];
    buf[m++] = v[jj][0];
    buf[m++] = v[jj][1];
    buf[m++] = v[jj][2];
    if(pairDPDE){
      buf[m++] = duCond[jj];
      buf[m++] = duMech[jj];
      buf[m++] = uCond[jj];
      buf[m++] = uMech[jj];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixShardlow::unpack_forward_comm(int n, int first, double *buf)
{
  int ii,m,last;
  double **v  = atom->v;
  double *duCond = atom->duCond;
  double *duMech = atom->duMech;
  double *uCond = atom->uCond;
  double *uMech = atom->uMech;

  m = 0;
  last = first + n ;
  for (ii = first; ii < last; ii++) {
    dvSSA[ii][0] = buf[m++];
    dvSSA[ii][1] = buf[m++];
    dvSSA[ii][2] = buf[m++];
    v[ii][0] = buf[m++];
    v[ii][1] = buf[m++];
    v[ii][2] = buf[m++];
    if(pairDPDE){
      duCond[ii] = buf[m++];
      duMech[ii] = buf[m++];
      uCond[ii] = buf[m++];
      uMech[ii] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixShardlow::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  double *duCond = atom->duCond;
  double *duMech = atom->duMech;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = dvSSA[i][0];
    buf[m++] = dvSSA[i][1];
    buf[m++] = dvSSA[i][2];
    if(pairDPDE){
      buf[m++] = duCond[i];
      buf[m++] = duMech[i];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixShardlow::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  double *duCond = atom->duCond;
  double *duMech = atom->duMech;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];

    dvSSA[j][0] += buf[m++];
    dvSSA[j][1] += buf[m++];
    dvSSA[j][2] += buf[m++];
    if(pairDPDE){
      duCond[j] += buf[m++];
      duMech[j] += buf[m++];
    }
  }
}
