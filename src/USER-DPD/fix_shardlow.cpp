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
    comm_forward = 3;
    comm_reverse = 5;
  } else {
    comm_forward = 3;
    comm_reverse = 3;
  }

  if(pairDPD == NULL && pairDPDE == NULL)
    error->all(FLERR,"Must use pair_style dpd/fdt or dpd/fdt/energy with fix shardlow");

  // Setup the ssaAIR array
  atom->ssaAIR = NULL;
  grow_arrays(atom->nmax);
  for (int i = 0; i < atom->nlocal; i++) {
    atom->ssaAIR[i] = 1; /* coord2ssaAIR(x[i]) */
  }

  // Setup callbacks for maintaining atom->ssaAIR[]
  atom->add_callback(0); // grow (aka exchange)
  atom->add_callback(1); // restart
  atom->add_callback(2); // border
}

/* ---------------------------------------------------------------------- */

FixShardlow::~FixShardlow()
{
  atom->delete_callback(id, 0);
  atom->delete_callback(id, 1);
  atom->delete_callback(id, 2);

  memory->destroy(atom->ssaAIR);
}

/* ---------------------------------------------------------------------- */

int FixShardlow::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
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
  bool fixShardlow = false;

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"nvt") == 0 || strcmp(modify->fix[i]->style,"npt") == 0)
      error->all(FLERR,"Cannot use constant temperature integration routines with DPD.");

  for (int i = 0; i < modify->nfix; i++){
    if (strcmp(modify->fix[i]->style,"shardlow") == 0) fixShardlow = true;
    if (strcmp(modify->fix[i]->style,"nve") == 0 || (strcmp(modify->fix[i]->style,"nph") == 0)){
      if(fixShardlow) break;
      else error->all(FLERR,"The deterministic integrator must follow fix shardlow in the input file.");
    }
    if (i == modify->nfix-1) error->all(FLERR,"A deterministic integrator (e.g. fix nve or fix nph) is required when using fix shardlow.");
  }
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
  int newton_pair = force->newton_pair;
  double randPair;

  int *ssaAIR = atom->ssaAIR;
  double *uCond = atom->uCond;
  double *uMech = atom->uMech;
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

  double rcut = 2.0*neighbor->cutneighmax;

  if (domain->triclinic)
    error->all(FLERR,"Fix shardlow does not yet support triclinic geometries");

  if(rcut >= bbx || rcut >= bby || rcut>= bbz )
    error->all(FLERR,"Shardlow algorithm requires sub-domain length > 2*(rcut+skin). Either reduce the number of processors requested, or change the cutoff/skin\n");

  // Allocate memory for v_t0 to hold the initial velocities for the ghosts
  v_t0 = (double (*)[3]) memory->smalloc(sizeof(double)*3*nghost, "FixShardlow:v_t0");

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

    if (idir > 1) {
      // Communicate the updated velocities to all nodes
      comm->forward_comm_fix(this);

      if(pairDPDE){
        // Zero out the ghosts' uCond & uMech to be used as delta accumulators
        memset(&(uCond[nlocal]), 0, sizeof(double)*nghost);
        memset(&(uMech[nlocal]), 0, sizeof(double)*nghost);
      }
    }

    // Loop over neighbors of my atoms
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];

      // load velocity for i from memory
      vxi = v[i][0];
      vyi = v[i][1];
      vzi = v[i][2];

      itype = type[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      // Loop over Directional Neighbors only
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        if (ssaAIR[j] != idir) continue;
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

        // if (rsq < pairDPD->cutsq[itype][jtype])
        if (rsq < cut2) {
          r = sqrt(rsq);
          if (r < EPSILON) continue;     // r can be 0.0 in DPD systems
          rinv = 1.0/r;

          // Keep a copy of the velocities from previous Shardlow step
          vx0i = vxi;
          vy0i = vyi;
          vz0i = vzi;

          vx0j = vxj = v[j][0];
          vy0j = vyj = v[j][1];
          vz0j = vzj = v[j][2];

          // Compute the velocity difference between atom i and atom j
          delvx = vx0i - vx0j;
          delvy = vy0i - vy0j;
          delvz = vz0i - vz0j;

          dot = (delx*delvx + dely*delvy + delz*delvz);
          // wr = 1.0 - r/pairDPD->cut[itype][jtype];
          wr = 1.0 - r/cut;
          wd = wr*wr;

          if(pairDPDE){
            // Compute the current temperature
            theta_ij = 0.5*(1.0/dpdTheta[i] + 1.0/dpdTheta[j]);
            theta_ij = 1.0/theta_ij;
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
          factor_dpd *= 0.5;

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
          massinv_i = 1.0 / mass_i;
          massinv_j = 1.0 / mass_j;

          // Update the velocity on i
          vxi += dpx*force->ftm2v*massinv_i;
          vyi += dpy*force->ftm2v*massinv_i;
          vzi += dpz*force->ftm2v*massinv_i;

          if (newton_pair || j < nlocal) {
            // Update the velocity on j
            vxj -= dpx*force->ftm2v*massinv_j;
            vyj -= dpy*force->ftm2v*massinv_j;
            vzj -= dpz*force->ftm2v*massinv_j;
          }

          //ii.   Compute the velocity diff
          delvx = vxi - vxj;
          delvy = vyi - vyj;
          delvz = vzi - vzj;

          dot = delx*delvx + dely*delvy + delz*delvz;

          //iii.    Compute dpi again
          mu_ij = massinv_i + massinv_j;
          denom = 1.0 + 0.5*mu_ij*gamma_ij*wd*dt*force->ftm2v;
          factor_dpd = -0.5*dt*gamma_ij*wd*force->ftm2v/denom;
          factor_dpd1 = factor_dpd*(mu_ij*randPair);
          factor_dpd1 += randPair;
          factor_dpd1 *= 0.5;

          // Compute the momentum change between t and t+dt
          dpx  = (factor_dpd*dot*rinv/force->ftm2v + factor_dpd1)*delx*rinv;
          dpy  = (factor_dpd*dot*rinv/force->ftm2v + factor_dpd1)*dely*rinv;
          dpz  = (factor_dpd*dot*rinv/force->ftm2v + factor_dpd1)*delz*rinv;

          // Update the velocity on i
          vxi += dpx*force->ftm2v*massinv_i;
          vyi += dpy*force->ftm2v*massinv_i;
          vzi += dpz*force->ftm2v*massinv_i;

          if (newton_pair || j < nlocal) {
            // Update the velocity on j
            vxj -= dpx*force->ftm2v*massinv_j;
            vyj -= dpy*force->ftm2v*massinv_j;
            vzj -= dpz*force->ftm2v*massinv_j;
            // Store updated velocity for j
            v[j][0] = vxj;
            v[j][1] = vyj;
            v[j][2] = vzj;
          }

          if(pairDPDE){
            // Compute uCond
            randnum = pairDPDE->random->gaussian();
            kappa_ij = pairDPDE->kappa[itype][jtype];
            alpha_ij = sqrt(2.0*force->boltz*kappa_ij);
            randPair = alpha_ij*wr*randnum*dtsqrt;

            factor_dpd = kappa_ij*(1.0/dpdTheta[i] - 1.0/dpdTheta[j])*wd*dt;
            factor_dpd += randPair;

            uCond[i] += factor_dpd;
            if (newton_pair || j < nlocal) {
              uCond[j] -= factor_dpd;
            }

            // Compute uMech
            dot1 = vxi*vxi + vyi*vyi + vzi*vzi;
            dot2 = vxj*vxj + vyj*vyj + vzj*vzj;
            dot3 = vx0i*vx0i + vy0i*vy0i + vz0i*vz0i;
            dot4 = vx0j*vx0j + vy0j*vy0j + vz0j*vz0j;

            dot1 = dot1*mass_i;
            dot2 = dot2*mass_j;
            dot3 = dot3*mass_i;
            dot4 = dot4*mass_j;

            factor_dpd = 0.25*(dot1+dot2-dot3-dot4)/force->ftm2v;
            uMech[i] -= factor_dpd;
            if (newton_pair || j < nlocal) {
              uMech[j] -= factor_dpd;
            }
          }
        }
      }
      // store updated velocity for i
      v[i][0] = vxi;
      v[i][1] = vyi;
      v[i][2] = vzi;
    }

    // Communicate the ghost deltas to the atom owners
    if (idir > 1) comm->reverse_comm_fix(this);

  }  //End Loop over all directions For idir = Top, Top-Right, Right, Bottom-Right, Back

  memory->sfree(v_t0);
  v_t0 = NULL;
}

/* ---------------------------------------------------------------------- */

int FixShardlow::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int ii,jj,m;
  double **v  = atom->v;

  m = 0;
  for (ii = 0; ii < n; ii++) {
    jj = list[ii];
    buf[m++] = v[jj][0];
    buf[m++] = v[jj][1];
    buf[m++] = v[jj][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixShardlow::unpack_forward_comm(int n, int first, double *buf)
{
  int ii,m,last;
  int nlocal  = atom->nlocal;
  double **v  = atom->v;

  m = 0;
  last = first + n ;
  for (ii = first; ii < last; ii++) {
    v_t0[ii - nlocal][0] = v[ii][0] = buf[m++];
    v_t0[ii - nlocal][1] = v[ii][1] = buf[m++];
    v_t0[ii - nlocal][2] = v[ii][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int FixShardlow::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  int nlocal  = atom->nlocal;
  double **v  = atom->v;
  double *uCond = atom->uCond;
  double *uMech = atom->uMech;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = v[i][0] - v_t0[i - nlocal][0];
    buf[m++] = v[i][1] - v_t0[i - nlocal][1];
    buf[m++] = v[i][2] - v_t0[i - nlocal][2];
    if(pairDPDE){
      buf[m++] = uCond[i]; // for ghosts, this is an accumulated delta
      buf[m++] = uMech[i]; // for ghosts, this is an accumulated delta
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixShardlow::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  double **v  = atom->v;
  double *uCond = atom->uCond;
  double *uMech = atom->uMech;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];

    v[j][0] += buf[m++];
    v[j][1] += buf[m++];
    v[j][2] += buf[m++];
    if(pairDPDE){
      uCond[j] += buf[m++]; // add in the accumulated delta
      uMech[j] += buf[m++]; // add in the accumulated delta
    }
  }
}

/* ----------------------------------------------------------------------
   convert atom coords into the ssa active interaction region number
------------------------------------------------------------------------- */

int FixShardlow::coord2ssaAIR(double *x)
{
  int ix, iy, iz;

  ix = iy = iz = 0;
  if (x[2] < domain->sublo[2]) iz = -1;
  if (x[2] >= domain->subhi[2]) iz = 1;
  if (x[1] < domain->sublo[1]) iy = -1;
  if (x[1] >= domain->subhi[1]) iy = 1;
  if (x[0] < domain->sublo[0]) ix = -1;
  if (x[0] >= domain->subhi[0]) ix = 1;

  if(iz < 0){
    return 0;
  } else if(iz == 0){
    if( iy<0 ) return 0; // bottom left/middle/right
    if( (iy==0) && (ix<0)  ) return 0; // left atoms
    if( (iy==0) && (ix==0) ) return 1; // Locally owned atoms
    if( (iy==0) && (ix>0)  ) return 3; // Right atoms
    if( (iy>0)  && (ix==0) ) return 2; // Top-middle atoms
    if( (iy>0)  && (ix!=0) ) return 4; // Top-right and top-left atoms
  } else { // iz > 0
    if((ix==0) && (iy==0)) return 5; // Back atoms
    if((ix==0) && (iy!=0)) return 6; // Top-back and bottom-back atoms
    if((ix!=0) && (iy==0)) return 7; // Left-back and right-back atoms
    if((ix!=0) && (iy!=0)) return 8; // Back corner atoms
  }

  return 0;
}

/* ---------------------------------------------------------------------- */

void FixShardlow::grow_arrays(int nmax)
{
  memory->grow(atom->ssaAIR,nmax,"fix_shardlow:ssaAIR");
}

void FixShardlow::copy_arrays(int i, int j, int delflag)
{
  atom->ssaAIR[j] = atom->ssaAIR[i];
}

void FixShardlow::set_arrays(int i)
{
  atom->ssaAIR[i] = 1; /* coord2ssaAIR(x[i]) */
}

int FixShardlow::unpack_border(int n, int first, double *buf)
{
  int i,last = first + n;
  for (i = first; i < last; i++) {
    atom->ssaAIR[i] = coord2ssaAIR(atom->x[i]);
  }
  return 0;
}

int FixShardlow::unpack_exchange(int i, double *buf)
{
  atom->ssaAIR[i] = 1; /* coord2ssaAIR(x[i]) */
  return 0;
}

void FixShardlow::unpack_restart(int i, int nth)
{
  atom->ssaAIR[i] = 1; /* coord2ssaAIR(x[i]) */
}

double FixShardlow::memory_usage()
{
  double bytes = 0.0;
  bytes += memory->usage(atom->ssaAIR,atom->nmax);
  bytes += sizeof(double)*3*atom->nghost; // v_t0[]
  return bytes;
}

