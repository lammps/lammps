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
#define EPSILON_SQUARED ((EPSILON) * (EPSILON))

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
  Fix(lmp, narg, arg), pairDPD(NULL), pairDPDE(NULL), v_t0(NULL)
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
  memset(atom->ssaAIR, 0, sizeof(int)*atom->nlocal);

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
  mask |= PRE_EXCHANGE | MIN_PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixShardlow::init()
{
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix  = 1;
  neighbor->requests[irequest]->ssa  = 1;
}

/* ---------------------------------------------------------------------- */

void FixShardlow::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixShardlow::pre_exchange()
{
  memset(atom->ssaAIR, 0, sizeof(int)*atom->nlocal);
}

/* ---------------------------------------------------------------------- */

void FixShardlow::setup_pre_exchange()
{
  memset(atom->ssaAIR, 0, sizeof(int)*atom->nlocal);
}

/* ---------------------------------------------------------------------- */

void FixShardlow::min_pre_exchange()
{
  memset(atom->ssaAIR, 0, sizeof(int)*atom->nlocal);
}

/* ---------------------------------------------------------------------- */

void FixShardlow::setup(int vflag)
{
  bool fixShardlow = false;

  for (int i = 0; i < modify->nfix; i++)
    if (strncmp(modify->fix[i]->style,"nvt",3) == 0 || strncmp(modify->fix[i]->style,"npt",3) == 0)
      error->all(FLERR,"Cannot use constant temperature integration routines with DPD.");

  for (int i = 0; i < modify->nfix; i++){
    if (strcmp(modify->fix[i]->style,"shardlow") == 0) fixShardlow = true;
    if (strncmp(modify->fix[i]->style,"nve",3) == 0 || (strncmp(modify->fix[i]->style,"nph",3) == 0)){
      if(fixShardlow) break;
      else error->all(FLERR,"The deterministic integrator must follow fix shardlow in the input file.");
    }
    if (i == modify->nfix-1) error->all(FLERR,"A deterministic integrator (e.g. fix nve or fix nph) is required when using fix shardlow.");
  }
}

/* ----------------------------------------------------------------------
   Perform the stochastic integration and Shardlow update for constant temperature
   Allow for both per-type and per-atom mass

   NOTE: only implemented for orthogonal boxes, not triclinic
------------------------------------------------------------------------- */
void FixShardlow::ssa_update_dpd(
  int i,
  int *jlist,
  int jlen
)
{
  class RanMars *pRNG;
  double **x = atom->x;
  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;

  double *cut_i, *cut2_i, *sigma_i;
  double theta_ij_inv;
  const double boltz_inv = 1.0/force->boltz;
  const double ftm2v = force->ftm2v;

  const double dt     = update->dt;

  const double xtmp = x[i][0];
  const double ytmp = x[i][1];
  const double ztmp = x[i][2];

  // load velocity for i from memory
  double vxi = v[i][0];
  double vyi = v[i][1];
  double vzi = v[i][2];

  int itype = type[i];

  pRNG = pairDPD->random;
  cut2_i = pairDPD->cutsq[itype];
  cut_i  = pairDPD->cut[itype];
  sigma_i = pairDPD->sigma[itype];
  theta_ij_inv = 1.0/pairDPD->temperature; // independent of i,j

  const double mass_i = (rmass) ? rmass[i] : mass[itype];
  const double massinv_i = 1.0 / mass_i;

  // Loop over Directional Neighbors only
  for (int jj = 0; jj < jlen; jj++) {
    int j = jlist[jj] & NEIGHMASK;
    int jtype = type[j];

    double delx = xtmp - x[j][0];
    double dely = ytmp - x[j][1];
    double delz = ztmp - x[j][2];
    double rsq = delx*delx + dely*dely + delz*delz;

    // NOTE: r can be 0.0 in DPD systems, so do EPSILON_SQUARED test
    if ((rsq < cut2_i[jtype]) && (rsq >= EPSILON_SQUARED)) {
      double r = sqrt(rsq);
      double rinv = 1.0/r;
      double delx_rinv = delx*rinv;
      double dely_rinv = dely*rinv;
      double delz_rinv = delz*rinv;

      double wr = 1.0 - r/cut_i[jtype];
      double wdt = wr*wr*dt;

      double halfsigma_ij = 0.5*sigma_i[jtype];
      double halfgamma_ij = halfsigma_ij*halfsigma_ij*boltz_inv*theta_ij_inv;

      double sigmaRand = halfsigma_ij*wr*dtsqrt*ftm2v * pRNG->gaussian();

      double mass_j = (rmass) ? rmass[j] : mass[jtype];
      double massinv_j = 1.0 / mass_j;

      double gammaFactor = halfgamma_ij*wdt*ftm2v;
      double inv_1p_mu_gammaFactor = 1.0/(1.0 + (massinv_i + massinv_j)*gammaFactor);

      double vxj = v[j][0];
      double vyj = v[j][1];
      double vzj = v[j][2];

      // Compute the initial velocity difference between atom i and atom j
      double delvx = vxi - vxj;
      double delvy = vyi - vyj;
      double delvz = vzi - vzj;
      double dot_rinv = (delx_rinv*delvx + dely_rinv*delvy + delz_rinv*delvz);

      // Compute momentum change between t and t+dt
      double factorA = sigmaRand - gammaFactor*dot_rinv;

      // Update the velocity on i
      vxi += delx_rinv*factorA*massinv_i;
      vyi += dely_rinv*factorA*massinv_i;
      vzi += delz_rinv*factorA*massinv_i;

      // Update the velocity on j
      vxj -= delx_rinv*factorA*massinv_j;
      vyj -= dely_rinv*factorA*massinv_j;
      vzj -= delz_rinv*factorA*massinv_j;

      //ii.   Compute the new velocity diff
      delvx = vxi - vxj;
      delvy = vyi - vyj;
      delvz = vzi - vzj;
      dot_rinv = delx_rinv*delvx + dely_rinv*delvy + delz_rinv*delvz;

      // Compute the new momentum change between t and t+dt
      double factorB = (sigmaRand - gammaFactor*dot_rinv)*inv_1p_mu_gammaFactor;

      // Update the velocity on i
      vxi += delx_rinv*factorB*massinv_i;
      vyi += dely_rinv*factorB*massinv_i;
      vzi += delz_rinv*factorB*massinv_i;

      // Update the velocity on j
      vxj -= delx_rinv*factorB*massinv_j;
      vyj -= dely_rinv*factorB*massinv_j;
      vzj -= delz_rinv*factorB*massinv_j;

      // Store updated velocity for j
      v[j][0] = vxj;
      v[j][1] = vyj;
      v[j][2] = vzj;
    }
  }
  // store updated velocity for i
  v[i][0] = vxi;
  v[i][1] = vyi;
  v[i][2] = vzi;
}

/* ----------------------------------------------------------------------
   Perform the stochastic integration and Shardlow update for constant energy
   Allow for both per-type and per-atom mass

   NOTE: only implemented for orthogonal boxes, not triclinic
------------------------------------------------------------------------- */
void FixShardlow::ssa_update_dpde(
  int i,
  int *jlist,
  int jlen
)
{
  class RanMars *pRNG;
  double **x = atom->x;
  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  double *uCond = atom->uCond;
  double *uMech = atom->uMech;
  double *dpdTheta = atom->dpdTheta;

  double *cut_i, *cut2_i, *sigma_i, *kappa_i;
  double theta_ij_inv, theta_i_inv;
  const double boltz2 = 2.0*force->boltz;
  const double boltz_inv = 1.0/force->boltz;
  const double ftm2v = force->ftm2v;

  const double dt     = update->dt;

  const double xtmp = x[i][0];
  const double ytmp = x[i][1];
  const double ztmp = x[i][2];

  // load velocity for i from memory
  double vxi = v[i][0];
  double vyi = v[i][1];
  double vzi = v[i][2];

  double uMech_i = uMech[i];
  double uCond_i = uCond[i];
  int itype = type[i];

  pRNG = pairDPDE->random;
  cut2_i = pairDPDE->cutsq[itype];
  cut_i  = pairDPDE->cut[itype];
  sigma_i = pairDPDE->sigma[itype];
  kappa_i = pairDPDE->kappa[itype];
  theta_i_inv = 1.0/dpdTheta[i];
  const double mass_i = (rmass) ? rmass[i] : mass[itype];
  const double massinv_i = 1.0 / mass_i;
  const double mass_i_div_neg4_ftm2v = mass_i*(-0.25)/ftm2v;

  // Loop over Directional Neighbors only
  for (int jj = 0; jj < jlen; jj++) {
    int j = jlist[jj] & NEIGHMASK;
    int jtype = type[j];

    double delx = xtmp - x[j][0];
    double dely = ytmp - x[j][1];
    double delz = ztmp - x[j][2];
    double rsq = delx*delx + dely*dely + delz*delz;

    // NOTE: r can be 0.0 in DPD systems, so do EPSILON_SQUARED test
    if ((rsq < cut2_i[jtype]) && (rsq >= EPSILON_SQUARED)) {
      double r = sqrt(rsq);
      double rinv = 1.0/r;
      double delx_rinv = delx*rinv;
      double dely_rinv = dely*rinv;
      double delz_rinv = delz*rinv;

      double wr = 1.0 - r/cut_i[jtype];
      double wdt = wr*wr*dt;

      // Compute the current temperature
      double theta_j_inv = 1.0/dpdTheta[j];
      theta_ij_inv = 0.5*(theta_i_inv + theta_j_inv);

      double halfsigma_ij = 0.5*sigma_i[jtype];
      double halfgamma_ij = halfsigma_ij*halfsigma_ij*boltz_inv*theta_ij_inv;

      double sigmaRand = halfsigma_ij*wr*dtsqrt*ftm2v * pRNG->gaussian();

      double mass_j = (rmass) ? rmass[j] : mass[jtype];
      double mass_ij_div_neg4_ftm2v = mass_j*mass_i_div_neg4_ftm2v;
      double massinv_j = 1.0 / mass_j;

      // Compute uCond
      double kappa_ij = kappa_i[jtype];
      double alpha_ij = sqrt(boltz2*kappa_ij);
      double del_uCond = alpha_ij*wr*dtsqrt * pRNG->gaussian();

      del_uCond += kappa_ij*(theta_i_inv - theta_j_inv)*wdt;
      uCond[j] -= del_uCond;
      uCond_i += del_uCond;

      double gammaFactor = halfgamma_ij*wdt*ftm2v;
      double inv_1p_mu_gammaFactor = 1.0/(1.0 + (massinv_i + massinv_j)*gammaFactor);

      double vxj = v[j][0];
      double vyj = v[j][1];
      double vzj = v[j][2];
      double dot4 = vxj*vxj + vyj*vyj + vzj*vzj;
      double dot3 = vxi*vxi + vyi*vyi + vzi*vzi;

      // Compute the initial velocity difference between atom i and atom j
      double delvx = vxi - vxj;
      double delvy = vyi - vyj;
      double delvz = vzi - vzj;
      double dot_rinv = (delx_rinv*delvx + dely_rinv*delvy + delz_rinv*delvz);

      // Compute momentum change between t and t+dt
      double factorA = sigmaRand - gammaFactor*dot_rinv;

      // Update the velocity on i
      vxi += delx_rinv*factorA*massinv_i;
      vyi += dely_rinv*factorA*massinv_i;
      vzi += delz_rinv*factorA*massinv_i;

      // Update the velocity on j
      vxj -= delx_rinv*factorA*massinv_j;
      vyj -= dely_rinv*factorA*massinv_j;
      vzj -= delz_rinv*factorA*massinv_j;

      //ii.   Compute the new velocity diff
      delvx = vxi - vxj;
      delvy = vyi - vyj;
      delvz = vzi - vzj;
      dot_rinv = delx_rinv*delvx + dely_rinv*delvy + delz_rinv*delvz;

      // Compute the new momentum change between t and t+dt
      double factorB = (sigmaRand - gammaFactor*dot_rinv)*inv_1p_mu_gammaFactor;

      // Update the velocity on i
      vxi += delx_rinv*factorB*massinv_i;
      vyi += dely_rinv*factorB*massinv_i;
      vzi += delz_rinv*factorB*massinv_i;
      double partial_uMech = (vxi*vxi + vyi*vyi + vzi*vzi - dot3)*massinv_j;

      // Update the velocity on j
      vxj -= delx_rinv*factorB*massinv_j;
      vyj -= dely_rinv*factorB*massinv_j;
      vzj -= delz_rinv*factorB*massinv_j;
      partial_uMech += (vxj*vxj + vyj*vyj + vzj*vzj - dot4)*massinv_i;

      // Store updated velocity for j
      v[j][0] = vxj;
      v[j][1] = vyj;
      v[j][2] = vzj;

      // Compute uMech
      double del_uMech = partial_uMech*mass_ij_div_neg4_ftm2v;
      uMech_i += del_uMech;
      uMech[j] += del_uMech;
    }
  }
  // store updated velocity for i
  v[i][0] = vxi;
  v[i][1] = vyi;
  v[i][2] = vzi;
  // store updated uMech and uCond for i
  uMech[i] = uMech_i;
  uCond[i] = uCond_i;
}

void FixShardlow::initial_integrate(int vflag)
{
  int i,ii,inum;
  int *ilist;

  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  int airnum;
  const bool useDPDE = (pairDPDE != NULL);

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

  inum = list->inum;
  ilist = list->ilist;

  dtsqrt = sqrt(update->dt);

  //Loop over all 14 directions (8 stages)
  for (airnum = 1; airnum <=8; airnum++){

    if (airnum > 1) {
      // Communicate the updated velocities to all nodes
      comm->forward_comm_fix(this);

      if(useDPDE){
        // Zero out the ghosts' uCond & uMech to be used as delta accumulators
        memset(&(atom->uCond[nlocal]), 0, sizeof(double)*nghost);
        memset(&(atom->uMech[nlocal]), 0, sizeof(double)*nghost);
      }
    }

    // Loop over neighbors of my atoms
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      int start = (airnum < 2) ? 0 : list->ndxAIR_ssa[i][airnum - 2];
      int len = list->ndxAIR_ssa[i][airnum - 1] - start;
      if (len > 0) {
        if (useDPDE) ssa_update_dpde(i, &(list->firstneigh[i][start]), len);
        else ssa_update_dpd(i, &(list->firstneigh[i][start]), len);
      }
    }

    // Communicate the ghost deltas to the atom owners
    if (airnum > 1) comm->reverse_comm_fix(this);

  }  //End Loop over all directions For airnum = Top, Top-Right, Right, Bottom-Right, Back

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
    return -1;
  } else if(iz == 0){
    if( iy<0 ) return -1; // bottom left/middle/right
    if( (iy==0) && (ix<0)  ) return -1; // left atoms
    if( (iy==0) && (ix==0) ) return 0; // Locally owned atoms
    if( (iy==0) && (ix>0)  ) return 3; // Right atoms
    if( (iy>0)  && (ix==0) ) return 2; // Top-middle atoms
    if( (iy>0)  && (ix!=0) ) return 4; // Top-right and top-left atoms
  } else { // iz > 0
    if((ix==0) && (iy==0)) return 5; // Back atoms
    if((ix==0) && (iy!=0)) return 6; // Top-back and bottom-back atoms
    if((ix!=0) && (iy==0)) return 7; // Left-back and right-back atoms
    if((ix!=0) && (iy!=0)) return 8; // Back corner atoms
  }

  return -2;
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
  atom->ssaAIR[i] = 0; /* coord2ssaAIR(x[i]) */
}

int FixShardlow::pack_border(int n, int *list, double *buf)
{
  for (int i = 0; i < n; i++) {
    int j = list[i];
    if (atom->ssaAIR[j] == 0) atom->ssaAIR[j] = 1; // not purely local anymore
  }
  return 0;
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
  atom->ssaAIR[i] = 0; /* coord2ssaAIR(x[i]) */
  return 0;
}

void FixShardlow::unpack_restart(int i, int nth)
{
  atom->ssaAIR[i] = 0; /* coord2ssaAIR(x[i]) */
}

double FixShardlow::memory_usage()
{
  double bytes = 0.0;
  bytes += memory->usage(atom->ssaAIR,atom->nmax);
  bytes += sizeof(double)*3*atom->nghost; // v_t0[]
  return bytes;
}

