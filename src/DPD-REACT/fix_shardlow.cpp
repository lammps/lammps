// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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

#include "fix_shardlow.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "npair.h"
#include "npair_half_bin_newton_ssa.h"
#include "pair_dpd_fdt.h"
#include "pair_dpd_fdt_energy.h"
#include "update.h"

#include <cstring>
#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace random_external_state;

static constexpr double EPSILON = 1.0e-10;
#define EPSILON_SQUARED ((EPSILON) * (EPSILON))

static const char cite_fix_shardlow[] =
  "fix shardlow command: doi:10.1016/j.cpc.2014.03.029, doi:10.1063/1.3660209\n\n"
  "@Article{Larentzos14,\n"
  " author = {J. P. Larentzos and J. K. Brennan and J. D. Moore and M. Lisal and W. D. Mattson},\n"
  " title = {Parallel Implementation of Isothermal and Isoenergetic Dissipative Particle Dynamics Using {S}hardlow-Like Splitting Algorithms},\n"
  " journal = {Comput.\\ Phys.\\ Commun.},\n"
  " year =    2014,\n"
  " volume =  185\n"
  " pages =   {1987--1998}\n"
  "}\n\n"
  "@Article{Lisal11,\n"
  " author = {M. Lisal and  J. K. Brennan and J. Bonet Avalos},\n"
  " title = {Dissipative Particle Dynamics at Isothermal, Isobaric, Isoenergetic, and Isoenthalpic Conditions Using {S}hardlow-Like Splitting Algorithms},\n"
  " journal = {J.~Chem.\\ Phys.},\n"
  " year =    2011,\n"
  " volume =  135,\n"
  " number =  20,\n"
  " pages =   {204105}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

FixShardlow::FixShardlow(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), pairDPD(nullptr), pairDPDE(nullptr), v_t0(nullptr)
  ,rand_state(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_shardlow);

  if (narg != 3) error->all(FLERR,"Illegal fix shardlow command");

  pairDPD = nullptr;
  pairDPDE = nullptr;
  pairDPD = dynamic_cast<PairDPDfdt *>(force->pair_match("dpd/fdt",1));
  pairDPDE = dynamic_cast<PairDPDfdtEnergy *>(force->pair_match("dpd/fdt/energy",1));
  if (pairDPDE == nullptr)
    pairDPDE = dynamic_cast<PairDPDfdtEnergy *>(force->pair_match("dpd/fdt/energy/kk",1));

  maxRNG = 0;
  if (pairDPDE) {
    comm_forward = 3;
    comm_reverse = 5;
  } else {
    comm_forward = 3;
    comm_reverse = 3;
  }

  if (pairDPD == nullptr && pairDPDE == nullptr)
    error->all(FLERR,"Must use pair_style dpd/fdt or dpd/fdt/energy with fix shardlow");

}

/* ---------------------------------------------------------------------- */

FixShardlow::~FixShardlow()
{
  memory->destroy(rand_state);
  maxRNG = 0;
}

/* ---------------------------------------------------------------------- */

int FixShardlow::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixShardlow::init()
{
  // SSA requires newton on
  neighbor->add_request(this, NeighConst::REQ_GHOST|NeighConst::REQ_NEWTON_ON|NeighConst::REQ_SSA);
}

/* ---------------------------------------------------------------------- */

void FixShardlow::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixShardlow::setup(int /*vflag*/)
{
  bool fixShardlow = false;

  for (int i = 0; i < modify->nfix; i++)
    if (strstr(modify->fix[i]->style,"nvt") || strstr(modify->fix[i]->style,"npt") ||
        strstr(modify->fix[i]->style,"gle") || strstr(modify->fix[i]->style,"gld"))
      error->all(FLERR,"Cannot use constant temperature integration routines with DPD-REACT.");

  for (int i = 0; i < modify->nfix; i++) {
    if (utils::strmatch(modify->fix[i]->style,"^shardlow")) fixShardlow = true;
    if (utils::strmatch(modify->fix[i]->style,"^nve") || utils::strmatch(modify->fix[i]->style,"^nph")) {
      if (fixShardlow) break;
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
  int start_ii,
  int count,
  int id
)
{
  es_RNG_t RNGstate = rand_state[id];
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

  int ct = count;
  int ii = start_ii;

while (ct-- > 0) {
  const int i = list->ilist[ii];
  const int *jlist = list->firstneigh[ii];
  const int jlen = list->numneigh[ii];
  ii++;
  if (jlen <= 0) continue;

  const double xtmp = x[i][0];
  const double ytmp = x[i][1];
  const double ztmp = x[i][2];

  // load velocity for i from memory
  double vxi = v[i][0];
  double vyi = v[i][1];
  double vzi = v[i][2];

  int itype = type[i];

  cut2_i = pairDPD->cutsq[itype];
  cut_i  = pairDPD->cut[itype];
  sigma_i = pairDPD->sigma[itype];
  theta_ij_inv = 1.0/pairDPD->temperature; // independent of i,j

  const double mass_i = (rmass) ? rmass[i] : mass[itype];
  const double massinv_i = 1.0 / mass_i;

#ifdef DEBUG_SSA_PAIR_CT
  const int nlocal = atom->nlocal;
#endif

  // Loop over Directional Neighbors only
  for (int jj = 0; jj < jlen; jj++) {
    int j = jlist[jj] & NEIGHMASK;
    int jtype = type[j];

    double delx = xtmp - x[j][0];
    double dely = ytmp - x[j][1];
    double delz = ztmp - x[j][2];
    double rsq = delx*delx + dely*dely + delz*delz;
#ifdef DEBUG_SSA_PAIR_CT
    if ((i < nlocal) && (j < nlocal)) ++(counters[0][0]);
    else ++(counters[0][1]);
    ++(counters[0][2]);
    int rsqi = rsq / 8;
    if (rsqi < 0) rsqi = 0;
    else if (rsqi > 31) rsqi = 31;
    ++(hist[rsqi]);
#endif

    // NOTE: r can be 0.0 in DPD systems, so do EPSILON_SQUARED test
    if ((rsq < cut2_i[jtype]) && (rsq >= EPSILON_SQUARED)) {
#ifdef DEBUG_SSA_PAIR_CT
      if ((i < nlocal) && (j < nlocal)) ++(counters[1][0]);
      else ++(counters[1][1]);
      ++(counters[1][2]);
#endif
      double r = sqrt(rsq);
      double rinv = 1.0/r;
      double delx_rinv = delx*rinv;
      double dely_rinv = dely*rinv;
      double delz_rinv = delz*rinv;

      double wr = 1.0 - r/cut_i[jtype];
      double wdt = wr*wr*dt;

      double halfsigma_ij = 0.5*sigma_i[jtype];
      double halfgamma_ij = halfsigma_ij*halfsigma_ij*boltz_inv*theta_ij_inv;

      double sigmaRand = halfsigma_ij*wr*dtsqrt*ftm2v * es_normal(RNGstate);

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

  rand_state[id] = RNGstate;
}

/* ----------------------------------------------------------------------
   Perform the stochastic integration and Shardlow update for constant energy
   Allow for both per-type and per-atom mass

   NOTE: only implemented for orthogonal boxes, not triclinic
------------------------------------------------------------------------- */
void FixShardlow::ssa_update_dpde(
  int start_ii,
  int count,
  int id
)
{
  es_RNG_t RNGstate = rand_state[id];
  double **x = atom->x;
  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  double *uCond = atom->uCond;
  double *uMech = atom->uMech;
  double *dpdTheta = atom->dpdTheta;

  double *cut_i, *cut2_i, *sigma_i, *kappa_i, *alpha_i;
  double theta_ij_inv, theta_i_inv;
  const double boltz_inv = 1.0/force->boltz;
  const double ftm2v = force->ftm2v;

  const double dt     = update->dt;

  int ct = count;
  int ii = start_ii;

while (ct-- > 0) {
  const int i = list->ilist[ii];
  const int *jlist = list->firstneigh[ii];
  const int jlen = list->numneigh[ii];
  ii++;
  if (jlen <= 0) continue;

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

  cut2_i = pairDPDE->cutsq[itype];
  cut_i  = pairDPDE->cut[itype];
  sigma_i = pairDPDE->sigma[itype];
  kappa_i = pairDPDE->kappa[itype];
  alpha_i = pairDPDE->alpha[itype];
  theta_i_inv = 1.0/dpdTheta[i];
  const double mass_i = (rmass) ? rmass[i] : mass[itype];
  const double massinv_i = 1.0 / mass_i;
  const double mass_i_div_neg4_ftm2v = mass_i*(-0.25)/ftm2v;

#ifdef DEBUG_SSA_PAIR_CT
  const int nlocal = atom->nlocal;
#endif

  // Loop over Directional Neighbors only
  for (int jj = 0; jj < jlen; jj++) {
    int j = jlist[jj] & NEIGHMASK;
    int jtype = type[j];

    double delx = xtmp - x[j][0];
    double dely = ytmp - x[j][1];
    double delz = ztmp - x[j][2];
    double rsq = delx*delx + dely*dely + delz*delz;
#ifdef DEBUG_SSA_PAIR_CT
    if ((i < nlocal) && (j < nlocal)) ++(counters[0][0]);
    else ++(counters[0][1]);
    ++(counters[0][2]);
    int rsqi = rsq / 8;
    if (rsqi < 0) rsqi = 0;
    else if (rsqi > 31) rsqi = 31;
    ++(hist[rsqi]);
#endif

    // NOTE: r can be 0.0 in DPD systems, so do EPSILON_SQUARED test
    if ((rsq < cut2_i[jtype]) && (rsq >= EPSILON_SQUARED)) {
#ifdef DEBUG_SSA_PAIR_CT
      if ((i < nlocal) && (j < nlocal)) ++(counters[1][0]);
      else ++(counters[1][1]);
      ++(counters[1][2]);
#endif
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

      double sigmaRand = halfsigma_ij*wr*dtsqrt*ftm2v * es_normal(RNGstate);

      double mass_j = (rmass) ? rmass[j] : mass[jtype];
      double mass_ij_div_neg4_ftm2v = mass_j*mass_i_div_neg4_ftm2v;
      double massinv_j = 1.0 / mass_j;

      // Compute uCond
      double kappa_ij = kappa_i[jtype];
      double alpha_ij = alpha_i[jtype];
      double del_uCond = alpha_ij*wr*dtsqrt * es_normal(RNGstate);

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

  rand_state[id] = RNGstate;
}

void FixShardlow::initial_integrate(int /*vflag*/)
{
  int ii;

  const int nlocal = atom->nlocal;
  const int nghost = atom->nghost;

  const bool useDPDE = (pairDPDE != nullptr);

  // NOTE: this logic is specific to orthogonal boxes, not triclinic

  // Enforce the constraint that ghosts must be contained in the nearest sub-domains
  double bbx = domain->subhi[0] - domain->sublo[0];
  double bby = domain->subhi[1] - domain->sublo[1];
  double bbz = domain->subhi[2] - domain->sublo[2];

  double rcut = 2.0*neighbor->cutneighmax;

  if (domain->triclinic)
    error->all(FLERR,"Fix shardlow does not yet support triclinic geometries");

  if (rcut >= bbx || rcut >= bby || rcut>= bbz )
    error->one(FLERR,"Shardlow algorithm requires sub-domain length > 2*(rcut+skin). "
                    "Either reduce the number of processors requested, or change the cutoff/skin: "
                    "rcut= {} bbx= {} bby= {} bbz= {}\n", rcut, bbx, bby, bbz);

  auto np_ssa = dynamic_cast<NPairHalfBinNewtonSSA*>(list->np);
  if (!np_ssa) error->one(FLERR, "NPair wasn't a NPairHalfBinNewtonSSA object");
  int ssa_phaseCt = np_ssa->ssa_phaseCt;
  int *ssa_phaseLen = np_ssa->ssa_phaseLen;
  int **ssa_itemLoc = np_ssa->ssa_itemLoc;
  int **ssa_itemLen = np_ssa->ssa_itemLen;
  int ssa_gphaseCt = np_ssa->ssa_gphaseCt;
  int *ssa_gphaseLen = np_ssa->ssa_gphaseLen;
  int **ssa_gitemLoc = np_ssa->ssa_gitemLoc;
  int **ssa_gitemLen = np_ssa->ssa_gitemLen;

  int maxWorkItemCt = np_ssa->ssa_maxPhaseLen;
  if (maxWorkItemCt > maxRNG) {
    uint64_t my_seed = comm->me + (useDPDE ? pairDPDE->seed : pairDPD->seed);
    es_RNG_t serial_rand_state;
    es_init(serial_rand_state, my_seed);

    memory->grow(rand_state, maxWorkItemCt, "FixShardlow:rand_state");
    for (int i = 0; i < maxWorkItemCt; ++i) {
      es_genNextParallelState(serial_rand_state, rand_state[i]);
    }

    maxRNG = maxWorkItemCt;
  }

#ifdef DEBUG_SSA_PAIR_CT
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j)
      counters[i][j] = 0;
  for (int i = 0; i < 32; ++i) hist[i] = 0;
#endif

  // Allocate memory for v_t0 to hold the initial velocities for the ghosts
  v_t0 = (double (*)[3]) memory->smalloc(sizeof(double)*3*nghost, "FixShardlow:v_t0");

  dtsqrt = sqrt(update->dt);

  // process neighbors in the local AIR
  for (int workPhase = 0; workPhase < ssa_phaseCt; ++workPhase) {
    int workItemCt = ssa_phaseLen[workPhase];

    for (int workItem = 0; workItem < workItemCt; ++workItem) {
      int ct = ssa_itemLen[workPhase][workItem];
      ii = ssa_itemLoc[workPhase][workItem];
      if (useDPDE) ssa_update_dpde(ii, ct, workItem);
      else ssa_update_dpd(ii, ct, workItem);
    }
  }

  //Loop over all 13 outward directions (7 stages)
  for (int workPhase = 0; workPhase < ssa_gphaseCt; ++workPhase) {
    int workItemCt = ssa_gphaseLen[workPhase];

    // Communicate the updated velocities to all nodes
    comm->forward_comm(this);

    if (useDPDE) {
      // Zero out the ghosts' uCond & uMech to be used as delta accumulators
      memset(&(atom->uCond[nlocal]), 0, sizeof(double)*nghost);
      memset(&(atom->uMech[nlocal]), 0, sizeof(double)*nghost);
    }

    for (int workItem = 0; workItem < workItemCt; ++workItem) {
      int ct = ssa_gitemLen[workPhase][workItem];
      ii = ssa_gitemLoc[workPhase][workItem];
      if (useDPDE) ssa_update_dpde(ii, ct, workItem);
      else ssa_update_dpd(ii, ct, workItem);
    }

    // Communicate the ghost deltas to the atom owners
    comm->reverse_comm(this);

  }  //End Loop over all directions For airnum = Top, Top-Right, Right, Bottom-Right, Back

#ifdef DEBUG_SSA_PAIR_CT
for (int i = 0; i < 32; ++i) fprintf(stdout, "%8d", hist[i]);
fprintf(stdout, "\n%6d %6d,%6d %6d: "
  ,counters[0][2]
  ,counters[1][2]
  ,counters[0][1]
  ,counters[1][1]
);
#endif

  memory->sfree(v_t0);
  v_t0 = nullptr;
}

/* ---------------------------------------------------------------------- */

int FixShardlow::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
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
    if (pairDPDE) {
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
    if (pairDPDE) {
      uCond[j] += buf[m++]; // add in the accumulated delta
      uMech[j] += buf[m++]; // add in the accumulated delta
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixShardlow::memory_usage()
{
  double bytes = 0.0;
  bytes += (double)sizeof(double)*3*atom->nghost; // v_t0[]
  bytes += (double)sizeof(*rand_state)*maxRNG; // rand_state[]
  return bytes;
}

