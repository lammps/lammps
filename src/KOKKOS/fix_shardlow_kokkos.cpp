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

#include "fix_shardlow_kokkos.h"

#include "atom.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "npair_ssa_kokkos.h"
#include "pair_dpd_fdt_energy_kokkos.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace random_external_state;

#define EPSILON 1.0e-10
#define EPSILON_SQUARED ((EPSILON) * (EPSILON))

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixShardlowKokkos<DeviceType>::FixShardlowKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixShardlow(lmp, narg, arg), k_pairDPDE(nullptr), ghostmax(0), nlocal(0) , nghost(0)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  if (narg != 3) error->all(FLERR,"Illegal fix shardlow command");

//  k_pairDPD = nullptr;
  k_pairDPDE = nullptr;
//  k_pairDPD = (PairDPDfdtKokkos *) force->pair_match("dpd/fdt",1);
  k_pairDPDE = dynamic_cast<PairDPDfdtEnergyKokkos<DeviceType> *>(force->pair_match("dpd/fdt/energy",0));

//   if (k_pairDPDE) {
    comm_forward = 3;
    comm_reverse = 5;
//   } else {
//     comm_forward = 3;
//     comm_reverse = 3;
//   }


  if (/* k_pairDPD == nullptr &&*/ k_pairDPDE == nullptr)
    error->all(FLERR,"Must use pair_style "/*"dpd/fdt/kk or "*/"dpd/fdt/energy/kk with fix shardlow/kk");

#ifdef DEBUG_SSA_PAIR_CT
  d_counters = typename AT::t_int_2d("FixShardlowKokkos::d_counters", 2, 3);
  d_hist = typename AT::t_int_1d("FixShardlowKokkos::d_hist", 32);
  h_counters = Kokkos::create_mirror_view(d_counters);
  h_hist = Kokkos::create_mirror_view(d_hist);
#endif
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixShardlowKokkos<DeviceType>::~FixShardlowKokkos()
{
  ghostmax = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixShardlowKokkos<DeviceType>::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE | PRE_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixShardlowKokkos<DeviceType>::init()
{
  FixShardlow::init();

  // adjust neighbor list request for KOKKOS

  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);

  int ntypes = atom->ntypes;
  k_params = Kokkos::DualView<params_ssa**,Kokkos::LayoutRight,DeviceType>
    ("FixShardlowKokkos::params",ntypes+1,ntypes+1);
  params = k_params.template view<DeviceType>();
  k_pairDPDE->k_cutsq.template sync<DeviceType>();
  d_cutsq = k_pairDPDE->k_cutsq.template view<DeviceType>();

  for (int i = 1; i <= ntypes; i++) {
    for (int j = i; j <= ntypes; j++) {
      F_FLOAT cutone = k_pairDPDE->cut[i][j];
      if (cutone > EPSILON) k_params.h_view(i,j).cutinv = 1.0/cutone;
      else k_params.h_view(i,j).cutinv = FLT_MAX;
      k_params.h_view(i,j).halfsigma = 0.5*k_pairDPDE->sigma[i][j];
      k_params.h_view(i,j).kappa = k_pairDPDE->kappa[i][j];
      k_params.h_view(i,j).alpha = k_pairDPDE->alpha[i][j];

      k_params.h_view(j,i) = k_params.h_view(i,j);

      if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
        m_params[i][j] = m_params[j][i] = k_params.h_view(i,j);
        m_cutsq[j][i] = m_cutsq[i][j] = k_pairDPDE->k_cutsq.h_view(i,j);
      }
    }
  }

  k_params.template modify<LMPHostType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixShardlowKokkos<DeviceType>::init_list(int id, NeighList *ptr)
{
  FixShardlow::init_list(id, ptr);
  k_list = static_cast<NeighListKokkos<DeviceType>*>(ptr);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixShardlowKokkos<DeviceType>::pre_neighbor()
{
  // NOTE: this logic is specific to orthogonal boxes, not triclinic

  // Enforce the constraint that ghosts must be contained in the nearest sub-domains
  double bbx = domain->subhi[0] - domain->sublo[0];
  double bby = domain->subhi[1] - domain->sublo[1];
  double bbz = domain->subhi[2] - domain->sublo[2];

  double rcut = 2.0*neighbor->cutneighmax;

  if (domain->triclinic)
    error->all(FLERR,"Fix shardlow does not yet support triclinic geometries");

  if (rcut >= bbx || rcut >= bby || rcut>= bbz )
  {
    char fmt[] = {"Shardlow algorithm requires sub-domain length > 2*(rcut+skin). Either reduce the number of processors requested, or change the cutoff/skin: rcut= %e bbx= %e bby= %e bbz= %e\n"};
    char *msg = (char *) malloc(sizeof(fmt) + 4*15);
    sprintf(msg, fmt, rcut, bbx, bby, bbz);
    error->one(FLERR, msg);
  }

  nlocal = atomKK->nlocal;
  nghost = atomKK->nghost;

  // Allocate memory for h_v_t0 to hold the initial velocities for the ghosts
  if (nghost > ghostmax) {
    ghostmax = nghost;
    k_v_t0 = DAT::tdual_v_array("FixShardlowKokkos:v_t0", ghostmax);
    // d_v_t0 = k_v_t0.template view<DeviceType>();
    h_v_t0 = k_v_t0.h_view;
  }

  // Setup views of relevant data
  x = atomKK->k_x.template view<DeviceType>();
  v = atomKK->k_v.template view<DeviceType>();
  h_v = atomKK->k_v.h_view;
  uCond = atomKK->k_uCond.template view<DeviceType>();
  h_uCond = atomKK->k_uCond.h_view;
  uMech = atomKK->k_uMech.template view<DeviceType>();
  h_uMech = atomKK->k_uMech.h_view;
  type = atomKK->k_type.view<DeviceType>();
  if (atomKK->rmass) {
    massPerI = true;
    masses = atomKK->k_rmass.view<DeviceType>();
  } else {
    massPerI = false;
    masses = atomKK->k_mass.view<DeviceType>();
  }
//   if (k_pairDPDE) {
  dpdTheta = atomKK->k_dpdTheta.view<DeviceType>();

//} else {
//}
}

template<class DeviceType>
void FixShardlowKokkos<DeviceType>::setup_pre_neighbor()
{
  pre_neighbor();
}

/* ---------------------------------------------------------------------- */

#ifdef ENABLE_KOKKOS_DPD_CONSTANT_TEMPERATURE
#error "FixShardlowKokkos::ssa_update_dpd() is not functional yet - TIM 20170830"
/* ----------------------------------------------------------------------
   Perform the stochastic integration and Shardlow update for constant temperature
   Allow for both per-type and per-atom mass

   NOTE: only implemented for orthogonal boxes, not triclinic
------------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS>
void FixShardlowKokkos<DeviceType>::ssa_update_dpd(
  int start_ii, int count, int id
)
{
  es_RNG_t RNGstate = d_rand_state(id);

  int ct = count;
  int ii = start_ii;

  while (ct-- > 0) {
    const int i = d_ilist(ii);
    const int jlen = d_numneigh(ii);

    const double xtmp = x(i, 0);
    const double ytmp = x(i, 1);
    const double ztmp = x(i, 2);

    // load velocity for i from memory
    double vxi = v(i, 0);
    double vyi = v(i, 1);
    double vzi = v(i, 2);

    const int itype = type(i);

    const double mass_i = masses(massPerI ? i : itype);
    const double massinv_i = 1.0 / mass_i;

    // Loop over Directional Neighbors only
    for (int jj = 0; jj < jlen; jj++) {
      const int j = d_neighbors(ii,jj) & NEIGHMASK;
      int jtype = type[j];

      const X_FLOAT delx = xtmp - x(j, 0);
      const X_FLOAT dely = ytmp - x(j, 1);
      const X_FLOAT delz = ztmp - x(j, 2);
      const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
#ifdef DEBUG_SSA_PAIR_CT
      if ((i < nlocal) && (j < nlocal)) Kokkos::atomic_increment(&(d_counters(0, 0)));
      else Kokkos::atomic_increment(&(d_counters(0, 1)));
      Kokkos::atomic_increment(&(d_counters(0, 2)));
      int rsqi = rsq / 8;
      if (rsqi < 0) rsqi = 0;
      else if (rsqi > 31) rsqi = 31;
      Kokkos::atomic_increment(&(d_hist(rsqi)));
#endif

      // NOTE: r can be 0.0 in DPD systems, so do EPSILON_SQUARED test
      if ((rsq < (STACKPARAMS?m_cutsq[itype][jtype]:d_cutsq(itype,jtype)))
        && (rsq >= EPSILON_SQUARED)) {
#ifdef DEBUG_SSA_PAIR_CT
        if ((i < nlocal) && (j < nlocal)) Kokkos::atomic_increment(&(d_counters(1, 0)));
        else Kokkos::atomic_increment(&(d_counters(1, 1)));
        Kokkos::atomic_increment(&(d_counters(1, 2)));
#endif
        double r = sqrt(rsq);
        double rinv = 1.0/r;
        double delx_rinv = delx*rinv;
        double dely_rinv = dely*rinv;
        double delz_rinv = delz*rinv;

        double wr = 1.0 - r*(STACKPARAMS?m_params[itype][jtype].cutinv:params(itype,jtype).cutinv);
        double wdt = wr*wr*dt;

        double halfsigma_ij = STACKPARAMS?m_params[itype][jtype].halfsigma:params(itype,jtype).halfsigma;
        double halfgamma_ij = halfsigma_ij*halfsigma_ij*boltz_inv*theta_ij_inv;

        double sigmaRand = halfsigma_ij*wr*dtsqrt*ftm2v * es_normal(RNGstate);

        const double mass_j = masses(massPerI ? j : jtype);
        double massinv_j = 1.0 / mass_j;

        double gammaFactor = halfgamma_ij*wdt*ftm2v;
        double inv_1p_mu_gammaFactor = 1.0/(1.0 + (massinv_i + massinv_j)*gammaFactor);

        double vxj = v(j, 0);
        double vyj = v(j, 1);
        double vzj = v(j, 2);

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
        v(j, 0) = vxj;
        v(j, 1) = vyj;
        v(j, 2) = vzj;
      }
    }
    // store updated velocity for i
    v(i, 0) = vxi;
    v(i, 1) = vyi;
    v(i, 2) = vzi;
  }

  d_rand_state(id) = RNGstate;
}
#endif

/* ----------------------------------------------------------------------
   Perform the stochastic integration and Shardlow update for constant energy
   Allow for both per-type and per-atom mass

   NOTE: only implemented for orthogonal boxes, not triclinic
------------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS>
KOKKOS_INLINE_FUNCTION
void FixShardlowKokkos<DeviceType>::ssa_update_dpde(
  int start_ii, int count, int id
) const
{
  es_RNG_t RNGstate = d_rand_state(id);

  int ct = count;
  int ii = start_ii;

  while (ct-- > 0) {
    const int i = d_ilist(ii);
    const int jlen = d_numneigh(ii);

    const double xtmp = x(i, 0);
    const double ytmp = x(i, 1);
    const double ztmp = x(i, 2);

    // load velocity for i from memory
    double vxi = v(i, 0);
    double vyi = v(i, 1);
    double vzi = v(i, 2);

    double uMech_i = uMech(i);
    double uCond_i = uCond(i);
    const int itype = type(i);

    const double theta_i_inv = 1.0/dpdTheta(i);
    const double mass_i = masses(massPerI ? i : itype);
    const double massinv_i = 1.0 / mass_i;
    const double mass_i_div_neg4_ftm2v = mass_i*(-0.25)/ftm2v;

    // Loop over Directional Neighbors only
    for (int jj = 0; jj < jlen; jj++) {
      const int j = d_neighbors(ii,jj) & NEIGHMASK;
      const int jtype = type(j);

      const X_FLOAT delx = xtmp - x(j, 0);
      const X_FLOAT dely = ytmp - x(j, 1);
      const X_FLOAT delz = ztmp - x(j, 2);
      const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
#ifdef DEBUG_SSA_PAIR_CT
      if ((i < nlocal) && (j < nlocal)) Kokkos::atomic_increment(&(d_counters(0, 0)));
      else Kokkos::atomic_increment(&(d_counters(0, 1)));
      Kokkos::atomic_increment(&(d_counters(0, 2)));
      int rsqi = rsq / 8;
      if (rsqi < 0) rsqi = 0;
      else if (rsqi > 31) rsqi = 31;
      Kokkos::atomic_increment(&(d_hist(rsqi)));
#endif

      // NOTE: r can be 0.0 in DPD systems, so do EPSILON_SQUARED test
      if ((rsq < (STACKPARAMS?m_cutsq[itype][jtype]:d_cutsq(itype,jtype)))
        && (rsq >= EPSILON_SQUARED)) {
#ifdef DEBUG_SSA_PAIR_CT
        if ((i < nlocal) && (j < nlocal)) Kokkos::atomic_increment(&(d_counters(1, 0)));
        else Kokkos::atomic_increment(&(d_counters(1, 1)));
        Kokkos::atomic_increment(&(d_counters(1, 2)));
#endif

        double r = sqrt(rsq);
        double rinv = 1.0/r;
        double delx_rinv = delx*rinv;
        double dely_rinv = dely*rinv;
        double delz_rinv = delz*rinv;

        double wr = 1.0 - r*(STACKPARAMS?m_params[itype][jtype].cutinv:params(itype,jtype).cutinv);
        double wdt = wr*wr*dt;

        // Compute the current temperature
        double theta_j_inv = 1.0/dpdTheta(j);
        double theta_ij_inv = 0.5*(theta_i_inv + theta_j_inv);

        double halfsigma_ij = STACKPARAMS?m_params[itype][jtype].halfsigma:params(itype,jtype).halfsigma;
        double halfgamma_ij = halfsigma_ij*halfsigma_ij*boltz_inv*theta_ij_inv;

        double sigmaRand = halfsigma_ij*wr*dtsqrt*ftm2v * es_normal(RNGstate);

        const double mass_j = masses(massPerI ? j : jtype);
        double mass_ij_div_neg4_ftm2v = mass_j*mass_i_div_neg4_ftm2v;
        double massinv_j = 1.0 / mass_j;

        // Compute uCond
        double kappa_ij = STACKPARAMS?m_params[itype][jtype].kappa:params(itype,jtype).kappa;
        double alpha_ij = STACKPARAMS?m_params[itype][jtype].alpha:params(itype,jtype).alpha;
        double del_uCond = alpha_ij*wr*dtsqrt * es_normal(RNGstate);

        del_uCond += kappa_ij*(theta_i_inv - theta_j_inv)*wdt;
        uCond[j] -= del_uCond;
        uCond_i += del_uCond;

        double gammaFactor = halfgamma_ij*wdt*ftm2v;
        double inv_1p_mu_gammaFactor = 1.0/(1.0 + (massinv_i + massinv_j)*gammaFactor);

        double vxj = v(j, 0);
        double vyj = v(j, 1);
        double vzj = v(j, 2);
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
        v(j, 0) = vxj;
        v(j, 1) = vyj;
        v(j, 2) = vzj;

        // Compute uMech
        double del_uMech = partial_uMech*mass_ij_div_neg4_ftm2v;
        uMech_i += del_uMech;
        uMech(j) += del_uMech;
      }
    }
    // store updated velocity for i
    v(i, 0) = vxi;
    v(i, 1) = vyi;
    v(i, 2) = vzi;
    // store updated uMech and uCond for i
    uMech(i) = uMech_i;
    uCond(i) = uCond_i;
    ii++;
  }

  d_rand_state(id) = RNGstate;
}


template<class DeviceType>
void FixShardlowKokkos<DeviceType>::initial_integrate(int /*vflag*/)
{
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  copymode = 1;

  dtsqrt = sqrt(update->dt);

  NPairSSAKokkos<DeviceType> *np_ssa = dynamic_cast<NPairSSAKokkos<DeviceType>*>(list->np);
  if (!np_ssa) error->one(FLERR, "NPair wasn't a NPairSSAKokkos object");
  ssa_phaseCt = np_ssa->ssa_phaseCt;
  ssa_phaseLen = np_ssa->ssa_phaseLen;
  ssa_itemLoc = np_ssa->ssa_itemLoc;
  ssa_itemLen = np_ssa->ssa_itemLen;
  ssa_gphaseCt = np_ssa->ssa_gphaseCt;
  ssa_gphaseLen = np_ssa->ssa_gphaseLen;
  ssa_gitemLoc = np_ssa->ssa_gitemLoc;
  ssa_gitemLen = np_ssa->ssa_gitemLen;

  np_ssa->k_ssa_itemLoc.template sync<DeviceType>();
  np_ssa->k_ssa_itemLen.template sync<DeviceType>();
  np_ssa->k_ssa_gitemLoc.template sync<DeviceType>();
  np_ssa->k_ssa_gitemLen.template sync<DeviceType>();

  np_ssa->k_ssa_phaseLen.template sync<LMPHostType>();
  np_ssa->k_ssa_gphaseLen.template sync<LMPHostType>();
  auto h_ssa_phaseLen = np_ssa->k_ssa_phaseLen.h_view;
  auto h_ssa_gphaseLen = np_ssa->k_ssa_gphaseLen.h_view;

  int maxWorkItemCt = (int) ssa_itemLoc.extent(1);
  if (maxWorkItemCt < (int) ssa_gitemLoc.extent(1)) {
    maxWorkItemCt = (int) ssa_gitemLoc.extent(1);
  }
  if (maxWorkItemCt > maxRNG) {
    es_RNG_t serial_rand_state;
    es_init(serial_rand_state, pairDPDE->seed + comm->me);

    d_rand_state = es_RNGs_type("Kokkos::fix_shardlow::rand_state",maxWorkItemCt);
    typename es_RNGs_type::HostMirror h_rand_state = create_mirror_view(d_rand_state);
    for (int i = 0; i < maxWorkItemCt; ++i) {
      es_genNextParallelState(serial_rand_state, h_rand_state(i));
    }
    Kokkos::deep_copy(d_rand_state,h_rand_state);

    maxRNG = maxWorkItemCt;
  }

#ifdef DEBUG_SSA_PAIR_CT
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j)
      h_counters(i,j) = 0;
  for (int i = 0; i < 32; ++i) h_hist[i] = 0;
  Kokkos::deep_copy(d_counters, h_counters);
  Kokkos::deep_copy(d_hist, h_hist);
#endif

  //theta_ij_inv = 1.0/k_pairDPD->temperature; // independent of i,j
  boltz_inv = 1.0/force->boltz;
  ftm2v = force->ftm2v;
  dt     = update->dt;

  k_params.template sync<DeviceType>();

  // process neighbors in the local AIR
  atomKK->sync(execution_space,X_MASK | V_MASK | TYPE_MASK | RMASS_MASK | UCOND_MASK | UMECH_MASK | DPDTHETA_MASK);
  for (workPhase = 0; workPhase < ssa_phaseCt; ++workPhase) {
    int workItemCt = h_ssa_phaseLen[workPhase];

    if (atom->ntypes > MAX_TYPES_STACKPARAMS)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixShardlowSSAUpdateDPDE<false> >(0,workItemCt),*this);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixShardlowSSAUpdateDPDE<true> >(0,workItemCt),*this);
  }
  atomKK->modified(execution_space,V_MASK | UCOND_MASK | UMECH_MASK);

  //Loop over all 13 outward directions (7 stages)
  for (workPhase = 0; workPhase < ssa_gphaseCt; ++workPhase) {
    // int airnum = workPhase + 1;
    int workItemCt = h_ssa_gphaseLen[workPhase];

    // Communicate the updated velocities to all nodes
    atomKK->sync(Host,V_MASK);
    comm->forward_comm(this);
    atomKK->modified(Host,V_MASK);

    if (k_pairDPDE) {
      // Zero out the ghosts' uCond & uMech to be used as delta accumulators
//      memset(&(atom->uCond[nlocal]), 0, sizeof(double)*nghost);
//      memset(&(atom->uMech[nlocal]), 0, sizeof(double)*nghost);

      // must capture local variables, not class variables
      atomKK->sync(execution_space,UCOND_MASK | UMECH_MASK);
      auto l_uCond = uCond;
      auto l_uMech = uMech;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(nlocal,nlocal+nghost),
       LAMMPS_LAMBDA (const int i) {
        l_uCond(i) = 0.0;
        l_uMech(i) = 0.0;
      });
      atomKK->modified(execution_space,UCOND_MASK | UMECH_MASK);
    }

    // process neighbors in this AIR
    atomKK->sync(execution_space,X_MASK | V_MASK | TYPE_MASK | RMASS_MASK | UCOND_MASK | UMECH_MASK | DPDTHETA_MASK);
    if (atom->ntypes > MAX_TYPES_STACKPARAMS)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixShardlowSSAUpdateDPDEGhost<false> >(0,workItemCt),*this);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixShardlowSSAUpdateDPDEGhost<true> >(0,workItemCt),*this);
    atomKK->modified(execution_space,V_MASK | UCOND_MASK | UMECH_MASK);

    // Communicate the ghost deltas to the atom owners
    atomKK->sync(Host,V_MASK | UCOND_MASK | UMECH_MASK);
    comm->reverse_comm(this);
    atomKK->modified(Host,V_MASK | UCOND_MASK | UMECH_MASK);

  }  //End Loop over all directions For airnum = Top, Top-Right, Right, Bottom-Right, Back

#ifdef DEBUG_SSA_PAIR_CT
deep_copy(h_counters, d_counters);
deep_copy(h_hist, d_hist);
for (int i = 0; i < 32; ++i) fprintf(stdout, "%8d", h_hist[i]);
fprintf(stdout, "\n%6d %6d,%6d %6d: "
  ,h_counters(0, 2)
  ,h_counters(1, 2)
  ,h_counters(0, 1)
  ,h_counters(1, 1)
);
#endif

  copymode = 0;
}

template<class DeviceType>
template<bool STACKPARAMS>
KOKKOS_INLINE_FUNCTION
void FixShardlowKokkos<DeviceType>::operator()(TagFixShardlowSSAUpdateDPDE<STACKPARAMS>, const int &workItem) const {
  const int ct = ssa_itemLen(workPhase, workItem);
  const int ii = ssa_itemLoc(workPhase, workItem);
  ssa_update_dpde<STACKPARAMS>(ii, ct, workItem);
}

template<class DeviceType>
template<bool STACKPARAMS>
KOKKOS_INLINE_FUNCTION
void FixShardlowKokkos<DeviceType>::operator()(TagFixShardlowSSAUpdateDPDEGhost<STACKPARAMS>, const int &workItem) const {
  const int ct = ssa_gitemLen(workPhase, workItem);
  const int ii = ssa_gitemLoc(workPhase, workItem);
  ssa_update_dpde<STACKPARAMS>(ii, ct, workItem);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixShardlowKokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int ii,jj,m;

  m = 0;
  for (ii = 0; ii < n; ii++) {
    jj = list[ii];
    buf[m++] = h_v(jj, 0);
    buf[m++] = h_v(jj, 1);
    buf[m++] = h_v(jj, 2);
  }
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixShardlowKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  int ii,m,last;

  m = 0;
  last = first + n ;
  for (ii = first; ii < last; ii++) {
    h_v_t0(ii - nlocal, 0) = h_v(ii, 0) = buf[m++];
    h_v_t0(ii - nlocal, 1) = h_v(ii, 1) = buf[m++];
    h_v_t0(ii - nlocal, 2) = h_v(ii, 2) = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixShardlowKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = h_v(i, 0) - h_v_t0(i - nlocal, 0);
    buf[m++] = h_v(i, 1) - h_v_t0(i - nlocal, 1);
    buf[m++] = h_v(i, 2) - h_v_t0(i - nlocal, 2);
    if (k_pairDPDE) {
      buf[m++] = h_uCond(i); // for ghosts, this is an accumulated delta
      buf[m++] = h_uMech(i); // for ghosts, this is an accumulated delta
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixShardlowKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];

    h_v(j, 0) += buf[m++];
    h_v(j, 1) += buf[m++];
    h_v(j, 2) += buf[m++];
    if (k_pairDPDE) {
      h_uCond(j) += buf[m++]; // add in the accumulated delta
      h_uMech(j) += buf[m++]; // add in the accumulated delta
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double FixShardlowKokkos<DeviceType>::memory_usage()
{
  double bytes = 0.0;
  bytes += (double)sizeof(double)*3*ghostmax; // v_t0[]
  return bytes;
}

namespace LAMMPS_NS {
template class FixShardlowKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixShardlowKokkos<LMPHostType>;
#endif
}
