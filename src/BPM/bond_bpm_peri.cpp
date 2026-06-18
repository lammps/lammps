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
   Contributing author: Claude Opus 4.8 (Anthropic), under the direction of
   Joel Clemmer (SNL) and Axel Kohlmeyer (Temple U).

   Re-implementation of the four PERI package constitutive models on the BPM
   framework, derived from pair_peri.cpp and pair_peri_{pmb,lps,ves,eps}.cpp.
   Original PERI models: PMB and LPS (and the shared dilatation/weighted-volume
   infrastructure) by Mike Parks (SNL); VES and EPS by Rezwanur Rahman and
   J. T. Foster (UTSA).
------------------------------------------------------------------------- */

#include "bond_bpm_peri.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_bond_history.h"
#include "force.h"
#include "lattice.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cfloat>
#include <cmath>
#include <cstring>

static constexpr double EPSILON = 1e-10;
static constexpr double NEAR_ZERO = 2.2204e-16;

using namespace LAMMPS_NS;

static const char cite_bpm_peri[] =
    "BPM peridynamics (bond_style bpm/peri), derived from the PERI package: "
    "https://doi.org/10.1016/j.cpc.2008.06.011\n\n"
    "@Article{Parks08,\n"
    " author = {M. L. Parks and R. B. Lehoucq and S. J. Plimpton and S. A. Silling},\n"
    " title = {Implementing Peridynamics Within a Molecular Dynamics Code},\n"
    " journal = {Comput.\\ Phys.\\ Commun.},\n"
    " year =    2008,\n"
    " volume =  179,\n"
    " number =  11,\n"
    " pages =   {777--783}\n"
    "}\n\n";

// the ves/eps models additionally derive from the state-based viscoplasticity
// theory of Foster, Silling, and Chen (the original PDLAMMPS ves/eps models were
// implemented by Rahman and Foster, UTSA)
static const char cite_bpm_peri_plastic[] =
    "BPM peridynamics eps/ves models: https://doi.org/10.1002/nme.2725\n\n"
    "@Article{Foster10,\n"
    " author = {J. T. Foster and S. A. Silling and W. W. Chen},\n"
    " title = {Viscoplasticity Using Peridynamics},\n"
    " journal = {Int.\\ J.\\ Numer.\\ Methods Eng.},\n"
    " year =    2010,\n"
    " volume =  81,\n"
    " number =  10,\n"
    " pages =   {1242--1258}\n"
    "}\n\n";

// Linear partial-volume taper (Silling 2007): a bond whose reference length r0
// is within half a lattice spacing of the horizon delta "sees" only part of the
// partner's nodal volume.  Returns 1 well inside the horizon and 0.5 at delta.
static inline double vfrac_taper(double r0, double delta, double half_lc)
{
  if (fabs(r0 - delta) <= half_lc)
    return (-1.0 / (2.0 * half_lc)) * r0 + (1.0 + (delta - half_lc) / (2.0 * half_lc));
  return 1.0;
}

/* ---------------------------------------------------------------------- */

BondBPMPeri::BondBPMPeri(LAMMPS *_lmp) :
    BondBPM(_lmp), model(nullptr), c(nullptr), kbulk(nullptr), gshear(nullptr), lambda(nullptr),
    tau(nullptr), yieldstress(nullptr), cut(nullptr), s00(nullptr), alpha(nullptr),
    id_fix_property_peri(nullptr), index_vfrac(-1), index_s0(-1), index_smin(-1), index_lambda(-1),
    index_vinter(-1), index_dtheta(-1),
    smin_new(nullptr), s0_new(nullptr), nmax(0), state_based(0), wvolume_setup(0), kbulk_rep(0.0),
    wvolume(nullptr), theta(nullptr), winv(nullptr), commflag(COMM_SMIN), plastic(0), tdnorm(nullptr),
    pointwise_yield(0.0), gshear_rep(0.0)
{
  if (lmp->citeme) lmp->citeme->add(cite_bpm_peri);

  partial_flag = 1;
  writedata = 0;

  // one reference value per bond (r0); static (does not evolve)
  nhistory = 1;
  update_flag = 0;
  id_fix_bond_history = utils::strdup("HISTORY_BPM_PERI");

  // forward-communicate smin each step so the break test can read a neighbor's
  // previous-step smin even when that neighbor is a ghost
  comm_forward = 1;

  single_extra = 1;
  svector = new double[1];
}

/* ---------------------------------------------------------------------- */

BondBPMPeri::~BondBPMPeri()
{
  delete[] svector;

  memory->destroy(smin_new);
  memory->destroy(s0_new);
  memory->destroy(wvolume);
  memory->destroy(theta);
  memory->destroy(winv);
  memory->destroy(tdnorm);

  if (id_fix_property_peri) {
    if (modify->nfix) modify->delete_fix(id_fix_property_peri);
    delete[] id_fix_property_peri;
  }

  if (allocated) {
    memory->destroy(model);
    memory->destroy(c);
    memory->destroy(kbulk);
    memory->destroy(gshear);
    memory->destroy(lambda);
    memory->destroy(tau);
    memory->destroy(yieldstress);
    memory->destroy(cut);
    memory->destroy(s00);
    memory->destroy(alpha);
    memory->destroy(setflag);
  }
}

/* ----------------------------------------------------------------------
   store the reference bond length r0 for every bond, called once
------------------------------------------------------------------------- */

void BondBPMPeri::store_data()
{
  int i, j, m, type;
  double delx, dely, delz, r;
  double **x = atom->x;
  int **bond_type = atom->bond_type;
  double *vfrac = atom->dvector[index_vfrac];
  double *vinter = atom->dvector[index_vinter];

  // vinter[i] = total reference interaction volume (sum of partner nodal volumes),
  // captured once here at the reference configuration for the damage diagnostic.
  // Under newton bond off every bond is stored at both endpoints, so each owned
  // atom sees all its bonds locally.
  for (i = 0; i < atom->nlocal; i++) vinter[i] = 0.0;

  for (i = 0; i < atom->nlocal; i++) {
    for (m = 0; m < atom->num_bond[i]; m++) {
      type = bond_type[i][m];

      // skip turned-off bonds
      if (type <= 0) continue;

      j = atom->map(atom->bond_atom[i][m]);
      if (j == -1) error->one(FLERR, Error::NOLASTLINE, "Atom missing in BPM bond");

      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];

      // closest image in case bonded with a ghost
      domain->minimum_image(FLERR, delx, dely, delz);
      r = sqrt(delx * delx + dely * dely + delz * delz);

      fix_bond_history->update_atom_value(i, m, 0, r);

      // zero any evolving history slots (e.g. the VES deviator extensions)
      for (int a = 1; a < nhistory; a++) fix_bond_history->update_atom_value(i, m, a, 0.0);

      vinter[i] += vfrac[j];
    }
  }
}

/* ----------------------------------------------------------------------
   PMB and LPS bond forces with per-type-pair (#984) breaking. State-based
   models (LPS) run a three-pass scheme inside this single compute(): build the
   weighted volume (once), accumulate and publish the dilatation theta, then the
   force loop reads theta/wvolume of both endpoints of each bond.
------------------------------------------------------------------------- */

void BondBPMPeri::compute(int eflag, int vflag)
{
  // rearrange stored bond values for hybrid bond styles, handled by parent
  pre_compute();

  // grow the per-step break scratch (and the state-based per-atom arrays) if needed
  if (atom->nmax > nmax) {
    memory->destroy(smin_new);
    memory->destroy(s0_new);
    nmax = atom->nmax;
    memory->create(smin_new, nmax, "bond:smin_new");
    memory->create(s0_new, nmax, "bond:s0_new");
    if (state_based) {
      memory->destroy(wvolume);
      memory->destroy(theta);
      memory->destroy(winv);
      memory->create(wvolume, nmax, "bond:wvolume");
      memory->create(theta, nmax, "bond:theta");
      memory->create(winv, nmax, "bond:winv");
      wvolume_setup = 0;    // arrays reallocated -> rebuild the weighted volume
    }
    if (plastic) {
      memory->destroy(tdnorm);
      memory->create(tdnorm, nmax, "bond:tdnorm");
    }
  }

  int i1, i2, itmp, n, type;
  double delx, dely, delz, rsq, r, r0, dr, stretch, fbond, ebond;
  double vfrac_scale, vfrac_eff;

  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  double **bondstore = fix_bond_history->bondstore;
  double *vfrac = atom->dvector[index_vfrac];
  double *smin = atom->dvector[index_smin];

  const double lc = domain->lattice->xlattice;
  const double half_lc = 0.5 * lc;
  const bool allow_breaks = (update->setupflag == 0) && break_flag;

  // sync ghost smin to last step's committed values so the break criterion can
  // read a neighbor's previous-step smin even when the neighbor is a ghost
  commflag = COMM_SMIN;
  comm->forward_comm(this);

  // state-based (LPS) prep: weighted volume (static) then this step's dilatation,
  // each published to ghosts so the force loop can read both bond endpoints
  if (state_based) {
    if (!wvolume_setup)
      compute_wvolume();
    else if (neighbor->ago == 0) {
      commflag = COMM_WVOLUME;
      comm->forward_comm(this);
    }
    compute_dilatation();

    // precompute 1/wvolume once per atom (incl. ghosts) so the per-bond force
    // and plastic-state loops use multiplies instead of repeated divisions
    int ntotal = nlocal + atom->nghost;
    for (int i = 0; i < ntotal; i++) winv[i] = (wvolume[i] > 0.0) ? 1.0 / wvolume[i] : 0.0;

    // volumetric part of the LPS energy (per atom); the deviatoric part is
    // tallied per bond in the force loop below
    if (eflag) {
      for (int i = 0; i < nlocal; i++) {
        double evol = 0.5 * kbulk_rep * theta[i] * theta[i];
        if (eflag_global) energy += evol;
        if (eflag_atom) eatom[i] += evol;
      }
    }

    // EPS plasticity: per-atom deviatoric force-state norm, yield check, and the
    // plastic-multiplier update, published to ghosts for the force loop
    if (plastic) compute_plastic_state();
  }

  // reset the per-step accumulators (MIN for smin, MAX for s0)
  for (int i = 0; i < nlocal; i++) {
    smin_new[i] = DBL_MAX;
    s0_new[i] = -DBL_MAX;
  }

  for (n = 0; n < nbondlist; n++) {

    // skip bond if already broken
    if (bondlist[n][2] <= 0) continue;

    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    // always process from the lower-tagged atom so a bond straddling an MPI
    // boundary (newton bond off) is handled identically on both procs
    if (tag[i2] < tag[i1]) {
      itmp = i1;
      i1 = i2;
      i2 = itmp;
    }

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];
    rsq = delx * delx + dely * dely + delz * delz;
    r = sqrt(rsq);

    // capture the reference length on the first visit
    r0 = bondstore[n][0];
    if ((r0 < EPSILON) || std::isnan(r0)) {
      r0 = bondstore[n][0] = r;
      process_new(n, i1, i2);
    }

    dr = r - r0;
    if (fabs(dr) < NEAR_ZERO) dr = 0.0;    // avoid roundoff noise
    stretch = dr / r0;

    // per-type-pair critical-stretch break test (#984): the critical stretch is
    // s00 - alpha*smin, evaluated per bond, using the more tensile (max) of the
    // two endpoints' previous-step smin. A broken bond exerts no force and does
    // not contribute to the new smin/s0.
    if (allow_breaks && (stretch > s00[type] - alpha[type] * MAX(smin[i1], smin[i2]))) {
      bondlist[n][2] = 0;
      process_broken(i1, i2);
      continue;
    }

    // partial-volume correction: taper the nodal volume for bonds whose
    // reference length is within half a lattice spacing of the horizon
    vfrac_scale = vfrac_taper(r0, cut[type], half_lc);

    // bond force and energy. The mean nodal volume keeps the bond Newton-third-
    // law balanced (equal and opposite); for uniform nodal volume this reduces to
    // legacy's vfrac[j] convention (which conserves momentum only when uniform).
    vfrac_eff = 0.5 * (vfrac[i1] + vfrac[i2]);
    if (model[type] == LPS) {
      const double winv1 = winv[i1], winv2 = winv[i2];
      const double rinv0 = 1.0 / r0;    // isotropic influence function omega = 1/r0
      double rk = 0.0;
      ebond = 0.0;
      if ((winv1 > 0.0) && (winv2 > 0.0)) {
        rk = (3.0 * kbulk[type] - 5.0 * gshear[type]) * vfrac_eff * vfrac_scale *
            (theta[i1] * winv1 + theta[i2] * winv2);
        rk += 15.0 * gshear[type] * vfrac_eff * vfrac_scale * rinv0 *
            (winv1 + winv2) * dr;
        // deviatoric energy per bond (vanishes for a pure dilatation)
        const double dev1 = dr - theta[i1] * r0 / 3.0;
        const double dev2 = dr - theta[i2] * r0 / 3.0;
        ebond = 0.25 * 15.0 * gshear[type] * vfrac_eff * vfrac_scale * rinv0 *
            (dev1 * dev1 * winv1 + dev2 * dev2 * winv2);
      }
      fbond = (r > 0.0) ? -(rk / r) : 0.0;
    } else if (model[type] == VES) {
      // viscoelastic: elastic dilatation (3K*theta) plus a relaxing deviatoric
      // term with a per-bond back-extension recurrence (one symmetric history)
      const double winv1 = winv[i1], winv2 = winv[i2];
      const double rinv0 = 1.0 / r0;
      double rk = 0.0;
      ebond = 0.0;
      if ((winv1 > 0.0) && (winv2 > 0.0)) {
        rk = 3.0 * kbulk[type] * vfrac_eff * vfrac_scale * (theta[i1] * winv1 + theta[i2] * winv2);

        const double dev = dr - 0.5 * (theta[i1] + theta[i2]) * r0 / 3.0;
        const double c1 = tau[type] / update->dt;
        const double decay = exp(-1.0 / c1);
        const double beta = 1.0 - c1 * (1.0 - decay);
        const double dev_prior = bondstore[n][1];
        const double devback_prior = bondstore[n][2];
        const double edb =
            dev_prior * (1.0 - decay) + devback_prior * decay + beta * (dev - dev_prior);

        const double lam = lambda[type];
        const double gfac =
            15.0 * gshear[type] * vfrac_eff * vfrac_scale * rinv0 * (winv1 + winv2);
        rk += gfac * ((1.0 - lam) * dev + lam * (dev - edb));
        ebond = 0.25 * gfac * ((1.0 - lam) * dev * dev + lam * (dev - edb) * (dev - edb));

        // write back the evolving deviator history for the next step
        bondstore[n][1] = dev;
        bondstore[n][2] = edb;
      }
      fbond = (r > 0.0) ? -(rk / r) : 0.0;
    } else if (model[type] == EPS) {
      // elastic-plastic: elastic dilatation (3K*theta) plus a deviatoric force
      // state radially returned onto the yield surface. The per-atom deviatoric
      // force-state norm and yield state are precomputed in compute_plastic_state;
      // the per-bond elastic/plastic split (edp) evolves each step.
      const double winv1 = winv[i1], winv2 = winv[i2];
      const double rinv0 = 1.0 / r0;
      double rk = 0.0;
      ebond = 0.0;
      if ((winv1 > 0.0) && (winv2 > 0.0)) {
        const double dev = dr - 0.5 * (theta[i1] + theta[i2]) * r0 / 3.0;
        const double edp = bondstore[n][1];
        rk = 3.0 * kbulk[type] * (theta[i1] * winv1 + theta[i2] * winv2);

        // trial deviatoric force state; radially return it onto the yield surface
        // when the (symmetric) bond force-state norm exceeds the yield norm
        const double tdtrial = 15.0 * gshear[type] * rinv0 * (winv1 + winv2) * (dev - edp);
        const double tdnorm_b = 0.5 * (tdnorm[i1] + tdnorm[i2]);
        const double yieldnorm = sqrt(2.0 * pointwise_yield);
        double rkdev = tdtrial;
        double edp_new = edp;
        if (tdnorm_b > yieldnorm) {
          const double radial = yieldnorm / tdnorm_b;    // < 1 in the plastic regime
          rkdev = radial * tdtrial;                       // returned deviatoric force state
          // dimensionally consistent plastic-extension increment: the elastic
          // deviatoric extension (dev - edp) is scaled by the same radial factor,
          // edp absorbs the remainder. Converges to the yield surface in one step,
          // unlike the legacy edp += rkNew*deltalambda update whose ~2/r0 gain
          // overshoots (a likely latent PERI bug -- report).
          edp_new = edp + (dev - edp) * (1.0 - radial);
        }
        rk += rkdev;
        // recoverable deviatoric elastic energy uses the returned elastic extension
        const double edev = dev - edp_new;
        ebond = 0.25 * 15.0 * gshear[type] * rinv0 * (winv1 + winv2) * edev * edev *
            vfrac_eff * vfrac_scale;

        // write back the evolving plastic deviator extension for the next step
        bondstore[n][1] = edp_new;
      }
      // elastic + returned-plastic force state, weighted by the nodal volume
      fbond = (r > 0.0) ? -((rk / r) * vfrac_eff * vfrac_scale) : 0.0;
    } else {    // PMB
      fbond = (r > 0.0) ? -(c[type] * vfrac_eff * vfrac_scale * stretch) / r : 0.0;
      ebond = 0.5 * c[type] * vfrac_eff * vfrac_scale * stretch * dr;
    }

    if (newton_bond || (i1 < nlocal)) {
      f[i1][0] += delx * fbond;
      f[i1][1] += dely * fbond;
      f[i1][2] += delz * fbond;
    }
    if (newton_bond || (i2 < nlocal)) {
      f[i2][0] -= delx * fbond;
      f[i2][1] -= dely * fbond;
      f[i2][2] -= delz * fbond;
    }

    // The peridynamic stress integrates over BOTH nodal volumes, so the virial
    // carries one more vfrac factor than the (one-vfrac) bond force and energy.
    // fbond already holds vfrac_eff (matching the force); the extra vfrac_eff
    // here reproduces legacy PERI's fbond*vfrac[i] virial (exact for uniform
    // nodal volume, symmetric otherwise). Energy keeps a single vfrac (ebond).
    if (evflag)
      ev_tally(i1, i2, nlocal, newton_bond, ebond, fbond * vfrac_eff, delx, dely, delz);

    // accumulate this step's minimum stretch (for breaking) and the diagnostic
    // critical stretch on each owned endpoint (newton bond off: each atom sees
    // all of its bonds locally)
    if (i1 < nlocal) {
      smin_new[i1] = MIN(smin_new[i1], stretch);
      s0_new[i1] = MAX(s0_new[i1], s00[type] - alpha[type] * stretch);
    }
    if (i2 < nlocal) {
      smin_new[i2] = MIN(smin_new[i2], stretch);
      s0_new[i2] = MAX(s0_new[i2], s00[type] - alpha[type] * stretch);
    }
  }

  // commit the new smin/s0. An atom with no surviving bonds keeps the no-breaking
  // sentinel (smin = -DBL_MAX) so the max() in the criterion cannot trigger a
  // neighbor's break (its implied critical stretch stays +infinity).
  smin = atom->dvector[index_smin];
  double *s0 = atom->dvector[index_s0];
  for (int i = 0; i < nlocal; i++) {
    smin[i] = (smin_new[i] == DBL_MAX) ? -DBL_MAX : smin_new[i];
    s0[i] = (s0_new[i] == -DBL_MAX) ? DBL_MAX : s0_new[i];
  }

  // revert changes for hybrid bond style, handled by parent
  post_compute();
}

/* ----------------------------------------------------------------------
   weighted volume of each atom from the reference configuration (static).
   wvolume[i] = sum_j omega * |xi|^2 * vfrac[j] * vfrac_scale = sum_j r0*vfrac[j]*scale
   for the isotropic influence function omega = 1/r0. Built once, then
   forward-communicated (refreshed to ghosts on each reneighbor).
------------------------------------------------------------------------- */

void BondBPMPeri::compute_wvolume()
{
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int ntotal = nlocal + atom->nghost;

  double **bondstore = fix_bond_history->bondstore;
  double *vfrac = atom->dvector[index_vfrac];
  const double half_lc = 0.5 * domain->lattice->xlattice;

  for (int i = 0; i < ntotal; i++) wvolume[i] = 0.0;

  for (int n = 0; n < nbondlist; n++) {
    if (bondlist[n][2] <= 0) continue;
    int i1 = bondlist[n][0];
    int i2 = bondlist[n][1];
    int type = bondlist[n][2];

    double r0 = bondstore[n][0];
    double scale = vfrac_taper(r0, cut[type], half_lc);

    if (i1 < nlocal) wvolume[i1] += r0 * vfrac[i2] * scale;
    if (i2 < nlocal) wvolume[i2] += r0 * vfrac[i1] * scale;
  }

  wvolume_setup = 1;
  commflag = COMM_WVOLUME;
  comm->forward_comm(this);
}

/* ----------------------------------------------------------------------
   dilatation theta of each atom (per step), then publish to ghosts.
   theta[i] = (3/wvolume[i]) * sum_j omega*r0*dr*vfrac[j]*scale
            = (3/wvolume[i]) * sum_j dr*vfrac[j]*scale   (omega = 1/r0)
------------------------------------------------------------------------- */

void BondBPMPeri::compute_dilatation()
{
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int ntotal = nlocal + atom->nghost;

  double **x = atom->x;
  double **bondstore = fix_bond_history->bondstore;
  double *vfrac = atom->dvector[index_vfrac];
  const double half_lc = 0.5 * domain->lattice->xlattice;

  for (int i = 0; i < ntotal; i++) theta[i] = 0.0;

  for (int n = 0; n < nbondlist; n++) {
    if (bondlist[n][2] <= 0) continue;
    int i1 = bondlist[n][0];
    int i2 = bondlist[n][1];
    int type = bondlist[n][2];
    // all state-based models (LPS/VES/EPS) carry a dilatation; PMB does not
    if ((model[type] != LPS) && (model[type] != VES) && (model[type] != EPS)) continue;

    double delx = x[i1][0] - x[i2][0];
    double dely = x[i1][1] - x[i2][1];
    double delz = x[i1][2] - x[i2][2];
    double r = sqrt(delx * delx + dely * dely + delz * delz);
    double r0 = bondstore[n][0];
    double dr = r - r0;
    if (fabs(dr) < NEAR_ZERO) dr = 0.0;

    double scale = vfrac_taper(r0, cut[type], half_lc);

    if (i1 < nlocal) theta[i1] += dr * vfrac[i2] * scale;
    if (i2 < nlocal) theta[i2] += dr * vfrac[i1] * scale;
  }

  for (int i = 0; i < nlocal; i++)
    theta[i] = (wvolume[i] > 0.0) ? (3.0 / wvolume[i]) * theta[i] : 0.0;

  // if the user opted in with a d_theta property/atom, publish theta on the
  // owned atoms there so it is visible (e.g. compute property/atom, dump)
  if (index_dtheta >= 0) {
    double *dtheta = atom->dvector[index_dtheta];
    for (int i = 0; i < nlocal; i++) dtheta[i] = theta[i];
  }

  commflag = COMM_THETA;
  comm->forward_comm(this);
}

/* ----------------------------------------------------------------------
   EPS per-atom plastic state: the deviatoric force-state norm tdnorm[i] and the
   yield check, accumulating the plastic multiplier into lambdaValue (a public
   diagnostic). Only tdnorm is published to ghosts; the force loop's radial return
   reads the norm of both endpoints. NOTE: the norm uses the same intensive trial
   force state the force loop applies; legacy PERI's norm carries an extra
   dilatation factor theta (a likely bug), deliberately dropped here.
------------------------------------------------------------------------- */

void BondBPMPeri::compute_plastic_state()
{
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int ntotal = nlocal + atom->nghost;

  double **x = atom->x;
  double **bondstore = fix_bond_history->bondstore;
  double *vfrac = atom->dvector[index_vfrac];
  double *lambdaValue = atom->dvector[index_lambda];
  const double half_lc = 0.5 * domain->lattice->xlattice;

  for (int i = 0; i < ntotal; i++) tdnorm[i] = 0.0;

  // accumulate the squared deviatoric force-state norm over each atom's EPS bonds
  for (int n = 0; n < nbondlist; n++) {
    if (bondlist[n][2] <= 0) continue;
    int i1 = bondlist[n][0];
    int i2 = bondlist[n][1];
    int type = bondlist[n][2];
    if (model[type] != EPS) continue;

    const double winv1 = winv[i1], winv2 = winv[i2];
    if ((winv1 <= 0.0) || (winv2 <= 0.0)) continue;

    double delx = x[i1][0] - x[i2][0];
    double dely = x[i1][1] - x[i2][1];
    double delz = x[i1][2] - x[i2][2];
    double r = sqrt(delx * delx + dely * dely + delz * delz);
    double r0 = bondstore[n][0];
    double dr = r - r0;
    if (fabs(dr) < NEAR_ZERO) dr = 0.0;

    double scale = vfrac_taper(r0, cut[type], half_lc);

    double dev = dr - 0.5 * (theta[i1] + theta[i2]) * r0 / 3.0;
    double edp = bondstore[n][1];
    double tdtrial = 15.0 * gshear[type] * (1.0 / r0) * (winv1 + winv2) * (dev - edp);

    if (i1 < nlocal) tdnorm[i1] += tdtrial * tdtrial * vfrac[i2] * scale;
    if (i2 < nlocal) tdnorm[i2] += tdtrial * tdtrial * vfrac[i1] * scale;
  }

  // per-atom yield-surface check; accumulate the plastic multiplier (lambdaValue,
  // a public diagnostic) and finalize the deviatoric force-state norm in place
  for (int i = 0; i < nlocal; i++) {
    double norm = sqrt(tdnorm[i]);
    tdnorm[i] = norm;
    if (0.5 * norm * norm - pointwise_yield > 0.0)
      lambdaValue[i] +=
          (norm / sqrt(2.0 * pointwise_yield) - 1.0) * wvolume[i] / (15.0 * gshear_rep);
  }

  // publish the finalized norm to ghosts so the force loop can read both endpoints
  commflag = COMM_TDNORM;
  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

void BondBPMPeri::allocate()
{
  allocated = 1;
  const int np1 = atom->nbondtypes + 1;

  memory->create(model, np1, "bond:model");
  memory->create(c, np1, "bond:c");
  memory->create(kbulk, np1, "bond:kbulk");
  memory->create(gshear, np1, "bond:gshear");
  memory->create(lambda, np1, "bond:lambda");
  memory->create(tau, np1, "bond:tau");
  memory->create(yieldstress, np1, "bond:yieldstress");
  memory->create(cut, np1, "bond:cut");
  memory->create(s00, np1, "bond:s00");
  memory->create(alpha, np1, "bond:alpha");

  memory->create(setflag, np1, "bond:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more bond types
   bond_coeff <type> <model> <model-specific coeffs...>
------------------------------------------------------------------------- */

void BondBPMPeri::coeff(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR, "Incorrect args for bond coefficients" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  int model_one;
  double c_one = 0.0, kbulk_one = 0.0, gshear_one = 0.0, lambda_one = 0.0, tau_one = 0.0;
  double yieldstress_one = 0.0, cut_one = 0.0, s00_one = 0.0, alpha_one = 0.0;

  if (strcmp(arg[1], "pmb") == 0) {
    // pmb <c> <horizon> <s00> <alpha>
    if (narg != 6)
      error->all(FLERR, 1,
                 "Incorrect args for bond_style bpm/peri pmb model (expected: <type> pmb c horizon "
                 "s00 alpha)");
    model_one = PMB;
    c_one = utils::numeric(FLERR, arg[2], false, lmp);
    cut_one = utils::numeric(FLERR, arg[3], false, lmp);
    s00_one = utils::numeric(FLERR, arg[4], false, lmp);
    alpha_one = utils::numeric(FLERR, arg[5], false, lmp);
  } else if (strcmp(arg[1], "lps") == 0) {
    // lps <K> <G> <horizon> <s00> <alpha>
    if (narg != 7)
      error->all(FLERR, 1,
                 "Incorrect args for bond_style bpm/peri lps model (expected: <type> lps K G horizon "
                 "s00 alpha)");
    model_one = LPS;
    kbulk_one = utils::numeric(FLERR, arg[2], false, lmp);
    gshear_one = utils::numeric(FLERR, arg[3], false, lmp);
    cut_one = utils::numeric(FLERR, arg[4], false, lmp);
    s00_one = utils::numeric(FLERR, arg[5], false, lmp);
    alpha_one = utils::numeric(FLERR, arg[6], false, lmp);
  } else if (strcmp(arg[1], "ves") == 0) {
    // ves <K> <G> <lambda> <tau> <horizon> <s00> <alpha>
    if (narg != 9)
      error->all(FLERR, 1,
                 "Incorrect args for bond_style bpm/peri ves model (expected: <type> ves K G lambda "
                 "tau horizon s00 alpha)");
    model_one = VES;
    kbulk_one = utils::numeric(FLERR, arg[2], false, lmp);
    gshear_one = utils::numeric(FLERR, arg[3], false, lmp);
    lambda_one = utils::numeric(FLERR, arg[4], false, lmp);
    tau_one = utils::numeric(FLERR, arg[5], false, lmp);
    cut_one = utils::numeric(FLERR, arg[6], false, lmp);
    s00_one = utils::numeric(FLERR, arg[7], false, lmp);
    alpha_one = utils::numeric(FLERR, arg[8], false, lmp);
    if (tau_one <= 0.0)
      error->all(FLERR, 5, "Invalid relaxation time {} for bond_style bpm/peri ves", tau_one);
  } else if (strcmp(arg[1], "eps") == 0) {
    // eps <K> <G> <yieldstress> <horizon> <s00> <alpha>
    if (narg != 8)
      error->all(FLERR, 1,
                 "Incorrect args for bond_style bpm/peri eps model (expected: <type> eps K G "
                 "yieldstress horizon s00 alpha)");
    model_one = EPS;
    kbulk_one = utils::numeric(FLERR, arg[2], false, lmp);
    gshear_one = utils::numeric(FLERR, arg[3], false, lmp);
    yieldstress_one = utils::numeric(FLERR, arg[4], false, lmp);
    cut_one = utils::numeric(FLERR, arg[5], false, lmp);
    s00_one = utils::numeric(FLERR, arg[6], false, lmp);
    alpha_one = utils::numeric(FLERR, arg[7], false, lmp);
    if (yieldstress_one <= 0.0)
      error->all(FLERR, 4, "Invalid yield stress {} for bond_style bpm/peri eps", yieldstress_one);
  } else {
    error->all(FLERR, 1, "Unknown bond_style bpm/peri model: {}", arg[1]);
  }
  if (cut_one <= 0.0) error->all(FLERR, "Invalid horizon value {} for bond_style bpm/peri", cut_one);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    model[i] = model_one;
    c[i] = c_one;
    kbulk[i] = kbulk_one;
    gshear[i] = gshear_one;
    lambda[i] = lambda_one;
    tau[i] = tau_one;
    yieldstress[i] = yieldstress_one;
    cut[i] = cut_one;
    s00[i] = s00_one;
    alpha[i] = alpha_one;
    setflag[i] = 1;
    count++;

    // rough stretch bound for the bond communication-cutoff estimate
    if (1.0 + s00[i] > max_stretch) max_stretch = 1.0 + s00[i];
  }

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients" + utils::errorurl(21));
}

/* ---------------------------------------------------------------------- */

void BondBPMPeri::settings(int narg, char **arg)
{
  BondBPM::settings(narg, arg);

  // bpm/peri has no extra style keywords yet; reject anything left over
  for (auto iarg : leftover_iarg)
    error->all(FLERR, iarg, "Illegal bond_style bpm/peri argument: {}", arg[iarg]);
}

/* ----------------------------------------------------------------------
   check settings, create internal per-atom storage, handshake with pair
------------------------------------------------------------------------- */

void BondBPMPeri::init_style()
{
  // resolve the bond-history size and evolution flag from the active models
  // before the base class creates the bond-history fix; VES stores two evolving
  // deviator values per bond in addition to the reference length r0
  nhistory = 1;
  update_flag = 0;
  for (int i = 1; i <= atom->nbondtypes; i++) {
    if (!setflag[i]) continue;
    if (model[i] == VES) {
      nhistory = MAX(nhistory, 3);
      update_flag = 1;
    } else if (model[i] == EPS) {
      nhistory = MAX(nhistory, 2);
      update_flag = 1;
    }
  }

  BondBPM::init_style();

  // the partial-volume correction needs a uniform cubic lattice spacing
  if (!domain->lattice)
    error->all(FLERR, Error::NOLASTLINE, "Bond style bpm/peri requires a lattice be defined");
  if ((domain->lattice->xlattice != domain->lattice->ylattice) ||
      (domain->lattice->xlattice != domain->lattice->zlattice))
    error->all(FLERR, Error::NOLASTLINE,
               "Bond style bpm/peri requires equal lattice spacing in x, y, and z");

  // detect a state-based model (LPS/VES/EPS) and pick a representative bulk
  // modulus for the per-atom volumetric energy; the weighted volume is rebuilt.
  // For EPS, also precompute the pointwise yield value and a representative shear
  // modulus used by the per-atom return-mapping.
  state_based = 0;
  plastic = 0;
  kbulk_rep = 0.0;
  pointwise_yield = 0.0;
  gshear_rep = 0.0;
  wvolume_setup = 0;
  for (int i = 1; i <= atom->nbondtypes; i++) {
    if (!setflag[i]) continue;
    if ((model[i] == LPS) || (model[i] == VES) || (model[i] == EPS)) {
      state_based = 1;
      kbulk_rep = kbulk[i];
    }
    if (model[i] == EPS) {
      plastic = 1;
      gshear_rep = gshear[i];
      pointwise_yield =
          25.0 / (8.0 * MathConst::MY_PI * pow(cut[i], 5.0)) * yieldstress[i] * yieldstress[i];
    }
  }

  // the plastic (eps) model derives from the Foster/Silling/Chen viscoplasticity
  // theory; surface that citation at runtime when an eps bond type is in use
  if (plastic && lmp->citeme) lmp->citeme->add(cite_bpm_peri_plastic);

  // peridynamics needs a per-atom nodal volume (vfrac); the user supplies it
  // via fix property/atom d_vfrac (uniform with set, or per-atom with read_data)
  int flag, cols;
  index_vfrac = atom->find_custom("vfrac", flag, cols);
  if ((index_vfrac < 0) || (flag != 1) || (cols != 0))
    error->all(FLERR, Error::NOLASTLINE,
               "Bond style bpm/peri requires a per-atom vfrac property; add "
               "'fix <ID> all property/atom d_vfrac ghost yes' before bond_style");

  // internal per-atom bookkeeping: s0 (diagnostic) and smin (break state), the
  // reference interaction volume vinter (for the damage diagnostic, all models),
  // and the EPS accumulated plastic multiplier lambdaValue (d_lambda). Auto-create
  // each field that does not already exist, so a user may pre-declare any of them
  // -- e.g. d_lambda for dumping -- and the bond style reuses it. On restart the
  // property/atom fix is restored before init_style, so the stored s0/smin/vinter
  // values are kept and only freshly created fields get the sentinels.
  bool need_s0 = (atom->find_custom("s0", flag, cols) < 0);
  bool need_smin = (atom->find_custom("smin", flag, cols) < 0);
  bool need_vinter = (atom->find_custom("vinter", flag, cols) < 0);
  bool need_lambda = (plastic && (atom->find_custom("lambda", flag, cols) < 0));
  if ((need_s0 || need_smin || need_vinter || need_lambda) && !id_fix_property_peri) {
    std::string fields;
    if (need_s0) fields += " d_s0";
    if (need_smin) fields += " d_smin";
    if (need_vinter) fields += " d_vinter";
    if (need_lambda) fields += " d_lambda";
    id_fix_property_peri = utils::strdup("BPM_PERI_PROPERTY_ATOM");
    modify->add_fix(fmt::format("{} all property/atom{} ghost yes writedata no",
                                id_fix_property_peri, fields));
  }

  index_s0 = atom->find_custom("s0", flag, cols);
  index_smin = atom->find_custom("smin", flag, cols);
  index_vinter = atom->find_custom("vinter", flag, cols);
  if ((index_s0 < 0) || (index_smin < 0) || (index_vinter < 0))
    error->all(FLERR, Error::NOLASTLINE,
               "Bond style bpm/peri internal error: missing d_s0/d_smin/d_vinter storage");
  if (plastic) {
    index_lambda = atom->find_custom("lambda", flag, cols);
    if (index_lambda < 0)
      error->all(FLERR, Error::NOLASTLINE,
                 "Bond style bpm/peri internal error: missing d_lambda storage");
  }
  // optional dilatation diagnostic: if the user pre-declares a per-atom theta
  // property ('fix <ID> all property/atom d_theta ghost yes' before bond_style),
  // the state-based models publish the per-step dilatation into it (only then,
  // so a compute property/atom defined in the input can read it from the start)
  index_dtheta = (state_based) ? atom->find_custom("theta", flag, cols) : -1;

  // initialize the no-breaking sentinels on the freshly created s0/smin (s0 =
  // +DBL_MAX and smin = -DBL_MAX, so the implied critical stretch is +infinity)
  if (need_s0) {
    double *s0 = atom->dvector[index_s0];
    for (int i = 0; i < atom->nlocal; i++) s0[i] = DBL_MAX;
  }
  if (need_smin) {
    double *smin = atom->dvector[index_smin];
    for (int i = 0; i < atom->nlocal; i++) smin[i] = -DBL_MAX;
  }

  // peridynamics needs the short-range contact pair; require pair_style bpm/peri.
  // A deliberate pair_style zero is also accepted so the bond style can be driven
  // in isolation by the test harness (and the explicit no-contact regime). This
  // blocks accidental misuse without forbidding a deliberate opt-out.
  if (!force->pair || (!force->pair_match("bpm/peri", 0) && !force->pair_match("zero", 1)))
    error->all(FLERR, Error::NOLASTLINE, "Bond style bpm/peri requires pair style bpm/peri");

  // these models do a single force pass and are several times faster than the
  // legacy peri/* styles, so a multi-timescale rRESPA split is never worthwhile;
  // refuse it rather than silently mis-integrate (no respa hooks are implemented)
  if (utils::strmatch(update->integrate_style, "^respa"))
    error->all(FLERR, Error::NOLASTLINE,
               "Bond style bpm/peri is not compatible with run_style respa");

  // consistency sanity check (typos happen): the bond horizon and the contact
  // pair horizon should be the same peridynamic horizon.  Compare against the
  // pair's largest cutoff when the bpm/peri contact pair is in use and warn (do
  // not error) on an obvious mismatch -- differing by more than half a lattice
  // spacing, the width of the partial-volume taper.
  if (force->pair_match("bpm/peri", 0)) {
    const double half_lc = 0.5 * domain->lattice->xlattice;
    for (int i = 1; i <= atom->nbondtypes; i++) {
      if (!setflag[i]) continue;
      if (fabs(force->pair->cutforce - cut[i]) > half_lc)
        error->warning(FLERR,
                       "Bond style bpm/peri horizon {} for bond type {} differs from the pair "
                       "style bpm/peri contact cutoff {}; check the pair_coeff and bond_coeff",
                       cut[i], i, force->pair->cutforce);
    }
  }
}

/* ---------------------------------------------------------------------- */

double BondBPMPeri::single(int type, double rsq, int i, int j, double &fforce)
{
  if (type <= 0) {
    fforce = 0.0;
    return 0.0;
  }

  double r0 = 0.0;
  for (int n = 0; n < atom->num_bond[i]; n++)
    if (atom->bond_atom[i][n] == atom->tag[j]) r0 = fix_bond_history->get_atom_value(i, n, 0);

  double r = sqrt(rsq);
  double dr = r - r0;
  if (fabs(dr) < NEAR_ZERO) dr = 0.0;
  double stretch = (r0 > 0.0) ? dr / r0 : 0.0;

  const double half_lc = 0.5 * domain->lattice->xlattice;
  double vfrac_scale = vfrac_taper(r0, cut[type], half_lc);

  double *vfrac = atom->dvector[index_vfrac];
  double vfrac_eff = 0.5 * (vfrac[i] + vfrac[j]);

  fforce = (r > 0.0) ? -(c[type] * vfrac_eff * vfrac_scale * stretch) / r : 0.0;
  svector[0] = r0;
  return 0.5 * c[type] * vfrac_eff * vfrac_scale * stretch * dr;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondBPMPeri::write_restart(FILE *fp)
{
  BondBPM::write_restart(fp);

  fwrite(&model[1], sizeof(int), atom->nbondtypes, fp);
  fwrite(&c[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&kbulk[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&gshear[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&lambda[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&tau[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&yieldstress[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&cut[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&s00[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&alpha[1], sizeof(double), atom->nbondtypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondBPMPeri::read_restart(FILE *fp)
{
  BondBPM::read_restart(fp);
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &model[1], sizeof(int), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &c[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &kbulk[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &gshear[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &lambda[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &tau[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &yieldstress[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &cut[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &s00[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &alpha[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
  }
  MPI_Bcast(&model[1], atom->nbondtypes, MPI_INT, 0, world);
  MPI_Bcast(&c[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&kbulk[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gshear[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&lambda[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&tau[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&yieldstress[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&s00[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&alpha[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   forward communication of one per-atom array to ghosts, selected by commflag:
   smin (break test), theta (LPS dilatation), or wvolume (LPS weighted volume)
------------------------------------------------------------------------- */

int BondBPMPeri::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  double *src;
  if (commflag == COMM_THETA) src = theta;
  else if (commflag == COMM_WVOLUME) src = wvolume;
  else if (commflag == COMM_TDNORM) src = tdnorm;
  else src = atom->dvector[index_smin];

  int m = 0;
  for (int i = 0; i < n; i++) buf[m++] = src[list[i]];
  return m;
}

/* ---------------------------------------------------------------------- */

void BondBPMPeri::unpack_forward_comm(int n, int first, double *buf)
{
  double *dst;
  if (commflag == COMM_THETA) dst = theta;
  else if (commflag == COMM_WVOLUME) dst = wvolume;
  else if (commflag == COMM_TDNORM) dst = tdnorm;
  else dst = atom->dvector[index_smin];

  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) dst[i] = buf[m++];
}
