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
   Re-implementation of the PERI package constitutive models on the BPM
   framework. This is the Stage 0 scaffold: it builds the peridynamic bond
   set, stores the reference configuration, and applies zero force. The PMB
   force law, the state-based (LPS/VES/EPS) models, and the #984 per-bond
   breaking criterion are added in subsequent stages.
------------------------------------------------------------------------- */

#include "bond_bpm_peri.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_bond_history.h"
#include "force.h"
#include "lattice.h"
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

/* ---------------------------------------------------------------------- */

BondBPMPeri::BondBPMPeri(LAMMPS *_lmp) :
    BondBPM(_lmp), model(nullptr), c(nullptr), kbulk(nullptr), gshear(nullptr), cut(nullptr),
    s00(nullptr), alpha(nullptr), id_fix_property_peri(nullptr), index_vfrac(-1), index_s0(-1),
    index_smin(-1), smin_new(nullptr), s0_new(nullptr), nmax(0), state_based(0), wvolume_setup(0),
    kbulk_rep(0.0), wvolume(nullptr), theta(nullptr), commflag(COMM_SMIN)
{
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

  if (id_fix_property_peri && modify->nfix) {
    modify->delete_fix(id_fix_property_peri);
    delete[] id_fix_property_peri;
  }

  if (allocated) {
    memory->destroy(model);
    memory->destroy(c);
    memory->destroy(kbulk);
    memory->destroy(gshear);
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
      memory->create(wvolume, nmax, "bond:wvolume");
      memory->create(theta, nmax, "bond:theta");
      wvolume_setup = 0;    // arrays reallocated -> rebuild the weighted volume
    }
  }

  int i1, i2, itmp, n, type;
  double delx, dely, delz, rsq, r, r0, dr, stretch, fbond, ebond;
  double delta, vfrac_scale, vfrac_eff;

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

    // volumetric part of the LPS energy (per atom); the deviatoric part is
    // tallied per bond in the force loop below
    if (eflag) {
      for (int i = 0; i < nlocal; i++) {
        double evol = 0.5 * kbulk_rep * theta[i] * theta[i];
        if (eflag_global) energy += evol;
        if (eflag_atom) eatom[i] += evol;
      }
    }
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
    if (r0 < EPSILON || std::isnan(r0)) {
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
    if (allow_breaks && stretch > s00[type] - alpha[type] * MAX(smin[i1], smin[i2])) {
      bondlist[n][2] = 0;
      process_broken(i1, i2);
      continue;
    }

    // partial-volume correction: taper the nodal volume for bonds whose
    // reference length is within half a lattice spacing of the horizon
    delta = cut[type];
    if (fabs(r0 - delta) <= half_lc)
      vfrac_scale = (-1.0 / (2.0 * half_lc)) * r0 + (1.0 + (delta - half_lc) / (2.0 * half_lc));
    else
      vfrac_scale = 1.0;

    // bond force and energy. The mean nodal volume keeps the bond Newton-third-
    // law balanced (equal and opposite); for uniform nodal volume this reduces to
    // legacy's vfrac[j] convention (which conserves momentum only when uniform).
    vfrac_eff = 0.5 * (vfrac[i1] + vfrac[i2]);
    if (model[type] == LPS) {
      const double wv1 = wvolume[i1], wv2 = wvolume[i2];
      const double rinv0 = 1.0 / r0;    // isotropic influence function omega = 1/r0
      double rk = 0.0;
      ebond = 0.0;
      if (wv1 > 0.0 && wv2 > 0.0) {
        rk = (3.0 * kbulk[type] - 5.0 * gshear[type]) * vfrac_eff * vfrac_scale *
            (theta[i1] / wv1 + theta[i2] / wv2);
        rk += 15.0 * gshear[type] * vfrac_eff * vfrac_scale * rinv0 *
            (1.0 / wv1 + 1.0 / wv2) * dr;
        // deviatoric energy per bond (vanishes for a pure dilatation)
        const double dev1 = dr - theta[i1] * r0 / 3.0;
        const double dev2 = dr - theta[i2] * r0 / 3.0;
        ebond = 0.25 * 15.0 * gshear[type] * vfrac_eff * vfrac_scale * rinv0 *
            (dev1 * dev1 / wv1 + dev2 * dev2 / wv2);
      }
      fbond = (r > 0.0) ? -(rk / r) : 0.0;
    } else {    // PMB
      fbond = (r > 0.0) ? -(c[type] * vfrac_eff * vfrac_scale * stretch) / r : 0.0;
      ebond = 0.5 * c[type] * vfrac_eff * vfrac_scale * stretch * dr;
    }

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += delx * fbond;
      f[i1][1] += dely * fbond;
      f[i1][2] += delz * fbond;
    }
    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= delx * fbond;
      f[i2][1] -= dely * fbond;
      f[i2][2] -= delz * fbond;
    }

    if (evflag) ev_tally(i1, i2, nlocal, newton_bond, ebond, fbond, delx, dely, delz);

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
    double delta = cut[type];
    double scale = 1.0;
    if (fabs(r0 - delta) <= half_lc)
      scale = (-1.0 / (2.0 * half_lc)) * r0 + (1.0 + (delta - half_lc) / (2.0 * half_lc));

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
    if (model[type] != LPS) continue;

    double delx = x[i1][0] - x[i2][0];
    double dely = x[i1][1] - x[i2][1];
    double delz = x[i1][2] - x[i2][2];
    double r = sqrt(delx * delx + dely * dely + delz * delz);
    double r0 = bondstore[n][0];
    double dr = r - r0;
    if (fabs(dr) < NEAR_ZERO) dr = 0.0;

    double delta = cut[type];
    double scale = 1.0;
    if (fabs(r0 - delta) <= half_lc)
      scale = (-1.0 / (2.0 * half_lc)) * r0 + (1.0 + (delta - half_lc) / (2.0 * half_lc));

    if (i1 < nlocal) theta[i1] += dr * vfrac[i2] * scale;
    if (i2 < nlocal) theta[i2] += dr * vfrac[i1] * scale;
  }

  for (int i = 0; i < nlocal; i++)
    theta[i] = (wvolume[i] > 0.0) ? (3.0 / wvolume[i]) * theta[i] : 0.0;

  commflag = COMM_THETA;
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
  double c_one = 0.0, kbulk_one = 0.0, gshear_one = 0.0;
  double cut_one = 0.0, s00_one = 0.0, alpha_one = 0.0;

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
  BondBPM::init_style();

  // the partial-volume correction needs a uniform cubic lattice spacing
  if (!domain->lattice)
    error->all(FLERR, Error::NOLASTLINE, "Bond style bpm/peri requires a lattice be defined");
  if (domain->lattice->xlattice != domain->lattice->ylattice ||
      domain->lattice->xlattice != domain->lattice->zlattice)
    error->all(FLERR, Error::NOLASTLINE,
               "Bond style bpm/peri requires equal lattice spacing in x, y, and z");

  // detect a state-based model (LPS) and pick a representative bulk modulus for
  // the per-atom volumetric energy; the weighted volume is rebuilt on next use
  state_based = 0;
  kbulk_rep = 0.0;
  wvolume_setup = 0;
  for (int i = 1; i <= atom->nbondtypes; i++) {
    if (setflag[i] && model[i] == LPS) {
      state_based = 1;
      kbulk_rep = kbulk[i];
    }
  }

  // peridynamics needs a per-atom nodal volume (vfrac); the user supplies it
  // via fix property/atom d_vfrac (uniform with set, or per-atom with read_data)
  int flag, cols;
  index_vfrac = atom->find_custom("vfrac", flag, cols);
  if (index_vfrac < 0 || flag != 1 || cols != 0)
    error->all(FLERR, Error::NOLASTLINE,
               "Bond style bpm/peri requires a per-atom vfrac property; add "
               "'fix <ID> all property/atom d_vfrac ghost yes' before bond_style");

  // internal critical-stretch bookkeeping: s0 (diagnostic) and smin (break
  // state); auto-create the fix property/atom if it does not already exist
  bool created = false;
  index_s0 = atom->find_custom("s0", flag, cols);
  if (index_s0 < 0) {
    if (!id_fix_property_peri) {
      id_fix_property_peri = utils::strdup("BPM_PERI_PROPERTY_ATOM");
      modify->add_fix(fmt::format("{} all property/atom d_s0 d_smin ghost yes writedata no",
                                  id_fix_property_peri));
      created = true;
    }
    index_s0 = atom->find_custom("s0", flag, cols);
  }
  index_smin = atom->find_custom("smin", flag, cols);
  if (index_smin < 0)
    error->all(FLERR, Error::NOLASTLINE,
               "Bond style bpm/peri internal error: missing d_smin storage");

  // initialize the no-breaking sentinels on first creation (s0 = +DBL_MAX and
  // smin = -DBL_MAX, so the implied critical stretch is +infinity); on restart
  // the stored values are restored, so only initialize when freshly created
  if (created) {
    double *s0 = atom->dvector[index_s0];
    double *smin = atom->dvector[index_smin];
    for (int i = 0; i < atom->nlocal; i++) {
      s0[i] = DBL_MAX;
      smin[i] = -DBL_MAX;
    }
  }

  // peridynamics needs the short-range contact pair; require pair_style bpm/peri.
  // A deliberate pair_style zero is also accepted so the bond style can be driven
  // in isolation by the test harness (and the explicit no-contact regime). This
  // blocks accidental misuse without forbidding a deliberate opt-out.
  if (!force->pair || (!force->pair_match("bpm/peri", 0) && !force->pair_match("zero", 1)))
    error->all(FLERR, Error::NOLASTLINE, "Bond style bpm/peri requires pair style bpm/peri");
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

  const double lc = domain->lattice->xlattice;
  const double half_lc = 0.5 * lc;
  double delta = cut[type];
  double vfrac_scale = 1.0;
  if (fabs(r0 - delta) <= half_lc)
    vfrac_scale = (-1.0 / (2.0 * half_lc)) * r0 + (1.0 + (delta - half_lc) / (2.0 * half_lc));

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
    utils::sfread(FLERR, &cut[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &s00[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &alpha[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
  }
  MPI_Bcast(&model[1], atom->nbondtypes, MPI_INT, 0, world);
  MPI_Bcast(&c[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&kbulk[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gshear[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
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
  else dst = atom->dvector[index_smin];

  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) dst[i] = buf[m++];
}
