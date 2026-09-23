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
   fix active

   Each model is a distinct dynamical scheme; none simply alias to ABP.

   Usage:
     fix ID group active model keyword value [keyword value ...]

   All model-specific parameters are specified as <keyword value> pairs
   (order-independent).  The list of keywords actually consumed by each
   model is given below; passing an unrelated keyword is harmless but is
   silently ignored.  Defaults are applied when a keyword is omitted.

   Models and recognised keywords:
     abp     v0 Dr [Dt] [gamma] [seed]
     aoup    Da tau [gamma] [seed]
     chiral  v0 Dr omega [gamma] [seed]
     rtp     v0 tumble_rate [Dt] [gamma] [seed]
     vicsek  v0 Dr Rcut [alignment] [gamma] [seed]
     rod     v0 Dr aspect [kT] [gamma] [seed]      (self-contained Langevin)
     spinner omega_spin eta_odd [seed]              (drag from fix langevin)
     iabp    v0 Dr gamma kT [seed]                  (self-contained Langevin)
     mips    v0 Dr Rcut rho_max [gamma] [seed]
     nematic activity K Rcut [Dr] [gamma] [seed]

   Combination notes:
     - rod and iabp supply their own Stokes drag and FDT noise; combine
       only with `fix nve`, NOT with `fix langevin` or `fix brownian`,
       otherwise the dissipation is double-counted.
     - All other models add only the active force; combine with
       `fix nve + fix langevin` to obtain a complete overdamped
       Langevin integrator with the desired k_B T.
     - spinner: the friction gamma that appears in the Hall-angle
       prediction tan(phi) = eta_odd * omega_spin / gamma is the drag
       supplied by fix langevin (gamma = m / Tdamp); the spinner kernel
       itself contributes only the odd-viscosity force eta_odd * (z x v).

   Parallel computing:
     Per-atom orientations are forward-communicated to ghost atoms before
     the Vicsek and nematic kernels iterate the neighbour list, so the
     alignment kernels are correct under arbitrary spatial decomposition.

   Spatial scope:
     The active force is currently applied in the xy-plane.  A 3D box is
     allowed but the active push has no z-component (a warning is emitted
     at init).  Quaternion storage is allocated as scaffolding for a
     future 3D extension but is not used by any kernel in this release.

   Contributing author: Dr. Ram Chand (The Begum Nusrat Bhutto Women University,
                        Sukkur, Sindh, Pakistan)
                        ram.chand@bnbwu.edu.pk, ram.chand2k11@yahoo.com
------------------------------------------------------------------------- */

#include "fix_active.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "random_mars.h"
#include "update.h"
#include "math_const.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixActive::FixActive(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), nmax_alloc(0),
    theta(nullptr), quat(nullptr), va(nullptr),
    random(nullptr), list(nullptr)
{
  if (narg < 4) error->all(FLERR, "Illegal fix active command");

  if (strcmp(arg[3], "abp") == 0) {
    model = ABP;
  } else if (strcmp(arg[3], "aoup") == 0) {
    model = AOUP;
  } else if (strcmp(arg[3], "chiral") == 0) {
    model = CHIRAL;
  } else if (strcmp(arg[3], "rtp") == 0) {
    model = RTP;
  } else if (strcmp(arg[3], "vicsek") == 0) {
    model = VICSEK;
  } else if (strcmp(arg[3], "rod") == 0) {
    model = ROD;
  } else if (strcmp(arg[3], "spinner") == 0) {
    model = SPINNER;
  } else if (strcmp(arg[3], "iabp") == 0) {
    model = IABP;
  } else if (strcmp(arg[3], "mips") == 0) {
    model = MIPS;
  } else if (strcmp(arg[3], "nematic") == 0) {
    model = NEMATIC;
  } else {
    error->all(FLERR, "Unknown active model: {}", arg[3]);
  }

  v0 = 1.0;
  Dr = 0.1;
  Dt = 0.0;
  seed = 12345;
  dimension = domain->dimension;
  omega = 0.0;
  tumble_rate = 0.1;
  Rcut = 1.0;
  alignment = 1.0;
  aspect = 2.0;
  omega_spin = 1.0;
  eta_odd = 0.0;
  gamma = 1.0;
  tau = 1.0;
  Da = 1.0;
  activity = 1.0;
  K = 1.0;
  kT = 1.0;
  rho_max = 1.0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "v0") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      v0 = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "Dr") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      Dr = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "Dt") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      Dt = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "seed") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      seed = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "omega") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      omega = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "tumble_rate") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      tumble_rate = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "Rcut") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      Rcut = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "alignment") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      alignment = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "aspect") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      aspect = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "omega_spin") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      omega_spin = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "eta_odd") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      eta_odd = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "gamma") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      gamma = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "tau") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      tau = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "Da") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      Da = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "activity") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      activity = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "K") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      K = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "kT") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      kT = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "rho_max") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      rho_max = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "dimension") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix active command");
      dimension = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix active command: unknown keyword {}", arg[iarg]);
    }
  }

  random = new RanMars(lmp, seed + comm->me);

  // Forward-communicate one double per atom (theta) to ghosts.  Only the
  // Vicsek and nematic kernels actually trigger this, but the size is
  // declared once here so LAMMPS's communication buffers are sized.
  comm_forward = 1;

  grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    if (dimension == 2) {
      theta[i][0] = 2.0 * MY_PI * random->uniform();
    } else {
      double u1 = random->uniform();
      double u2 = random->uniform();
      double u3 = random->uniform();
      quat[i][0] = sqrt(1.0 - u1) * sin(2.0 * MY_PI * u2);
      quat[i][1] = sqrt(1.0 - u1) * cos(2.0 * MY_PI * u2);
      quat[i][2] = sqrt(u1) * sin(2.0 * MY_PI * u3);
      quat[i][3] = sqrt(u1) * cos(2.0 * MY_PI * u3);
    }
  }

  if (model == AOUP) {
    for (int i = 0; i < nlocal; i++) {
      va[i][0] = sqrt(Da) * random->gaussian();
      va[i][1] = sqrt(Da) * random->gaussian();
      if (dimension == 3) va[i][2] = sqrt(Da) * random->gaussian();
    }
  }
}

/* ---------------------------------------------------------------------- */

FixActive::~FixActive()
{
  atom->delete_callback(id, Atom::GROW);
  memory->destroy(theta);
  memory->destroy(quat);
  memory->destroy(va);
  delete random;
}

/* ----------------------------------------------------------------------
   Allocate / reallocate per-atom storage.  Called from the constructor
   and from the GROW atom callback whenever atom->nmax changes.
------------------------------------------------------------------------- */

void FixActive::grow_arrays(int nmax)
{
  int old = nmax_alloc;
  memory->grow(theta, nmax, 1, "fix_active:theta");
  memory->grow(quat,  nmax, 4, "fix_active:quat");
  memory->grow(va,    nmax, 3, "fix_active:va");

  // Initialise newly allocated slots so atoms created after construction
  // (e.g. via fix deposit) start with well-defined orientations.
  for (int i = old; i < nmax; i++) {
    theta[i][0] = 0.0;
    quat[i][0] = 1.0; quat[i][1] = 0.0; quat[i][2] = 0.0; quat[i][3] = 0.0;
    va[i][0] = va[i][1] = va[i][2] = 0.0;
  }
  nmax_alloc = nmax;
}

/* ----------------------------------------------------------------------
   Per-atom data movement on sort and MPI exchange.  Without these the
   theta/quat/va arrays become misaligned with x/v whenever LAMMPS
   reorders local atoms, which silently scrambles the propulsion direction.
------------------------------------------------------------------------- */

void FixActive::copy_arrays(int i, int j, int /*delflag*/)
{
  theta[j][0] = theta[i][0];
  quat[j][0] = quat[i][0];
  quat[j][1] = quat[i][1];
  quat[j][2] = quat[i][2];
  quat[j][3] = quat[i][3];
  va[j][0] = va[i][0];
  va[j][1] = va[i][1];
  va[j][2] = va[i][2];
}

int FixActive::pack_exchange(int i, double *buf)
{
  int n = 0;
  buf[n++] = theta[i][0];
  buf[n++] = quat[i][0];
  buf[n++] = quat[i][1];
  buf[n++] = quat[i][2];
  buf[n++] = quat[i][3];
  buf[n++] = va[i][0];
  buf[n++] = va[i][1];
  buf[n++] = va[i][2];
  return n;
}

int FixActive::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  theta[nlocal][0] = buf[n++];
  quat[nlocal][0] = buf[n++];
  quat[nlocal][1] = buf[n++];
  quat[nlocal][2] = buf[n++];
  quat[nlocal][3] = buf[n++];
  va[nlocal][0] = buf[n++];
  va[nlocal][1] = buf[n++];
  va[nlocal][2] = buf[n++];
  return n;
}

/* ----------------------------------------------------------------------
   Ghost-atom forward communication of the per-atom orientation.

   Called from compute_vicsek() and compute_nematic() before the kernel
   iterates over the neighbour list, so that theta[j] is correct for j in
   [nlocal, nlocal+nghost).  Without this, neighbours owned by another
   MPI rank would be skipped by the alignment kernels (their theta would
   be uninitialised storage rather than the owner's current value),
   biasing the local mean field along every subdomain boundary.
------------------------------------------------------------------------- */

int FixActive::pack_forward_comm(int n, int *list, double *buf,
                                  int /*pbc_flag*/, int * /*pbc*/)
{
  for (int i = 0; i < n; i++) buf[i] = theta[list[i]][0];
  return n;
}

void FixActive::unpack_forward_comm(int n, int first, double *buf)
{
  for (int i = 0; i < n; i++) theta[first + i][0] = buf[i];
}

/* ---------------------------------------------------------------------- */

int FixActive::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixActive::init()
{
  // The compute kernels only update f[i][0] and f[i][1].  In a 3D box
  // the active force has no z-component; warn the user instead of
  // failing silently.
  if (domain->dimension == 3 && comm->me == 0)
    error->warning(FLERR,
                   "fix active currently propels in the xy-plane only; "
                   "3D simulations will see no active force in z");

  // VICSEK / MIPS / NEMATIC need the full list (all j around i) for
  // mean-field summations.  REQ_DEFAULT (half list) would only see j>i.
  if (model == VICSEK || model == NEMATIC || model == MIPS) {
    neighbor->add_request(this, NeighConst::REQ_FULL);

    // Warn if the requested alignment / density cutoff exceeds the
    // global neighbor cutoff: such neighbors will be missed.
    if (Rcut > neighbor->cutneighmax && comm->me == 0)
      error->warning(FLERR,
                     "fix active Rcut exceeds neighbor cutoff; increase "
                     "the pair cutoff or neighbor skin");
  }
}

/* ---------------------------------------------------------------------- */

void FixActive::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixActive::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixActive::post_force(int /*vflag*/)
{
  switch (model) {
    case ABP:     compute_abp();     break;
    case AOUP:    compute_aoup();    break;
    case CHIRAL:  compute_chiral();  break;
    case RTP:     compute_rtp();     break;
    case VICSEK:  compute_vicsek();  break;
    case ROD:     compute_rod();     break;
    case SPINNER: compute_spinner(); break;
    case IABP:    compute_iabp();    break;
    case MIPS:    compute_mips();    break;
    case NEMATIC: compute_nematic(); break;
  }
}

/* ----------------------------------------------------------------------
   ABP: dr/dt = v0 e_hat + mu F + sqrt(2 Dt) xi,    dtheta/dt = sqrt(2 Dr) xi
------------------------------------------------------------------------- */

void FixActive::compute_abp()
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dt = update->dt;
  double noise_rot = sqrt(2.0 * Dr * dt);

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      theta[i][0] += noise_rot * random->gaussian();

      f[i][0] += gamma * v0 * cos(theta[i][0]);
      f[i][1] += gamma * v0 * sin(theta[i][0]);

      if (Dt > 0.0) {
        double noise_amp = gamma * sqrt(2.0 * Dt / dt);
        f[i][0] += noise_amp * random->gaussian();
        f[i][1] += noise_amp * random->gaussian();
      }
    }
  }
}

/* ----------------------------------------------------------------------
   AOUP: dr/dt = mu F + v_a,    dv_a/dt = -v_a/tau + sqrt(2 Da/tau) eta
   Exact integration of OU process keeps the equilibrium variance Da.
------------------------------------------------------------------------- */

void FixActive::compute_aoup()
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dt = update->dt;
  double decay = exp(-dt / tau);
  double noise_amp = sqrt(Da * (1.0 - decay * decay));

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      va[i][0] = decay * va[i][0] + noise_amp * random->gaussian();
      va[i][1] = decay * va[i][1] + noise_amp * random->gaussian();
      f[i][0] += gamma * va[i][0];
      f[i][1] += gamma * va[i][1];
    }
  }
}

/* ----------------------------------------------------------------------
   Chiral ABP: same as ABP plus deterministic angular drift omega.
------------------------------------------------------------------------- */

void FixActive::compute_chiral()
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dt = update->dt;
  double noise_rot = sqrt(2.0 * Dr * dt);

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      theta[i][0] += omega * dt + noise_rot * random->gaussian();
      f[i][0] += gamma * v0 * cos(theta[i][0]);
      f[i][1] += gamma * v0 * sin(theta[i][0]);
    }
  }
}

/* ----------------------------------------------------------------------
   Run-and-Tumble: ballistic runs of v0 punctuated by Poisson tumbles
   that reset orientation uniformly.  No continuous rotational diffusion.
------------------------------------------------------------------------- */

void FixActive::compute_rtp()
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dt = update->dt;
  double p_tumble = 1.0 - exp(-tumble_rate * dt);

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (random->uniform() < p_tumble)
        theta[i][0] = 2.0 * MY_PI * random->uniform();

      f[i][0] += gamma * v0 * cos(theta[i][0]);
      f[i][1] += gamma * v0 * sin(theta[i][0]);
    }
  }
}

/* ----------------------------------------------------------------------
   Vicsek: each particle aligns its orientation toward the circular mean
   of the orientations of all neighbors within Rcut, plus angular noise.
   The parameter `alignment` in [0,1] interpolates between no alignment
   and full Vicsek alignment by linearly blending the unit orientation
   vectors of self and the neighbour mean (NOT the angles themselves,
   which would behave incorrectly across the +/- pi branch cut).
------------------------------------------------------------------------- */

void FixActive::compute_vicsek()
{
  if (!list) return;

  // Make theta[j] valid for j in [nlocal, nlocal+nghost) so that
  // alignment partners across MPI subdomain boundaries are not dropped.
  comm->forward_comm(this);

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dt = update->dt;
  double Rcut2 = Rcut * Rcut;
  double noise_rot = sqrt(2.0 * Dr * dt);

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // Snapshot current orientations so all updates see the same t-state.
  // We size the buffer to atom->nmax so that ilist entries (which can equal
  // any local atom index) are always in range.
  int nmax = atom->nmax;
  double *new_theta = new double[nmax];
  for (int i = 0; i < nmax; i++) new_theta[i] = theta[i][0];

  // Build a flag of which atoms received an alignment update so the second
  // loop can fall back to plain rotational diffusion for stragglers (atoms
  // whose neighbor list was not iterated this step).
  char *updated = new char[nlocal];
  for (int i = 0; i < nlocal; i++) updated[i] = 0;

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    if (i >= nlocal) continue;

    double sin_sum = sin(theta[i][0]);
    double cos_sum = cos(theta[i][0]);

    int *jlist = firstneigh[i];
    int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj] & NEIGHMASK;
      // j may be a ghost; its theta was set by the forward_comm above.

      double dx = x[i][0] - x[j][0];
      double dy = x[i][1] - x[j][1];
      if (dx * dx + dy * dy < Rcut2) {
        sin_sum += sin(theta[j][0]);
        cos_sum += cos(theta[j][0]);
      }
    }

    // Circular mean of self+neighbours.
    double mean_cos = cos_sum;
    double mean_sin = sin_sum;
    // Blend on the unit circle so the interpolation is well-defined
    // at the +/- pi branch cut.  With alignment=1 we recover the standard
    // Vicsek mean; with alignment=0 we recover the particle's own axis.
    double cos_self = cos(theta[i][0]);
    double sin_self = sin(theta[i][0]);
    double cos_blend = (1.0 - alignment) * cos_self + alignment * mean_cos;
    double sin_blend = (1.0 - alignment) * sin_self + alignment * mean_sin;
    double blend_angle = atan2(sin_blend, cos_blend);

    new_theta[i] = blend_angle + noise_rot * random->gaussian();
    updated[i] = 1;
  }

  // Atoms that did not appear in the iterated list still need rotational
  // diffusion so that their orientation does not freeze.
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if (!updated[i])
      new_theta[i] = theta[i][0] + noise_rot * random->gaussian();
  }

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      theta[i][0] = new_theta[i];
      f[i][0] += gamma * v0 * cos(theta[i][0]);
      f[i][1] += gamma * v0 * sin(theta[i][0]);
    }
  }

  delete[] new_theta;
  delete[] updated;
}

/* ----------------------------------------------------------------------
   Active rod (Baskaran-Marchetti style, self-contained inertial Langevin)

   We integrate the rod's translational dynamics ourselves so the
   anisotropic mobility is consistent: do NOT combine this fix with
   fix langevin or fix brownian - they would re-add isotropic drag/noise
   and double-count the dissipation.

   gamma_par   = gamma / aspect       (parallel friction)
   gamma_perp  = gamma                 (perpendicular friction)

   F_active    = gamma_par v0 e_hat   (so terminal v_par = v0)
   F_drag      = -gamma_par v_par e_hat - gamma_perp v_perp_vec
   F_noise     = sqrt(2 gamma_alpha kT / dt) xi_alpha   (alpha in {par, perp})

   Aligning torque from EXTERNAL forces (assumed to be in atom->f at entry):
     tau_z = (e_hat x F_ext)_z = ex Fy - ey Fx
------------------------------------------------------------------------- */

void FixActive::compute_rod()
{
  double **f = atom->f;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dt = update->dt;
  double noise_rot = sqrt(2.0 * Dr * dt);

  double gamma_par  = gamma / aspect;
  double gamma_perp = gamma;
  double noise_par  = sqrt(2.0 * gamma_par  * kT / dt);
  double noise_perp = sqrt(2.0 * gamma_perp * kT / dt);

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double ex = cos(theta[i][0]);
      double ey = sin(theta[i][0]);

      // External force on entry (e.g., pair forces).
      double Fext_x = f[i][0];
      double Fext_y = f[i][1];

      // Aligning torque from the external force only.
      double tau_z = ex * Fext_y - ey * Fext_x;
      theta[i][0] += (tau_z / gamma_perp) * dt + noise_rot * random->gaussian();

      // Anisotropic Langevin drag and FDT noise.
      double vpar = v[i][0] * ex + v[i][1] * ey;
      double vperp_x = v[i][0] - vpar * ex;
      double vperp_y = v[i][1] - vpar * ey;

      double xi_par  = random->gaussian();
      double xi_perp = random->gaussian();

      // Recompute body axis after rotation update.
      double exu = cos(theta[i][0]);
      double eyu = sin(theta[i][0]);
      double pxu = -eyu, pyu = exu;     // unit perpendicular

      f[i][0] += -gamma_par * vpar * exu - gamma_perp * vperp_x
                 + noise_par * xi_par * exu + noise_perp * xi_perp * pxu
                 + gamma_par * v0 * exu;
      f[i][1] += -gamma_par * vpar * eyu - gamma_perp * vperp_y
                 + noise_par * xi_par * eyu + noise_perp * xi_perp * pyu
                 + gamma_par * v0 * eyu;
    }
  }
}

/* ----------------------------------------------------------------------
   Chiral spinner with odd (Hall) viscosity:
       F_odd = eta_odd * (z_hat x v),
   producing a transverse Hall response.  Combined with frictional drag,
   tan(phi) = eta_odd * omega_spin / gamma.
------------------------------------------------------------------------- */

void FixActive::compute_spinner()
{
  double **f = atom->f;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double vx = v[i][0];
      double vy = v[i][1];
      f[i][0] += -eta_odd * omega_spin * vy;
      f[i][1] +=  eta_odd * omega_spin * vx;
    }
  }
}

/* ----------------------------------------------------------------------
   Inertial ABP (Lowen 2020, Mandal 2020):
       m r_dd = -gamma r_d + gamma v0 e_hat + sqrt(2 gamma kT) xi.
   Use this fix together with `fix nve` (NOT brownian/langevin) so that
   the Newtonian acceleration F/m drives the motion.  We add the drag,
   propulsion and thermal noise here; the orientation diffuses with Dr.
------------------------------------------------------------------------- */

void FixActive::compute_iabp()
{
  double **f = atom->f;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dt = update->dt;
  double noise_rot = sqrt(2.0 * Dr * dt);
  double noise_trans = sqrt(2.0 * gamma * kT / dt);

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      // Stokes drag.
      f[i][0] -= gamma * v[i][0];
      f[i][1] -= gamma * v[i][1];

      // Active self-propulsion expressed as a constant force.
      f[i][0] += gamma * v0 * cos(theta[i][0]);
      f[i][1] += gamma * v0 * sin(theta[i][0]);

      // Thermal noise (Einstein relation, fluctuation-dissipation).
      f[i][0] += noise_trans * random->gaussian();
      f[i][1] += noise_trans * random->gaussian();

      // Orientation diffusion.
      theta[i][0] += noise_rot * random->gaussian();
    }
  }
}

/* ----------------------------------------------------------------------
   MIPS (Cates-Tailleur quorum-sensing realisation): the local density
   slows the swimmer,  v(rho) = v0 * max(1 - rho/rho_max, 0).
   The local density is estimated from the neighbor list as the count of
   neighbors within Rcut divided by the disk/sphere volume.
------------------------------------------------------------------------- */

void FixActive::compute_mips()
{
  if (!list) return;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dt = update->dt;
  double Rcut2 = Rcut * Rcut;
  double noise_rot = sqrt(2.0 * Dr * dt);

  double volume;
  if (dimension == 2)
    volume = MY_PI * Rcut2;
  else
    volume = (4.0 / 3.0) * MY_PI * Rcut * Rcut * Rcut;

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    if (i >= nlocal) continue;

    int *jlist = firstneigh[i];
    int jnum = numneigh[i];

    int count = 1;          // include the central particle
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj] & NEIGHMASK;
      double dx = x[i][0] - x[j][0];
      double dy = x[i][1] - x[j][1];
      double dz = (dimension == 3) ? x[i][2] - x[j][2] : 0.0;
      double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 < Rcut2) count++;
    }
    double rho_local = count / volume;

    double v_eff = v0 * (1.0 - rho_local / rho_max);
    if (v_eff < 0.0) v_eff = 0.0;

    theta[i][0] += noise_rot * random->gaussian();
    f[i][0] += gamma * v_eff * cos(theta[i][0]);
    f[i][1] += gamma * v_eff * sin(theta[i][0]);
  }
}

/* ----------------------------------------------------------------------
   Active nematic (dry, particle-based)

   Model definition.
     Each particle carries a 2D body axis e_i = (cos theta_i, sin theta_i)
     that is head-tail symmetric (theta_i and theta_i + pi describe the
     same physical state).  Over a neighborhood of radius Rcut we
     construct the local 2D nematic Q-tensor

         Q_xx(i) = (1/N_i) sum_{j in B_i} ( cos^2 theta_j - 1/2 ),
         Q_xy(i) = (1/N_i) sum_{j in B_i} ( cos theta_j sin theta_j ),

     where B_i = { j : |r_j - r_i| < Rcut } u {i} and N_i = |B_i|.  The
     local director angle is

         phi_i = (1/2) atan2( 2 Q_xy(i) , 2 Q_xx(i) ),

     defined modulo pi.  The equations of motion are

         dtheta_i/dt = -K sin( 2 (theta_i - phi_i) ) + sqrt(2 D_r) eta_i ,
         d r_i  / dt =  alpha e_i  +  mu F_i ,

     where K is the apolar alignment stiffness (`K`), alpha the active
     stress amplitude (`activity`), D_r the rotational diffusion (`Dr`),
     and F_i the conservative inter-particle force from any pair_style
     in use.  The torque is pi-periodic, encoding the head-tail symmetry.

   Numerical scheme.
     Theta is integrated with the explicit Euler-Maruyama rule

         theta_i(t+dt) = theta_i(t)
                         - K sin(2 (theta_i(t) - phi_i(t))) dt
                         + sqrt(2 D_r dt) Z_i ,    Z_i ~ N(0,1) ,

     using the orientations at the start of the step (computed into a
     scratch buffer to keep the update synchronous).  The active force
     gamma * alpha * e_i(t+dt) is added to atom->f, so when combined
     with `fix brownian`/`fix langevin` the centroidal velocity is
     alpha e_i + mu F_i as required.
------------------------------------------------------------------------- */

void FixActive::compute_nematic()
{
  if (!list) return;

  // Same rationale as in compute_vicsek: forward-communicate theta to
  // ghosts so that the local Q-tensor sees all neighbours within Rcut,
  // including those owned by other MPI ranks.
  comm->forward_comm(this);

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dt = update->dt;
  double Rcut2 = Rcut * Rcut;
  double noise_rot = sqrt(2.0 * Dr * dt);

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  double *new_theta = new double[nlocal];
  for (int i = 0; i < nlocal; i++) new_theta[i] = theta[i][0];

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    if (i >= nlocal) continue;

    double exi = cos(theta[i][0]);
    double eyi = sin(theta[i][0]);
    double Qxx = exi * exi - 0.5;
    double Qxy = exi * eyi;
    int count = 1;

    int *jlist = firstneigh[i];
    int jnum = numneigh[i];
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj] & NEIGHMASK;
      // j may be a ghost; its theta was populated by the forward_comm above.

      double dx = x[i][0] - x[j][0];
      double dy = x[i][1] - x[j][1];
      if (dx * dx + dy * dy < Rcut2) {
        double exj = cos(theta[j][0]);
        double eyj = sin(theta[j][0]);
        Qxx += exj * exj - 0.5;
        Qxy += exj * eyj;
        count++;
      }
    }
    Qxx /= count;
    Qxy /= count;

    // Director angle (defined modulo pi).
    double phi = 0.5 * atan2(2.0 * Qxy, 2.0 * Qxx);

    // Apolar aligning torque; sin(2 dtheta) is automatically pi-periodic
    // so explicit wrapping is not required.
    double dtheta_ang = theta[i][0] - phi;
    double tau_align = -K * sin(2.0 * dtheta_ang);

    new_theta[i] = theta[i][0] + tau_align * dt + noise_rot * random->gaussian();

    // Polar self-propulsion along the body axis.
    f[i][0] += gamma * activity * cos(new_theta[i]);
    f[i][1] += gamma * activity * sin(new_theta[i]);
  }

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) theta[i][0] = new_theta[i];

  delete[] new_theta;
}

/* ---------------------------------------------------------------------- */

double FixActive::memory_usage()
{
  double bytes = 0.0;
  bytes += (double) nmax_alloc * 1 * sizeof(double);    // theta
  bytes += (double) nmax_alloc * 4 * sizeof(double);    // quat
  bytes += (double) nmax_alloc * 3 * sizeof(double);    // va
  return bytes;
}
