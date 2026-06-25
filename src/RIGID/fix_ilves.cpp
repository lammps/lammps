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
   ILVES bond/angle constraint solver.

   This fix is the LAMMPS interface to a port of the GROMACS 2021 reference
   implementation of the ILVES algorithm
   (https://github.com/LorienLV/_PAPER_ILVES, LGPL-2.1).  ILVES enforces
   holonomic distance constraints with Newton's method on a sparse system,
   solved by a parallel Schur-complement direct solver.  Unlike fix shake,
   it handles arbitrarily large connected constraint clusters.

   Reference:
     L. Lopez-Villellas, C. C. K. Mikkelsen, et al., "ILVES: Accurate and
     Efficient Bond Length and Angle Constraints in Molecular Dynamics",
     J. Chem. Theory Comput. 21, 8711-8719 (2025),
     doi:10.1021/acs.jctc.5c01376

   Contributing author: Axel Kohlmeyer (Temple U), with Claude Code (Opus 4.8),
   under the direction of the ILVES authors and following fix shake / fix
   rattle / fix restrain as LAMMPS-side references.
------------------------------------------------------------------------- */

#include "fix_ilves.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "ilves_asym.h"
#include "ilves_sym.h"
#include "math_const.h"
#include "memory.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

namespace {
enum { ILVES_FAST, ILVES_FULL };
enum { LINEAR_ERROR, LINEAR_SKIP, LINEAR_RESTRAIN };
constexpr double MASSDELTA = 0.1;
// default force constant (kcal/mol/Angstrom^2-equivalent, in the active unit
// system) for the near-linear-angle harmonic A-C restraint when the user does
// not set one with the kbond keyword
constexpr double KBOND_AUTO = 1000.0;

// keywords that terminate the b/a/t/m selector lists
int is_keyword(const char *s)
{
  return ((strcmp(s, "variant") == 0) || (strcmp(s, "store") == 0) ||
          (strcmp(s, "mode") == 0) || (strcmp(s, "kbond") == 0) ||
          (strcmp(s, "linearangle") == 0));
}

int is_selector(const char *s)
{
  return ((strcmp(s, "b") == 0) || (strcmp(s, "a") == 0) || (strcmp(s, "t") == 0) ||
          (strcmp(s, "m") == 0));
}
}    // namespace

static const char cite_fix_ilves[] =
    "fix ilves command: doi:10.1021/acs.jctc.5c01376\n\n"
    "@Article{LopezVillellas2025,\n"
    " author = {L{\\'o}pez-Villellas, Lori{\\'e}n and Mikkelsen, Carl Christian Kjelgaard and "
    "Galano-Frutos, Juan Jos{\\'e} and Marco-Sola, Santiago and Alastruey-Bened{\\'e}, Jes{\\'u}s "
    "and Ib{\\'a}{\\~n}ez, Pablo and Echenique, Pablo and Moret{\\'o}, Miquel and {De Rosa}, Maria "
    "Cristina and Garc{\\'i}a-Risue{\\~n}o, Pablo},\n"
    " title = {ILVES: Accurate and Efficient Bond Length and Angle Constraints in Molecular "
    "Dynamics},\n"
    " journal = {J.~Chem.\\ Theory Comput.},\n"
    " volume = 21,\n"
    " pages = {8711--8719},\n"
    " year = 2025\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

FixIlves::FixIlves(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), tolerance(1.0e-4), max_iter(25), output_every(0), next_output(0),
    variant(ILVES_FAST), fixed_iter(0), linear_mode(LINEAR_ERROR), linear_threshold(175.0),
    kbond(-1.0), molecular(0), types_negated(0), store_flag(0),
    fstore(nullptr),
    maxstore(0), niter_max(0), nconstraints(0), ilves_solver(nullptr),
    xpred(nullptr), xpred0(nullptr), dx(nullptr), maxatom(0), commstage(0), dtv(0.0), dtfsq(0.0),
    inv_dtfsq(0.0), x(nullptr), v(nullptr), f(nullptr), mass(nullptr), rmass(nullptr),
    type(nullptr), mask(nullptr), nlocal(0)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_ilves);

  // predicted positions and per-iteration increments are communicated to/from
  // ghost atoms each Newton iteration (forward: positions, reverse: increments)

  comm_forward = 3;
  comm_reverse = 3;

  // this fix removes degrees of freedom (one per constraint); tell the
  // temperature computes to query dof()
  dof_flag = 1;

  // the constraint forces can contribute to the global pressure virial; this is
  // opt-in via fix_modify virial yes (the post-restart and per-atom-mass virial
  // are not yet reproducible enough to enable by default)
  virial_global_flag = 1;
  for (int i = 0; i < 6; ++i) virial[i] = 0.0;

  // the linearangle restrain substitute is a real potential; expose its energy
  // as the fix's global scalar (opt-in to the thermodynamic output via
  // fix_modify energy yes, as for other restraint-style fixes)
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  energy_global_flag = 1;
  thermo_energy = 1;    // count the restraint / minimization energy in the PE by default
  erestraint = 0.0;

  if (narg < 7) utils::missing_cmd_args(FLERR, "fix ilves", error);

  molecular = atom->molecular;
  if (molecular == Atom::ATOMIC)
    error->all(FLERR, "Fix ilves requires a molecular atom style with bond topology");
  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, "Fix ilves requires an atom map, see atom_modify");

  tolerance = utils::numeric(FLERR, arg[3], false, lmp);
  max_iter = utils::inumeric(FLERR, arg[4], false, lmp);
  output_every = utils::inumeric(FLERR, arg[5], false, lmp);
  if (tolerance <= 0.0) error->all(FLERR, "Fix ilves tolerance must be > 0");
  if (max_iter <= 0) error->all(FLERR, "Fix ilves iteration count must be > 0");

  // size selector membership tables (1-based type indexing)

  bond_flag.assign(atom->nbondtypes + 1, 0);
  angle_flag.assign(atom->nangletypes + 1, 0);
  type_flag.assign(atom->ntypes + 1, 0);

  // parse one or more b/a/t/m selector lists, then optional keyword/value pairs

  int iarg = 6;
  while ((iarg < narg) && !is_keyword(arg[iarg])) {
    if (!is_selector(arg[iarg]))
      error->all(FLERR, "Unknown fix ilves selector or keyword: {}", arg[iarg]);
    const char sel = arg[iarg][0];
    ++iarg;
    int nvalues = 0;
    while ((iarg < narg) && !is_selector(arg[iarg]) && !is_keyword(arg[iarg])) {
      if (sel == 'm') {
        mass_list.push_back(utils::numeric(FLERR, arg[iarg], false, lmp));
      } else {
        int v = utils::inumeric(FLERR, arg[iarg], false, lmp);
        if (sel == 'b') {
          if ((v < 1) || (v > atom->nbondtypes))
            error->all(FLERR, "Invalid fix ilves bond type {}", v);
          bond_flag[v] = 1;
        } else if (sel == 'a') {
          if ((v < 1) || (v > atom->nangletypes))
            error->all(FLERR, "Invalid fix ilves angle type {}", v);
          angle_flag[v] = 1;
        } else if (sel == 't') {
          if ((v < 1) || (v > atom->ntypes))
            error->all(FLERR, "Invalid fix ilves atom type {}", v);
          type_flag[v] = 1;
        }
      }
      ++iarg;
      ++nvalues;
    }
    if (nvalues == 0) error->all(FLERR, "Fix ilves selector '{}' needs one or more values", sel);
  }

  // optional keyword/value pairs

  while (iarg < narg) {
    if (strcmp(arg[iarg], "variant") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix ilves variant", error);
      if (strcmp(arg[iarg + 1], "fast") == 0)
        variant = ILVES_FAST;
      else if (strcmp(arg[iarg + 1], "full") == 0)
        variant = ILVES_FULL;
      else
        error->all(FLERR, "Unknown fix ilves variant: {}", arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "store") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix ilves store", error);
      store_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "mode") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix ilves mode", error);
      if (strcmp(arg[iarg + 1], "converge") == 0)
        fixed_iter = 0;
      else if (strcmp(arg[iarg + 1], "fixed") == 0)
        fixed_iter = 1;
      else
        error->all(FLERR, "Unknown fix ilves mode: {}", arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "linearangle") == 0) {
      if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "fix ilves linearangle", error);
      if (strcmp(arg[iarg + 1], "error") == 0)
        linear_mode = LINEAR_ERROR;
      else if (strcmp(arg[iarg + 1], "skip") == 0)
        linear_mode = LINEAR_SKIP;
      else if (strcmp(arg[iarg + 1], "restrain") == 0)
        linear_mode = LINEAR_RESTRAIN;
      else
        error->all(FLERR, "Unknown fix ilves linearangle mode: {}", arg[iarg + 1]);
      linear_threshold = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      if ((linear_threshold <= 0.0) || (linear_threshold >= 180.0))
        error->all(FLERR, "Fix ilves linearangle threshold must be between 0 and 180 degrees");
      iarg += 3;
    } else if (strcmp(arg[iarg], "kbond") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix ilves kbond", error);
      kbond = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (kbond <= 0.0) error->all(FLERR, "Fix ilves kbond must be > 0");
      iarg += 2;
    } else {
      error->all(FLERR, "Unknown fix ilves keyword: {}", arg[iarg]);
    }
  }

  // with store yes, expose the per-atom constraint forces as a 3-column array

  if (store_flag) {
    peratom_flag = 1;
    size_peratom_cols = 3;
    peratom_freq = 1;
  }

  // require at least one selector

  int any = (int) mass_list.size();
  for (int i = 1; i <= atom->nbondtypes; ++i) any += bond_flag[i];
  for (int i = 1; i <= atom->ntypes; ++i) any += type_flag[i];
  if (any == 0) error->all(FLERR, "Fix ilves requires at least one b/t/m selector");

  next_output = 0;
}

/* ---------------------------------------------------------------------- */

FixIlves::~FixIlves()
{
  // restore the bond/angle types we negated so the bonded styles act on them again

  if (types_negated) {
    if (atom->bond_type) negate_bond_types(1);
    if (atom->angle_type) negate_angle_types(1);
  }

  delete ilves_solver;
  memory->destroy(xpred);
  memory->destroy(xpred0);
  memory->destroy(dx);
  memory->destroy(fstore);
}

/* ---------------------------------------------------------------------- */

int FixIlves::setmask()
{
  int mask = 0;
  mask |= PRE_NEIGHBOR;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  mask |= MIN_PRE_NEIGHBOR;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixIlves::init()
{
  if (!force->bond)
    error->all(FLERR, "Fix ilves requires a bond style to define equilibrium bond lengths");

  // r-RESPA is not supported yet
  if (utils::strmatch(update->integrate_style, "^respa"))
    error->all(FLERR, "Fix ilves does not support run_style respa");

  // equilibrium bond lengths per bond type, as in fix shake

  bond_distance.assign(atom->nbondtypes + 1, 0.0);
  for (int i = 1; i <= atom->nbondtypes; ++i)
    bond_distance[i] = force->bond->equilibrium_distance(i);

  // equilibrium A-C "virtual bond" distances for the selected angle types,
  // from the two leg bond lengths and the angle via the law of cosines (as in
  // fix shake).  computed per angle type from a representative angle and
  // reduced across ranks so every rank agrees.

  angle_distance.assign(atom->nangletypes + 1, 0.0);
  angle_linear.assign(atom->nangletypes + 1, 0);
  int any_angle = 0;
  for (int i = 1; i <= atom->nangletypes; ++i) any_angle += angle_flag[i];
  if (any_angle) {
    if (!force->angle)
      error->all(FLERR, "Fix ilves angle constraints require an angle style");

    // classify near-linear angle types (theta0 >= linear_threshold).  the A-C
    // virtual bond becomes rank-deficient near 180 degrees, so these are handled
    // per linear_mode instead of by an ordinary distance constraint.
    const double thresh_rad = linear_threshold * MathConst::MY_PI / 180.0;
    std::string linear_types;
    for (int i = 1; i <= atom->nangletypes; ++i) {
      if (!angle_flag[i]) continue;
      const double th = force->angle->equilibrium_angle(i);
      if (th >= thresh_rad) {
        angle_linear[i] = 1;
        linear_types += " " + std::to_string(i);
      }
    }
    if (!linear_types.empty()) {
      if (linear_mode == LINEAR_ERROR)
        error->all(FLERR, "Fix ilves angle type(s){} have an equilibrium angle at or above the "
                          "linearangle threshold of {} degrees and cannot be rigidly constrained; "
                          "use the linearangle keyword to skip or restrain them",
                   linear_types, linear_threshold);
      else if (comm->me == 0)
        error->warning(FLERR, "Fix ilves treating near-linear angle type(s){} with linearangle "
                              "mode {}", linear_types,
                       (linear_mode == LINEAR_SKIP) ? "skip" : "restrain");
    }

    std::vector<double> ad(atom->nangletypes + 1, 0.0);
    int *num_angle = atom->num_angle;
    int **angle_type = atom->angle_type;
    tagint **angle_atom1 = atom->angle_atom1;
    tagint **angle_atom2 = atom->angle_atom2;
    tagint **angle_atom3 = atom->angle_atom3;
    tagint *tag = atom->tag;
    int n = atom->nlocal;

    if (num_angle) {
      for (int i = 0; i < n; ++i) {
        for (int m = 0; m < num_angle[i]; ++m) {
          const int at = abs(angle_type[i][m]);
          if ((at == 0) || !angle_flag[at] || (ad[at] > 0.0)) continue;
          if (angle_atom2[i][m] != tag[i]) continue;    // process at the center atom
          const int a0 = atom->map(angle_atom1[i][m]);
          const int c0 = atom->map(angle_atom3[i][m]);
          if ((a0 < 0) || (c0 < 0)) continue;
          const int tab = find_bond_type(i, a0);
          const int tbc = find_bond_type(i, c0);
          if ((tab == 0) || (tbc == 0)) continue;
          const double th = force->angle->equilibrium_angle(at);
          const double b1 = bond_distance[tab], b2 = bond_distance[tbc];
          ad[at] = sqrt(b1 * b1 + b2 * b2 - 2.0 * b1 * b2 * cos(th));
        }
      }
    }
    MPI_Allreduce(ad.data(), angle_distance.data(), atom->nangletypes + 1, MPI_DOUBLE, MPI_MAX,
                  world);
  }

  // timestep factors for predicting unconstrained positions and converting
  // the Lagrange multipliers to forces, identical to fix shake (SHAKE form).

  dtv = update->dt;
  dtfsq = update->dt * update->dt * force->ftm2v;
  inv_dtfsq = 1.0 / dtfsq;
}

/* ---------------------------------------------------------------------- */

void FixIlves::setup_pre_neighbor()
{
  // negate the constrained bond types once, before the first neighbor list
  // (and thus the bond list) is built, so the bonded styles skip them from
  // the very first force evaluation.  the negated sign travels with migrating
  // atoms, so this is done only once.  build the initial constraint list here
  // too, since the atom map and ghosts are already current at this point.

  if (!types_negated) {
    negate_bond_types(-1);
    negate_angle_types(-1);
    types_negated = 1;
  }

  build_constraint_list();
}

/* ---------------------------------------------------------------------- */

void FixIlves::setup(int vflag)
{
  bigint nb = 0, na = 0;
  for (int k = 0; k < nconstraints; ++k)
    if (clist_btype[k] > 0) ++nb; else ++na;
  bigint nc[2] = {nb, na}, nctot[2] = {0, 0};
  MPI_Allreduce(nc, nctot, 2, MPI_LMP_BIGINT, MPI_SUM, world);
  if (comm->me == 0)
    utils::logmesg(lmp, "Fix ilves: constraining {} bond(s) and {} angle(s)\n", nctot[0], nctot[1]);

  // schedule the next statistics output

  const bigint ntimestep = update->ntimestep;
  if (output_every) {
    next_output = ntimestep + output_every;
    if (ntimestep % output_every != 0)
      next_output = (ntimestep / output_every) * output_every + output_every;
  } else {
    next_output = -1;
  }
  if (output_every) stats();

  // project the initial velocities onto the constraint manifold (remove the
  // component along each bond), so the run starts constraint-consistent
  project_velocities();

  post_force(vflag);
}

/* ----------------------------------------------------------------------
   minimization setup.  The constraint list and the bond/angle type negation
   were prepared in setup_pre_neighbor() (called before the neighbor build in
   minimization as well as dynamics), so here we only apply the harmonic
   restraint substitute used in place of the holonomic constraints, which have
   no meaning without time integration.
------------------------------------------------------------------------- */

void FixIlves::min_setup(int vflag)
{
  bigint nb = 0, na = 0;
  for (int k = 0; k < nconstraints; ++k)
    if (clist_btype[k] > 0) ++nb; else ++na;
  bigint nc[2] = {nb, na}, nctot[2] = {0, 0};
  MPI_Allreduce(nc, nctot, 2, MPI_LMP_BIGINT, MPI_SUM, world);
  if (comm->me == 0)
    utils::logmesg(lmp, "Fix ilves: replacing {} bond and {} angle constraint(s) with harmonic "
                        "restraints for minimization\n", nctot[0], nctot[1]);
  min_post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixIlves::end_of_step()
{
  // RATTLE-style velocity constraint: remove the relative velocity along each
  // bond after the final velocity update
  project_velocities();
}

/* ---------------------------------------------------------------------- */

void FixIlves::pre_neighbor()
{
  build_constraint_list();
}

/* ---------------------------------------------------------------------- */

void FixIlves::post_force(int vflag)
{
  x = atom->x;
  v = atom->v;
  f = atom->f;
  type = atom->type;
  mass = atom->mass;
  rmass = atom->rmass;
  nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;

  if (output_every && (update->ntimestep == next_output)) {
    stats();
    next_output += output_every;
  }

  v_init(vflag);

  grow_arrays_local();

  // (re)allocate and zero the per-atom constraint-force store

  if (store_flag) {
    if (nlocal > maxstore) {
      memory->destroy(fstore);
      maxstore = atom->nmax;
      memory->create(fstore, maxstore, 3, "ilves:fstore");
      array_atom = fstore;
    }
    for (int i = 0; i < nlocal; ++i) fstore[i][0] = fstore[i][1] = fstore[i][2] = 0.0;
  }

  // add the near-linear-angle restraint forces (linearangle restrain) to
  // atom->f before predicting positions, so they enter the dynamics like any
  // other force.  uses the dx buffer as scratch (zeroed and reused by the solve)
  if (linear_mode == LINEAR_RESTRAIN) apply_linear_restraint();

  // predict the unconstrained positions for owned atoms (fix shake style) and
  // save them, so the net constrained displacement can be turned into a force

  for (int i = 0; i < nlocal; ++i) {
    const double m = rmass ? rmass[i] : mass[type[i]];
    if (m > 0.0) {
      const double dtfm = dtfsq / m;
      xpred[i][0] = x[i][0] + dtv * v[i][0] + dtfm * f[i][0];
      xpred[i][1] = x[i][1] + dtv * v[i][1] + dtfm * f[i][1];
      xpred[i][2] = x[i][2] + dtv * v[i][2] + dtfm * f[i][2];
    } else {
      xpred[i][0] = x[i][0];
      xpred[i][1] = x[i][1];
      xpred[i][2] = x[i][2];
    }
    xpred0[i][0] = xpred[i][0];
    xpred0[i][1] = xpred[i][1];
    xpred0[i][2] = xpred[i][2];
  }

  // zero the per-iteration increment buffer (home + ghost)

  for (int i = 0; i < nall; ++i) dx[i][0] = dx[i][1] = dx[i][2] = 0.0;

  // Global Newton iteration.  The loop is driven uniformly on all ranks: the
  // convergence test uses the all-reduced maximum residual, so every rank takes
  // the same number of iterations and participates in the same collective
  // communication even if it owns no constraints.  Cross-rank coupling is a
  // block-Jacobi sweep: each iteration reverse-sums the per-atom position
  // increments to their owners, applies them, and forward-communicates the
  // predicted positions back to the ghosts.

  commstage = 0;    // forward-comm predicted positions (with PBC shift)
  comm->forward_comm(this);
  double local = ilves_solver ? ilves_solver->prepare(x, xpred) : 0.0;
  double ptau = 0.0;
  if (!fixed_iter) MPI_Allreduce(&local, &ptau, 1, MPI_DOUBLE, MPI_MAX, world);

  int numit = 0;
  for (int i = 0; i < max_iter; ++i) {
    // in convergence mode (the default) stop once the global maximum relative
    // violation is below the tolerance; in fixed mode always run max_iter steps
    // (which avoids the per-iteration MPI reduction)
    if (!fixed_iter && (!std::isfinite(ptau) || (ptau <= tolerance))) break;
    ++numit;
    if (ilves_solver) ilves_solver->step(dx);

    comm->reverse_comm(this);
    for (int k = 0; k < nlocal; ++k) {
      xpred[k][0] += dx[k][0];
      xpred[k][1] += dx[k][1];
      xpred[k][2] += dx[k][2];
    }
    for (int k = 0; k < nall; ++k) dx[k][0] = dx[k][1] = dx[k][2] = 0.0;
    comm->forward_comm(this);

    local = ilves_solver ? ilves_solver->recompute(x, xpred, i == 0) : 0.0;
    if (!fixed_iter) MPI_Allreduce(&local, &ptau, 1, MPI_DOUBLE, MPI_MAX, world);
  }

  if (numit > niter_max) niter_max = numit;

  // convert the net constrained displacement of each owned atom into a force,
  // identical to the multiplier-to-force coupling of fix shake (f += m*dx/dtfsq)

  if (numit > 0) {
    for (int i = 0; i < nlocal; ++i) {
      const double m = rmass ? rmass[i] : mass[type[i]];
      if (m <= 0.0) continue;
      const double fac = m * inv_dtfsq;
      const double fcx = fac * (xpred[i][0] - xpred0[i][0]);
      const double fcy = fac * (xpred[i][1] - xpred0[i][1]);
      const double fcz = fac * (xpred[i][2] - xpred0[i][2]);
      f[i][0] += fcx;
      f[i][1] += fcy;
      f[i][2] += fcz;
      if (store_flag) {
        fstore[i][0] = fcx;
        fstore[i][1] = fcy;
        fstore[i][2] = fcz;
      }
    }
    if (vflag_global) ilves_solver->add_global_virial(virial, inv_dtfsq);
  }
}

/* ----------------------------------------------------------------------
   RATTLE-style velocity projection: remove the component of relative velocity
   along each constrained bond.  The velocity constraint is linear, so this is a
   block-Jacobi sweep (exact in one pass for an isolated bond) driven uniformly
   across ranks by the all-reduced residual.  Reuses the xpred buffer as the
   velocity work array (forward-communicated without a PBC shift) and dx as the
   per-atom increment accumulator.
------------------------------------------------------------------------- */

void FixIlves::project_velocities()
{
  x = atom->x;
  v = atom->v;
  nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;

  grow_arrays_local();

  for (int i = 0; i < nlocal; ++i) {
    xpred[i][0] = v[i][0];
    xpred[i][1] = v[i][1];
    xpred[i][2] = v[i][2];
  }
  commstage = 1;
  comm->forward_comm(this);

  for (int iter = 0; iter < max_iter; ++iter) {
    for (int i = 0; i < nall; ++i) dx[i][0] = dx[i][1] = dx[i][2] = 0.0;

    double local = 0.0;
    for (int k = 0; k < nconstraints; ++k) {
      const int a = clist_a[k], b = clist_b[k];
      const double rx = x[b][0] - x[a][0];
      const double ry = x[b][1] - x[a][1];
      const double rz = x[b][2] - x[a][2];
      const double rr = rx * rx + ry * ry + rz * rz;
      const double vrel = (xpred[b][0] - xpred[a][0]) * rx + (xpred[b][1] - xpred[a][1]) * ry +
          (xpred[b][2] - xpred[a][2]) * rz;
      const double res = fabs(vrel) / rr;
      if (res > local) local = res;
      const double mu = vrel / ((invmass[a] + invmass[b]) * rr);
      dx[a][0] += mu * invmass[a] * rx;
      dx[a][1] += mu * invmass[a] * ry;
      dx[a][2] += mu * invmass[a] * rz;
      dx[b][0] -= mu * invmass[b] * rx;
      dx[b][1] -= mu * invmass[b] * ry;
      dx[b][2] -= mu * invmass[b] * rz;
    }

    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, world);
    if (global < tolerance) break;

    comm->reverse_comm(this);
    for (int i = 0; i < nlocal; ++i) {
      xpred[i][0] += dx[i][0];
      xpred[i][1] += dx[i][1];
      xpred[i][2] += dx[i][2];
    }
    commstage = 1;
    comm->forward_comm(this);
  }

  for (int i = 0; i < nlocal; ++i) {
    v[i][0] = xpred[i][0];
    v[i][1] = xpred[i][1];
    v[i][2] = xpred[i][2];
  }
}

/* ----------------------------------------------------------------------
   (re)allocate the predicted-position / increment buffers to hold local+ghost
------------------------------------------------------------------------- */

void FixIlves::grow_arrays_local()
{
  if (atom->nmax > maxatom) {
    memory->destroy(xpred);
    memory->destroy(xpred0);
    memory->destroy(dx);
    maxatom = atom->nmax;
    memory->create(xpred, maxatom, 3, "ilves:xpred");
    memory->create(xpred0, maxatom, 3, "ilves:xpred0");
    memory->create(dx, maxatom, 3, "ilves:dx");
  }
}

/* ----------------------------------------------------------------------
   forward communication of predicted positions to ghosts (with PBC shift)
------------------------------------------------------------------------- */

int FixIlves::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int m = 0;
  // commstage 2 carries the per-atom inverse mass (one value, zero-padded to
  // the 3-wide buffer); never shifted
  if (commstage == 2) {
    for (int i = 0; i < n; ++i) {
      buf[m++] = invmass[list[i]];
      buf[m++] = 0.0;
      buf[m++] = 0.0;
    }
    return m;
  }
  // commstage 1 carries velocities (no periodic-image shift); commstage 0
  // carries predicted positions (shifted across periodic boundaries)
  if ((pbc_flag == 0) || (commstage == 1)) {
    for (int i = 0; i < n; ++i) {
      const int j = list[i];
      buf[m++] = xpred[j][0];
      buf[m++] = xpred[j][1];
      buf[m++] = xpred[j][2];
    }
  } else {
    double dxs, dys, dzs;
    if (domain->triclinic == 0) {
      dxs = pbc[0] * domain->xprd;
      dys = pbc[1] * domain->yprd;
      dzs = pbc[2] * domain->zprd;
    } else {
      dxs = pbc[0] * domain->xprd + pbc[5] * domain->xy + pbc[4] * domain->xz;
      dys = pbc[1] * domain->yprd + pbc[3] * domain->yz;
      dzs = pbc[2] * domain->zprd;
    }
    for (int i = 0; i < n; ++i) {
      const int j = list[i];
      buf[m++] = xpred[j][0] + dxs;
      buf[m++] = xpred[j][1] + dys;
      buf[m++] = xpred[j][2] + dzs;
    }
  }
  return m;
}

void FixIlves::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0;
  const int last = first + n;
  if (commstage == 2) {
    for (int i = first; i < last; ++i) {
      invmass[i] = buf[m++];
      m += 2;
    }
    return;
  }
  for (int i = first; i < last; ++i) {
    xpred[i][0] = buf[m++];
    xpred[i][1] = buf[m++];
    xpred[i][2] = buf[m++];
  }
}

/* ----------------------------------------------------------------------
   reverse communication of position increments (summed into the owners).
   increments are displacement vectors, so no PBC shift is applied.
------------------------------------------------------------------------- */

int FixIlves::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  const int last = first + n;
  for (int i = first; i < last; ++i) {
    buf[m++] = dx[i][0];
    buf[m++] = dx[i][1];
    buf[m++] = dx[i][2];
  }
  return m;
}

void FixIlves::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m = 0;
  for (int i = 0; i < n; ++i) {
    const int j = list[i];
    dx[j][0] += buf[m++];
    dx[j][1] += buf[m++];
    dx[j][2] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   rebuild the local list of constrained bonds from local + ghost topology.
   selection uses abs(bond_type) so it is independent of the negation state.
------------------------------------------------------------------------- */

void FixIlves::build_constraint_list()
{
  type = atom->type;
  mask = atom->mask;
  mass = atom->mass;
  rmass = atom->rmass;
  nlocal = atom->nlocal;

  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  tagint *tag = atom->tag;
  const int newton_bond = force->newton_bond;

  clist_a.clear();
  clist_b.clear();
  clist_d.clear();
  clist_btype.clear();
  rlist_a.clear();
  rlist_c.clear();
  rlist_d.clear();

  if (num_bond) {
    for (int i = 0; i < nlocal; ++i) {
      if (!(mask[i] & groupbit)) continue;
      for (int m = 0; m < num_bond[i]; ++m) {
        const int btype = abs(bond_type[i][m]);
        if (btype == 0) continue;
        int j = atom->map(bond_atom[i][m]);
        if (j < 0)
          error->one(FLERR, "Fix ilves bond atom missing on this processor; increase the "
                            "communication cutoff with comm_modify cutoff");
        j = domain->closest_image(i, j);
        // with newton_bond off the bond is stored on both atoms; keep one copy
        if (!newton_bond && (tag[i] > tag[j])) continue;
        if (!(mask[j] & groupbit)) continue;
        if (!bond_selected(i, j, btype)) continue;
        clist_a.push_back(i);
        clist_b.push_back(j);
        clist_d.push_back(bond_distance[btype]);
        clist_btype.push_back(btype);
      }
    }
  }

  // angle "virtual bond" constraints: the A-C distance of a selected angle
  // A-B-C whose two flanking bonds are themselves constrained, which makes the
  // triangle (and hence the angle) rigid.

  int *num_angle = atom->num_angle;
  if (num_angle) {
    for (int i = 0; i < nlocal; ++i) {
      if (!(mask[i] & groupbit)) continue;
      for (int m = 0; m < num_angle[i]; ++m) {
        int a, c, atype;
        if (!angle_selected(i, m, a, c, atype)) continue;
        if (angle_linear[atype]) {
          // near-linear angle: the A-C virtual bond is rank-deficient.  with
          // linearangle skip do nothing (the angle is left to the bonded style);
          // with restrain record the A-C pair for the stiff harmonic bond
          // substitute applied in post_force.
          if (linear_mode == LINEAR_RESTRAIN) {
            rlist_a.push_back(a);
            rlist_c.push_back(c);
            rlist_d.push_back(angle_distance[atype]);
          }
          continue;
        }
        clist_a.push_back(a);
        clist_b.push_back(c);
        clist_d.push_back(angle_distance[atype]);
        clist_btype.push_back(-atype);    // negative marks an A-C angle constraint
      }
    }
  }

  nconstraints = (int) clist_a.size();

  // per-atom inverse mass (1/m) for the constrained atoms, handed to the solver.
  // sized for local+ghost atoms since a constraint partner may be a ghost.

  const int nall = atom->nlocal + atom->nghost;
  invmass.assign(nall, 0.0);
  for (int i = 0; i < nlocal; ++i)
    invmass[i] = rmass ? 1.0 / rmass[i] : 1.0 / mass[type[i]];

  // communicate the inverse mass to the ghosts, so a constraint partner owned by
  // another rank (or a periodic image) has the correct value regardless of
  // whether the atom style communicates per-atom mass to ghosts.

  if (atom->nghost > 0) {
    commstage = 2;
    comm->forward_comm(this);
  }

  // (re)build the ILVES solver for the current constraint topology.
  // Phase 1 is serial single-thread: always one OpenMP thread.

  delete ilves_solver;
  ilves_solver = nullptr;
  if (nconstraints > 0) {
    if (variant == ILVES_FAST)
      ilves_solver = new ILVES::IlvesSym(lmp, nconstraints, clist_a.data(), clist_b.data(),
                                         clist_d.data(), invmass.data(), 1);
    else
      ilves_solver = new ILVES::IlvesAsym(lmp, nconstraints, clist_a.data(), clist_b.data(),
                                          clist_d.data(), invmass.data(), 1);
  }
}

/* ----------------------------------------------------------------------
   add the stiff harmonic A-C "virtual bond" restraint for near-linear angle
   types (linearangle restrain mode) to atom->f.  E = k (r - d)^2 with r the
   current A-C distance and d the law-of-cosines target; this replaces the
   rank-deficient holonomic A-C constraint with a well-behaved distance
   restraint.  forces are accumulated into the dx buffer (home + ghost) and
   reverse-summed to the owning ranks, then added to atom->f, so the restraint
   is part of the force used to predict positions for the constraint solve.
------------------------------------------------------------------------- */

void FixIlves::apply_linear_restraint()
{
  const double k = (kbond > 0.0) ? kbond : KBOND_AUTO;
  const int nall = atom->nlocal + atom->nghost;

  for (int i = 0; i < nall; ++i) dx[i][0] = dx[i][1] = dx[i][2] = 0.0;

  erestraint = 0.0;
  const int nr = (int) rlist_a.size();
  for (int kk = 0; kk < nr; ++kk)
    erestraint += min_harmonic_bond(rlist_a[kk], rlist_c[kk], rlist_d[kk], k);

  comm->reverse_comm(this);    // sum ghost restraint forces into their owners

  double **ff = atom->f;
  for (int i = 0; i < atom->nlocal; ++i) {
    ff[i][0] += dx[i][0];
    ff[i][1] += dx[i][1];
    ff[i][2] += dx[i][2];
  }
}

/* ----------------------------------------------------------------------
   energy minimization: a holonomic constraint cannot be enforced without time
   integration, so during minimization each constraint is replaced by a stiff
   harmonic bond E = k (r - d)^2 (as in fix shake).  The bond constraints use
   the kbond force constant (default 1e9*boltz, very stiff, as fix shake) and
   the near-linear-angle A-C restraints use the softer restrain force constant.
   Forces are accumulated for home + ghost atoms in the dx buffer, reverse-summed
   to the owners, and added to atom->f; the total energy is the fix scalar.
------------------------------------------------------------------------- */

void FixIlves::min_post_force(int /*vflag*/)
{
  x = atom->x;
  type = atom->type;
  nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;

  grow_arrays_local();

  if (store_flag) {
    if (nlocal > maxstore) {
      memory->destroy(fstore);
      maxstore = atom->nmax;
      memory->create(fstore, maxstore, 3, "ilves:fstore");
      array_atom = fstore;
    }
    for (int i = 0; i < nlocal; ++i) fstore[i][0] = fstore[i][1] = fstore[i][2] = 0.0;
  }

  const double k_bond = (kbond > 0.0) ? kbond : 1.0e9 * force->boltz;
  const double k_lin = (kbond > 0.0) ? kbond : KBOND_AUTO;

  for (int i = 0; i < nall; ++i) dx[i][0] = dx[i][1] = dx[i][2] = 0.0;
  erestraint = 0.0;

  // stiff harmonic bonds replacing the holonomic constraints, plus the
  // (softer) near-linear-angle A-C restraints
  const int nc = nconstraints;
  for (int kk = 0; kk < nc; ++kk)
    erestraint += min_harmonic_bond(clist_a[kk], clist_b[kk], clist_d[kk], k_bond);
  const int nr = (int) rlist_a.size();
  for (int kk = 0; kk < nr; ++kk)
    erestraint += min_harmonic_bond(rlist_a[kk], rlist_c[kk], rlist_d[kk], k_lin);

  comm->reverse_comm(this);    // sum ghost restraint forces into their owners

  double **ff = atom->f;
  for (int i = 0; i < nlocal; ++i) {
    ff[i][0] += dx[i][0];
    ff[i][1] += dx[i][1];
    ff[i][2] += dx[i][2];
    if (store_flag) {
      fstore[i][0] = dx[i][0];
      fstore[i][1] = dx[i][1];
      fstore[i][2] = dx[i][2];
    }
  }
}

/* ----------------------------------------------------------------------
   accumulate the harmonic-bond force E = k (r - d)^2 between atoms a and b
   (local or ghost) into the dx buffer and return the bond energy.
------------------------------------------------------------------------- */

double FixIlves::min_harmonic_bond(int a, int b, double d, double k)
{
  double **xx = atom->x;
  const double ux = xx[b][0] - xx[a][0];
  const double uy = xx[b][1] - xx[a][1];
  const double uz = xx[b][2] - xx[a][2];
  const double r = sqrt(ux * ux + uy * uy + uz * uz);
  if (r < 1.0e-10) return 0.0;
  const double dr = r - d;
  const double fac = 2.0 * k * dr / r;    // F_a = fac*(x_b - x_a), F_b = -F_a
  dx[a][0] += fac * ux; dx[a][1] += fac * uy; dx[a][2] += fac * uz;
  dx[b][0] -= fac * ux; dx[b][1] -= fac * uy; dx[b][2] -= fac * uz;
  return k * dr * dr;
}

/* ----------------------------------------------------------------------
   negate (sign<0) or restore (sign>0) the bond_type of selected bonds.
   only flips bonds that are currently the wrong sign, so it is idempotent.
------------------------------------------------------------------------- */

void FixIlves::negate_bond_types(int sign)
{
  int *atype = atom->type;
  int *amask = atom->mask;
  int n = atom->nlocal;
  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;

  if (!num_bond) return;

  type = atype;    // bond_selected() reads these
  mask = amask;
  mass = atom->mass;
  rmass = atom->rmass;

  for (int i = 0; i < n; ++i) {
    if (!(amask[i] & groupbit)) continue;
    for (int m = 0; m < num_bond[i]; ++m) {
      const int btype = abs(bond_type[i][m]);
      if (btype == 0) continue;
      const int j = atom->map(bond_atom[i][m]);
      if (j < 0) continue;
      if (!(amask[j] & groupbit)) continue;
      if (!bond_selected(i, j, btype)) continue;
      if ((sign < 0) && (bond_type[i][m] > 0))
        bond_type[i][m] = -bond_type[i][m];
      else if ((sign > 0) && (bond_type[i][m] < 0))
        bond_type[i][m] = -bond_type[i][m];
    }
  }
}

/* ----------------------------------------------------------------------
   return the (positive) type of the bond between atoms i and j, or 0 if there
   is no such bond available locally.  searches the bond list of whichever of
   the two atoms is a local (owned) atom.
------------------------------------------------------------------------- */

int FixIlves::find_bond_type(int i, int j)
{
  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  tagint *tag = atom->tag;
  const int n = atom->nlocal;

  if (!num_bond) return 0;

  if (i < n) {
    const tagint tj = tag[j];
    for (int m = 0; m < num_bond[i]; ++m)
      if (bond_atom[i][m] == tj) return abs(bond_type[i][m]);
  }
  if (j < n) {
    const tagint ti = tag[i];
    for (int m = 0; m < num_bond[j]; ++m)
      if (bond_atom[j][m] == ti) return abs(bond_type[j][m]);
  }
  return 0;
}

/* ----------------------------------------------------------------------
   decide whether angle m of (local center) atom i is constrained.  it is when
   its type is selected, all three atoms are in the group, both flanking bonds
   are themselves selected (constrained), and the A-C distance is known.  fills
   the (local/ghost, closest-image) outer-atom indices a, c and the angle type.
------------------------------------------------------------------------- */

int FixIlves::angle_selected(int i, int m, int &a, int &c, int &atype)
{
  int **angle_type = atom->angle_type;
  tagint **angle_atom1 = atom->angle_atom1;
  tagint **angle_atom2 = atom->angle_atom2;
  tagint **angle_atom3 = atom->angle_atom3;
  tagint *tag = atom->tag;

  atype = abs(angle_type[i][m]);
  if ((atype == 0) || !angle_flag[atype] || (angle_distance[atype] <= 0.0)) return 0;
  if (angle_atom2[i][m] != tag[i]) return 0;    // only process at the center atom

  const int a0 = atom->map(angle_atom1[i][m]);
  const int c0 = atom->map(angle_atom3[i][m]);
  if ((a0 < 0) || (c0 < 0)) return 0;
  a = domain->closest_image(i, a0);
  c = domain->closest_image(i, c0);

  if (!(mask[a] & groupbit) || !(mask[c] & groupbit)) return 0;

  const int tab = find_bond_type(i, a);
  const int tbc = find_bond_type(i, c);
  if ((tab == 0) || (tbc == 0)) return 0;
  if (!bond_selected(i, a, tab) || !bond_selected(i, c, tbc)) return 0;

  return 1;
}

/* ----------------------------------------------------------------------
   negate (sign<0) or restore (sign>0) the angle_type of selected angles.
------------------------------------------------------------------------- */

void FixIlves::negate_angle_types(int sign)
{
  int n = atom->nlocal;
  int *num_angle = atom->num_angle;
  int **angle_type = atom->angle_type;

  if (!num_angle) return;

  type = atom->type;    // angle_selected() / bond_selected() read these
  mask = atom->mask;
  mass = atom->mass;
  rmass = atom->rmass;

  for (int i = 0; i < n; ++i) {
    if (!(mask[i] & groupbit)) continue;
    for (int m = 0; m < num_angle[i]; ++m) {
      int a, c, atype;
      if (!angle_selected(i, m, a, c, atype)) continue;
      // a near-linear angle handled by linearangle skip keeps its bonded-style
      // term (it is not constrained), so leave its type sign alone
      if (angle_linear[atype] && (linear_mode == LINEAR_SKIP)) continue;
      if ((sign < 0) && (angle_type[i][m] > 0))
        angle_type[i][m] = -angle_type[i][m];
      else if ((sign > 0) && (angle_type[i][m] < 0))
        angle_type[i][m] = -angle_type[i][m];
    }
  }
}

/* ----------------------------------------------------------------------
   a bond between local i and (local/ghost) j of (positive) type btype is
   selected when its type, either atom type, or either atom mass matches.
------------------------------------------------------------------------- */

int FixIlves::bond_selected(int i, int j, int btype)
{
  if (bond_flag[btype]) return 1;
  if (type_flag[type[i]] || type_flag[type[j]]) return 1;
  if (!mass_list.empty()) {
    const double mi = rmass ? rmass[i] : mass[type[i]];
    const double mj = rmass ? rmass[j] : mass[type[j]];
    if (masscheck(mi) || masscheck(mj)) return 1;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

int FixIlves::masscheck(double massone)
{
  for (const double mv : mass_list)
    if (fabs(mv - massone) <= MASSDELTA) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   number of degrees of freedom removed by the constraints for atoms in igroup.
   each distance constraint removes one DOF; count each once (clist_a is local,
   so a constraint is owned by exactly one rank).
------------------------------------------------------------------------- */

bigint FixIlves::dof(int igroup)
{
  const int igroupbit = group->bitmask[igroup];
  int *amask = atom->mask;

  bigint n = 0;
  for (int k = 0; k < nconstraints; ++k)
    if ((amask[clist_a[k]] & igroupbit) && (amask[clist_b[k]] & igroupbit)) ++n;

  bigint nall = 0;
  MPI_Allreduce(&n, &nall, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  return nall;
}

/* ----------------------------------------------------------------------
   global potential energy of the near-linear-angle restrain substitute,
   summed over all ranks (each restrained angle is counted once at its local
   center).  zero unless linearangle restrain is active.
------------------------------------------------------------------------- */

double FixIlves::compute_scalar()
{
  double all = 0.0;
  MPI_Allreduce(&erestraint, &all, 1, MPI_DOUBLE, MPI_SUM, world);
  return all;
}

/* ----------------------------------------------------------------------
   print per-bond-type constraint statistics: count, average length, and the
   spread (max - min) of the constrained bond lengths across all ranks.
------------------------------------------------------------------------- */

void FixIlves::stats()
{
  const int nb = atom->nbondtypes + 1;
  const int na = atom->nangletypes + 1;
  std::vector<bigint> bcount(nb, 0), bcount_all(nb, 0), acount(na, 0), acount_all(na, 0);
  std::vector<double> bsum(nb, 0.0), bmin(nb, 1.0e20), bmax(nb, 0.0);
  std::vector<double> bsum_all(nb, 0.0), bmin_all(nb, 0.0), bmax_all(nb, 0.0);
  std::vector<double> asum(na, 0.0), amin(na, 1.0e20), amax(na, 0.0);
  std::vector<double> asum_all(na, 0.0), amin_all(na, 0.0), amax_all(na, 0.0);

  double **xx = atom->x;
  for (int k = 0; k < nconstraints; ++k) {
    const int a = clist_a[k], b = clist_b[k], t = clist_btype[k];
    const double dx0 = xx[b][0] - xx[a][0];
    const double dy0 = xx[b][1] - xx[a][1];
    const double dz0 = xx[b][2] - xx[a][2];
    const double r = sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0);
    if (t > 0) {    // bond constraint
      ++bcount[t];
      bsum[t] += r;
      if (r < bmin[t]) bmin[t] = r;
      if (r > bmax[t]) bmax[t] = r;
    } else {    // A-C angle constraint, type encoded as -t
      const int at = -t;
      ++acount[at];
      asum[at] += r;
      if (r < amin[at]) amin[at] = r;
      if (r > amax[at]) amax[at] = r;
    }
  }

  MPI_Allreduce(bcount.data(), bcount_all.data(), nb, MPI_LMP_BIGINT, MPI_SUM, world);
  MPI_Allreduce(bsum.data(), bsum_all.data(), nb, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(bmin.data(), bmin_all.data(), nb, MPI_DOUBLE, MPI_MIN, world);
  MPI_Allreduce(bmax.data(), bmax_all.data(), nb, MPI_DOUBLE, MPI_MAX, world);
  MPI_Allreduce(acount.data(), acount_all.data(), na, MPI_LMP_BIGINT, MPI_SUM, world);
  MPI_Allreduce(asum.data(), asum_all.data(), na, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(amin.data(), amin_all.data(), na, MPI_DOUBLE, MPI_MIN, world);
  MPI_Allreduce(amax.data(), amax_all.data(), na, MPI_DOUBLE, MPI_MAX, world);

  if (comm->me == 0) {
    utils::logmesg(lmp, "Fix ilves constraint statistics at step {} (max {} Newton iterations):\n",
                   update->ntimestep, niter_max);
    for (int t = 1; t < nb; ++t) {
      if (bcount_all[t] == 0) continue;
      utils::logmesg(lmp, "  bond type {}: count {}  ave {:.6g}  spread {:.3g}\n", t, bcount_all[t],
                     bsum_all[t] / (double) bcount_all[t], bmax_all[t] - bmin_all[t]);
    }
    for (int t = 1; t < na; ++t) {
      if (acount_all[t] == 0) continue;
      utils::logmesg(lmp, "  angle type {} (A-C): count {}  ave {:.6g}  spread {:.3g}\n", t,
                     acount_all[t], asum_all[t] / (double) acount_all[t], amax_all[t] - amin_all[t]);
    }
  }

  // reset the iteration counter for the next reporting interval
  niter_max = 0;
}
