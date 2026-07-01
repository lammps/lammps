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
#include "fix_respa.h"
#include "force.h"
#include "group.h"
#include "ilves.h"
#include "label_map.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

namespace {
enum { LINEAR_ERROR, LINEAR_SKIP, LINEAR_RESTRAIN };
constexpr double MASSDELTA = 0.1;
// Default harmonic-restraint force constant.  It is expressed as a multiple of
// the Boltzmann constant so the default scales with the unit system's energy
// scale instead of being tied to one unit choice.  KBOND_KB_MULT * boltz
// reproduces the value (~1000) found empirically for real units in the
// MD/restrain regime; the minimization regime has no timestep stability limit
// and is KBOND_MIN_FACTOR times stiffer, so its default works out to the
// 1e9*boltz used by fix shake.  An explicit kbond keyword sets the MD value and
// is scaled by the same factor for minimization.
constexpr double KBOND_KB_MULT = 5.0e5;
constexpr double KBOND_MIN_FACTOR = 2.0e3;

// keywords that terminate the b/a/t/m selector lists
int is_keyword(const char *s)
{
  return ((strcmp(s, "store") == 0) || (strcmp(s, "mode") == 0) || (strcmp(s, "kbond") == 0) ||
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
    fixed_iter_flag(false), linear_mode(LINEAR_ERROR), linear_threshold(175.0), kbond(-1.0),
    molecular(0), types_negated(0), store_flag(0), fstore(nullptr), maxstore(0), niter_max(0),
    nconstraints(0), ilves_solver(nullptr), xpred(nullptr), xpred0(nullptr), dx(nullptr),
    maxatom(0), commstage(0), dtv(0.0), dtfsq(0.0), inv_dtfsq(0.0), respa(0), nlevels_respa(0),
    loop_respa(nullptr), step_respa(nullptr), fix_respa(nullptr), dtf_inner(0.0),
    dtf_innerhalf(0.0), x(nullptr), v(nullptr), f(nullptr), mass(nullptr), rmass(nullptr),
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

  // the constraint forces contribute to the global pressure virial, on by
  // default as for fix shake.  the contribution is computed from the converged
  // Lagrange multipliers and the start-of-step bond vectors, so it is
  // reproducible after a restart and does not depend on whether the masses are
  // per-type or per-atom.  use fix_modify virial no to exclude it.
  virial_global_flag = 1;
  thermo_virial = 1;
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
    error->all(FLERR, Error::COMMAND,
               "Fix ilves requires a molecular atom style with bond topology");
  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, Error::COMMAND, "Fix ilves requires an atom map, see atom_modify");

  tolerance = utils::numeric(FLERR, arg[3], false, lmp);
  max_iter = utils::inumeric(FLERR, arg[4], false, lmp);
  output_every = utils::inumeric(FLERR, arg[5], false, lmp);
  if (tolerance <= 0.0) error->all(FLERR, 3, "Fix ilves tolerance must be > 0");
  if (max_iter <= 0) error->all(FLERR, 4, "Fix ilves iteration count must be > 0");

  // size selector membership tables (1-based type indexing)

  bond_flag.assign(atom->nbondtypes + 1, 0);
  angle_flag.assign(atom->nangletypes + 1, 0);
  type_flag.assign(atom->ntypes + 1, 0);

  // allow the b/a/t selectors to take symbolic type labels (as fix shake does),
  // unless a selector keyword (b/a/t/m) collides with an actual type label
  bool allow_typelabels = (atom->labelmapflag != 0);
  if (allow_typelabels) {
    for (int i = Atom::ATOM; i < Atom::DIHEDRAL; ++i) {
      if ((atom->lmap->find_type("b", i) >= 0) || (atom->lmap->find_type("a", i) >= 0) ||
          (atom->lmap->find_type("t", i) >= 0) || (atom->lmap->find_type("m", i) >= 0))
        allow_typelabels = false;
    }
    if (!allow_typelabels && (comm->me == 0))
      error->warning(FLERR,
                     "At least one typelabel conflicts with a fix ilves selector: "
                     "support for typelabels is disabled");
  }

  // parse one or more b/a/t/m selector lists, then optional keyword/value pairs

  int iarg = 6;
  while ((iarg < narg) && !is_keyword(arg[iarg])) {
    if (!is_selector(arg[iarg]))
      error->all(FLERR, iarg, "Unknown fix ilves selector or keyword: {}", arg[iarg]);
    const char sel = arg[iarg][0];
    ++iarg;
    int nvalues = 0;
    while ((iarg < narg) && !is_selector(arg[iarg]) && !is_keyword(arg[iarg])) {
      if (sel == 'm') {
        mass_list.push_back(utils::numeric(FLERR, arg[iarg], false, lmp));
      } else {
        // type may be given as an integer or, with a labelmap, a symbolic label
        const int kind = (sel == 'b') ? Atom::BOND : (sel == 'a') ? Atom::ANGLE : Atom::ATOM;
        const int v = allow_typelabels ? utils::expand_type_int(FLERR, arg[iarg], kind, lmp)
                                       : utils::inumeric(FLERR, arg[iarg], false, lmp);
        if (sel == 'b') {
          if ((v < 1) || (v > atom->nbondtypes))
            error->all(FLERR, iarg, "Invalid fix ilves bond type {}", arg[iarg]);
          bond_flag[v] = 1;
        } else if (sel == 'a') {
          if ((v < 1) || (v > atom->nangletypes))
            error->all(FLERR, iarg, "Invalid fix ilves angle type {}", arg[iarg]);
          angle_flag[v] = 1;
        } else if (sel == 't') {
          if ((v < 1) || (v > atom->ntypes))
            error->all(FLERR, iarg, "Invalid fix ilves atom type {}", arg[iarg]);
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
    if (strcmp(arg[iarg], "store") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix ilves store", error);
      store_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "mode") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix ilves mode", error);
      if (strcmp(arg[iarg + 1], "converge") == 0)
        fixed_iter_flag = false;
      else if (strcmp(arg[iarg + 1], "fixed") == 0)
        fixed_iter_flag = true;
      else
        error->all(FLERR, iarg + 1, "Unknown fix ilves mode: {}", arg[iarg + 1]);
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
        error->all(FLERR, iarg + 1, "Unknown fix ilves linearangle mode: {}", arg[iarg + 1]);
      linear_threshold = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      if ((linear_threshold <= 0.0) || (linear_threshold >= 180.0))
        error->all(FLERR, iarg + 2,
                   "Fix ilves linearangle threshold must be between 0 and 180 degrees");
      iarg += 3;
    } else if (strcmp(arg[iarg], "kbond") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix ilves kbond", error);
      kbond = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (kbond <= 0.0) error->all(FLERR, iarg + 1, "Fix ilves kbond must be > 0");
      iarg += 2;
    } else {
      error->all(FLERR, iarg, "Unknown fix ilves keyword: {}", arg[iarg]);
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
  if (any == 0) error->all(FLERR, Error::COMMAND, "Fix ilves requires at least one b/t/m selector");

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
  mask |= POST_FORCE_RESPA;
  mask |= MIN_PRE_NEIGHBOR;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixIlves::init()
{
  if (modify->get_fix_by_style("^ilves").size() > 1)
    error->all(FLERR,"More than one fix ilves instance");

  if (!force->bond)
    error->all(FLERR, Error::NOLASTLINE,
               "Fix ilves requires a bond style to define equilibrium bond lengths");

  // detect r-RESPA and cache its level structure (as in fix shake).  the
  // per-level timestep factors and the FixRespa pointer are set in setup().
  respa = 0;
  if (utils::strmatch(update->integrate_style, "^respa")) {
    auto *respa_ptr = dynamic_cast<Respa *>(update->integrate);
    if (!respa_ptr)
      error->all(FLERR, Error::NOLASTLINE, "Failure to access Respa style {}",
                 update->integrate_style);
    respa = 1;
    nlevels_respa = respa_ptr->nlevels;
    loop_respa = respa_ptr->loop;
    step_respa = respa_ptr->step;
  }

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
      error->all(FLERR, Error::NOLASTLINE, "Fix ilves angle constraints require an angle style");

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
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix ilves angle type(s){} have an equilibrium angle at or above the "
                   "linearangle threshold of {} degrees and cannot be rigidly constrained; "
                   "use the linearangle keyword to skip or restrain them",
                   linear_types, linear_threshold);
      else if (comm->me == 0)
        error->warning(FLERR,
                       "Fix ilves treating near-linear angle type(s){} with linearangle "
                       "mode {}",
                       linear_types, (linear_mode == LINEAR_SKIP) ? "skip" : "restrain");
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
  // the constraint multipliers to forces, identical to fix shake (SHAKE form).
  // for r-RESPA these are level-dependent and set in setup()/post_force_respa.

  if (!respa) {
    dtv = update->dt;
    dtfsq = update->dt * update->dt * force->ftm2v;
    inv_dtfsq = 1.0 / dtfsq;
  }
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
    if (clist_btype[k] > 0)
      ++nb;
    else
      ++na;
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

  // remove the velocity component along each bond so the run starts with
  // constraint-consistent velocities (as fix shake correct_velocities does)
  project_velocities();

  // precompute the constraint forces for the first integration step
  if (!respa) {
    post_force(vflag);
  } else {
    // find the FixRespa that stores the per-level forces, and set the SHAKE-form
    // (full-step) timestep factors used to predict the unconstrained positions
    if (update->whichflag > 0) {
      auto fixes = modify->get_fix_by_style("^RESPA");
      if (fixes.size() > 0)
        fix_respa = dynamic_cast<FixRespa *>(fixes.front());
      else
        error->all(FLERR, Error::NOLASTLINE, "Run style respa did not create fix RESPA");
    }
    dtf_innerhalf = 0.5 * step_respa[0] * force->ftm2v;
    dtf_inner = step_respa[0] * force->ftm2v;

    // step-0 virial: one Verlet-style solve on the total force (atom->f still
    // holds the summed per-level forces at this point), so the reported initial
    // pressure equals the Verlet value.  the per-level precompute below would
    // instead count the initial constraint reaction once per level.
    dtv = update->dt;
    dtfsq = update->dt * update->dt * force->ftm2v;
    inv_dtfsq = 1.0 / dtfsq;
    post_force(vflag);

    // precompute the per-level constraint forces for the first step (so the
    // first r-RESPA half-kick is constrained), swapping each level's stored
    // force into atom->f around the call (as fix shake does).  the virial is
    // disabled here so it is not changed from the Verlet-style value above.
    dtv = step_respa[0];
    const int saved_vflag_global = vflag_global;
    vflag_global = 0;
    auto *respa_ptr = dynamic_cast<Respa *>(update->integrate);
    for (int ilevel = 0; ilevel < nlevels_respa; ++ilevel) {
      respa_ptr->copy_flevel_f(ilevel);
      post_force_respa(0, ilevel, loop_respa[ilevel] - 1);
      respa_ptr->copy_f_flevel(ilevel);
    }
    vflag_global = saved_vflag_global;
  }
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
    if (clist_btype[k] > 0)
      ++nb;
    else
      ++na;
  bigint nc[2] = {nb, na}, nctot[2] = {0, 0};
  MPI_Allreduce(nc, nctot, 2, MPI_LMP_BIGINT, MPI_SUM, world);
  if (comm->me == 0)
    utils::logmesg(lmp,
                   "Fix ilves: replacing {} bond and {} angle constraint(s) with harmonic "
                   "restraints for minimization\n",
                   nctot[0], nctot[1]);
  min_post_force(vflag);
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

  // solve the constraint system on the prepared prediction and add the forces

  const int numit = solve_constraints();

  // add this step's constraint contribution to the global pressure virial

  if (numit > 0 && vflag_global) ilves_solver->add_global_virial(virial, inv_dtfsq);
}

/* ----------------------------------------------------------------------
   Run the global Newton iteration that drives the predicted positions xpred
   onto the constraint manifold, and return the number of iterations taken.
   xpred must already be prepared by the caller; on return it holds the
   constrained positions (home + ghost).  solve_constraints then turns the
   net displacement into a force.

   The loop is driven uniformly on all ranks: the convergence test uses the
   all-reduced maximum residual, so every rank takes the same number of
   iterations and joins the same collective communication even if it owns no
   constraints.  Cross-rank coupling is a block-Jacobi sweep: each iteration
   reverse-sums the per-atom position increments to their owners, applies them,
   and forward-communicates the predicted positions back to the ghosts.
------------------------------------------------------------------------- */

int FixIlves::run_newton()
{
  const int nall = nlocal + atom->nghost;

  // zero the per-iteration increment buffer (home + ghost)

  for (int i = 0; i < nall; ++i) dx[i][0] = dx[i][1] = dx[i][2] = 0.0;

  commstage = 0;    // forward-comm predicted positions (with PBC shift)
  comm->forward_comm(this);
  double local = ilves_solver ? ilves_solver->prepare(x, xpred) : 0.0;
  double ptau = 0.0;
  if (!fixed_iter_flag) MPI_Allreduce(&local, &ptau, 1, MPI_DOUBLE, MPI_MAX, world);

  int numit = 0;
  for (int i = 0; i < max_iter; ++i) {
    // in convergence mode (the default) stop once the global maximum relative
    // violation is below the tolerance; in fixed mode always run max_iter steps
    // (which avoids the per-iteration MPI reduction)
    if (!fixed_iter_flag && (!std::isfinite(ptau) || (ptau <= tolerance))) break;
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
    if (!fixed_iter_flag) MPI_Allreduce(&local, &ptau, 1, MPI_DOUBLE, MPI_MAX, world);
  }

  if (numit > niter_max) niter_max = numit;
  return numit;
}

/* ----------------------------------------------------------------------
   Drive the predicted positions xpred onto the constraint manifold (relative to
   xpred0, the saved unconstrained prediction) with the Newton iteration, then
   convert the net constrained displacement of each owned atom into a force on
   atom->f (and fstore when store yes), exactly as the multiplier-to-force
   coupling of fix shake (f += m*dx/dtfsq), using the current inv_dtfsq.  Returns
   the number of Newton iterations taken.  Shared by post_force and
   post_force_respa; the caller sets the timestep factors and prepares xpred.
------------------------------------------------------------------------- */

int FixIlves::solve_constraints()
{
  const int numit = run_newton();

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
  }

  return numit;
}

/* ----------------------------------------------------------------------
   enforce the constraints from within r-RESPA.  Like post_force but the
   unconstrained-position prediction and the multiplier-to-force conversion use
   the level-dependent effective timestep, exactly as fix shake does:
     xpred = x + dt0*v + (dt0*dtN/m) fN
                       + sum_{j<N} (1/2 dt0*dtj/m) f_level[j]
   with dt0 = step_respa[0] (innermost) and dtN = step_respa[ilevel].
------------------------------------------------------------------------- */

void FixIlves::post_force_respa(int vflag, int ilevel, int iloop)
{
  x = atom->x;
  v = atom->v;
  f = atom->f;
  type = atom->type;
  mass = atom->mass;
  rmass = atom->rmass;
  nlocal = atom->nlocal;

  // statistics output only on the outermost level

  if (output_every && (ilevel == nlevels_respa - 1) && (update->ntimestep == next_output)) {
    stats();
    next_output += output_every;
  }

  // effective timestep for this level (SHAKE/full-step form)

  dtfsq = dtf_inner * step_respa[ilevel];
  inv_dtfsq = 1.0 / dtfsq;

  // the global virial accumulates the per-level contributions: zero it at the
  // innermost level's last sub-iteration (as fix shake does), then add each
  // level's contribution at that level's last sub-iteration

  const bool last_iloop = (iloop == loop_respa[ilevel] - 1);
  if ((ilevel == 0) && last_iloop && vflag) v_init(vflag);

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

  // the near-linear-angle restraint is a force-field-like force; add it once
  // per step at the innermost level only, so it is not counted at every level

  if ((linear_mode == LINEAR_RESTRAIN) && (ilevel == 0)) apply_linear_restraint();

  // predict the unconstrained positions, including this and all inner levels

  double ***f_level = fix_respa->f_level;
  for (int i = 0; i < nlocal; ++i) {
    const double m = rmass ? rmass[i] : mass[type[i]];
    if (m > 0.0) {
      const double invm = 1.0 / m;
      const double dtfm = dtfsq * invm;
      xpred[i][0] = x[i][0] + dtv * v[i][0] + dtfm * f[i][0];
      xpred[i][1] = x[i][1] + dtv * v[i][1] + dtfm * f[i][1];
      xpred[i][2] = x[i][2] + dtv * v[i][2] + dtfm * f[i][2];
      for (int j = 0; j < ilevel; ++j) {
        const double c = dtf_innerhalf * step_respa[j] * invm;
        xpred[i][0] += c * f_level[i][j][0];
        xpred[i][1] += c * f_level[i][j][1];
        xpred[i][2] += c * f_level[i][j][2];
      }
    } else {
      xpred[i][0] = x[i][0];
      xpred[i][1] = x[i][1];
      xpred[i][2] = x[i][2];
    }
    xpred0[i][0] = xpred[i][0];
    xpred0[i][1] = xpred[i][1];
    xpred0[i][2] = xpred[i][2];
  }

  const int numit = solve_constraints();

  if (numit > 0 && last_iloop && vflag_global) ilves_solver->add_global_virial(virial, inv_dtfsq);
}

/* ----------------------------------------------------------------------
   Project the velocities onto the constraint manifold once at the start of a
   run (remove the relative velocity along each constrained bond), the analogue
   of fix shake correct_velocities, so the run begins constraint-consistent.
   This is not repeated every step: the position-constraint force already keeps
   the bonds rigid, exactly as fix shake (which likewise does not project
   velocities during the run).  The velocity constraint is linear, so this is a
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
  clist_node_a.clear();
  clist_node_b.clear();
  clist_d.clear();
  clist_btype.clear();
  clist_vertex.clear();
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
          error->one(FLERR,
                     "Fix ilves bond atom missing on this processor; increase the "
                     "communication cutoff with comm_modify cutoff");
        // GEOMETRY index: the nearest periodic image, so the raw bond vector in
        // the solver is correct at any box size (no minimum-image restriction).
        j = domain->closest_image(i, j);
        // with newton_bond off the bond is stored on both atoms; keep one copy
        if (!newton_bond && (tag[i] > tag[j])) continue;
        if (!(mask[j] & groupbit)) continue;
        if (!bond_selected(i, j, btype)) continue;
        // NODE id: the canonical owner of the (possibly ghost) geometry atom, so
        // an atom and its periodic image share one node and a wrapped bond stays
        // a single graph edge.  atom i is local, hence already its own owner.  A
        // genuinely off-rank partner has no local index, so map() returns its
        // ghost -- the intended cross-rank block-Jacobi halo node.
        clist_a.push_back(i);
        clist_b.push_back(j);
        clist_node_a.push_back(i);
        clist_node_b.push_back(atom->map(atom->tag[j]));
        clist_d.push_back(bond_distance[btype]);
        clist_btype.push_back(btype);
        clist_vertex.push_back(-1);    // bonds have no angle vertex
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
        // canonical node ids of the two outer atoms (owner of each geometry
        // image), so the A-C virtual bond joins the same graph nodes as its two
        // flanking bonds even when the angle spans a periodic boundary.
        clist_node_a.push_back(atom->map(atom->tag[a]));
        clist_node_b.push_back(atom->map(atom->tag[c]));
        clist_d.push_back(angle_distance[atype]);
        clist_btype.push_back(-atype);    // negative marks an A-C angle constraint
        clist_vertex.push_back(i);        // center atom, for reporting the angle
      }
    }
  }

  nconstraints = (int) clist_a.size();

  // per-atom inverse mass (1/m) for the constrained atoms, handed to the solver.
  // sized for local+ghost atoms since a constraint partner may be a ghost.

  const int nall = atom->nlocal + atom->nghost;
  invmass.assign(nall, 0.0);
  for (int i = 0; i < nlocal; ++i) invmass[i] = rmass ? 1.0 / rmass[i] : 1.0 / mass[type[i]];

  // communicate the inverse mass to the ghosts, so a constraint partner owned by
  // another rank (or a periodic image) has the correct value regardless of
  // whether the atom style communicates per-atom mass to ghosts.

  if (atom->nghost > 0) {
    commstage = 2;
    comm->forward_comm(this);
  }

  // (re)build the ILVES solver for the current constraint topology

  delete ilves_solver;
  ilves_solver = nullptr;
  if (nconstraints > 0)
    ilves_solver = new ILVES::Ilves(lmp, clist_a, clist_b, clist_node_a, clist_node_b, clist_d,
                                    invmass);
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
  const double k = (kbond > 0.0) ? kbond : KBOND_KB_MULT * force->boltz;
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

  // base (MD-regime) force constant, then stiffen the true constraints by the
  // minimization factor; the near-linear-angle A-C restraints stay at the base
  // value since they are an inherently soft substitute in either regime
  const double k_md = (kbond > 0.0) ? kbond : KBOND_KB_MULT * force->boltz;
  const double k_bond = k_md * KBOND_MIN_FACTOR;

  for (int i = 0; i < nall; ++i) dx[i][0] = dx[i][1] = dx[i][2] = 0.0;
  erestraint = 0.0;

  const int nc = nconstraints;
  for (int kk = 0; kk < nc; ++kk)
    erestraint += min_harmonic_bond(clist_a[kk], clist_b[kk], clist_d[kk], k_bond);
  const int nr = (int) rlist_a.size();
  for (int kk = 0; kk < nr; ++kk)
    erestraint += min_harmonic_bond(rlist_a[kk], rlist_c[kk], rlist_d[kk], k_md);

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
  dx[a][0] += fac * ux;
  dx[a][1] += fac * uy;
  dx[a][2] += fac * uz;
  dx[b][0] -= fac * ux;
  dx[b][1] -= fac * uy;
  dx[b][2] -= fac * uz;
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

  const tagint ti = tag[i];
  const tagint tj = tag[j];

  // the bond joining i and j is stored (newton_bond on) at only one endpoint's
  // LOCAL copy; look it up at the canonical (owner) index of each atom so a
  // partner reached as a periodic ghost still resolves a bond stored at its
  // owner.  atom->map(tag) is the owner's local index, or a ghost only when the
  // atom is genuinely off-rank.
  const int ic = atom->map(ti);
  const int jc = atom->map(tj);

  if ((ic >= 0) && (ic < n)) {
    for (int m = 0; m < num_bond[ic]; ++m)
      if (bond_atom[ic][m] == tj) return abs(bond_type[ic][m]);
  }
  if ((jc >= 0) && (jc < n)) {
    for (int m = 0; m < num_bond[jc]; ++m)
      if (bond_atom[jc][m] == ti) return abs(bond_type[jc][m]);
  }
  return 0;
}

/* ----------------------------------------------------------------------
   decide whether angle m of (local center) atom i is constrained.  it is when
   its type is selected, all three atoms are in the group, both flanking bonds
   are themselves selected (constrained), and the A-C distance is known.  fills
   the closest-image (geometry) outer-atom indices a, c and the angle type.
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
  // GEOMETRY indices: nearest periodic image of the two outer atoms, so the A-C
  // bond vector is a correct raw subtraction at any box size.  The caller records
  // the canonical node ids separately for graph connectivity.
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
   estimate the memory used by this fix (per MPI rank).  the term that grows
   with the size of the connected constraint clusters is the ported solver's
   factored sparse matrix, reported via ilves_solver->memory_usage().
------------------------------------------------------------------------- */

double FixIlves::memory_usage()
{
  double bytes = 0.0;

  // per-atom buffers grown to nmax: xpred, xpred0, dx (plus fstore with store yes)
  bytes += 3.0 * (double) maxatom * 3 * sizeof(double);
  if (store_flag) bytes += (double) maxstore * 3 * sizeof(double);
  bytes += (double) invmass.size() * sizeof(double);

  // local constraint and near-linear-restraint lists, rebuilt every reneighbor
  bytes += (double) (clist_a.size() + clist_b.size() + clist_btype.size() + clist_vertex.size() +
                     rlist_a.size() + rlist_c.size()) *
      sizeof(int);
  bytes += (double) (clist_d.size() + rlist_d.size()) * sizeof(double);

  // per-type selector and equilibrium-distance lookup tables
  bytes +=
      (double) (bond_flag.size() + angle_flag.size() + type_flag.size() + angle_linear.size()) *
      sizeof(int);
  bytes +=
      (double) (bond_distance.size() + angle_distance.size() + mass_list.size()) * sizeof(double);

  // the ported ILVES solver (topology graphs + factored sparse matrix)
  if (ilves_solver) bytes += ilves_solver->memory_usage();

  return bytes;
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

  // tally bond lengths and (like fix shake) the actual bend angle in degrees of
  // each constrained angle, computed from its two legs at the center atom.

  double **xx = atom->x;
  for (int k = 0; k < nconstraints; ++k) {
    const int a = clist_a[k], b = clist_b[k], t = clist_btype[k];
    if (t > 0) {    // bond constraint: the bond length
      const double dx0 = xx[b][0] - xx[a][0];
      const double dy0 = xx[b][1] - xx[a][1];
      const double dz0 = xx[b][2] - xx[a][2];
      const double r = sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0);
      ++bcount[t];
      bsum[t] += r;
      if (r < bmin[t]) bmin[t] = r;
      if (r > bmax[t]) bmax[t] = r;
    } else {    // A-C angle constraint (type -t): the A-vertex-C angle in degrees
      const int at = -t;
      const int vtx = clist_vertex[k];
      const double r1x = xx[a][0] - xx[vtx][0], r1y = xx[a][1] - xx[vtx][1],
                   r1z = xx[a][2] - xx[vtx][2];
      const double r2x = xx[b][0] - xx[vtx][0], r2y = xx[b][1] - xx[vtx][1],
                   r2z = xx[b][2] - xx[vtx][2];
      const double r1 = sqrt(r1x * r1x + r1y * r1y + r1z * r1z);
      const double r2 = sqrt(r2x * r2x + r2y * r2y + r2z * r2z);
      double cosv = (r1x * r2x + r1y * r2y + r1z * r2z) / (r1 * r2);
      if (cosv > 1.0) cosv = 1.0;
      if (cosv < -1.0) cosv = -1.0;
      const double angle = acos(cosv) * 180.0 / MathConst::MY_PI;
      ++acount[at];
      asum[at] += angle;
      if (angle < amin[at]) amin[at] = angle;
      if (angle > amax[at]) amax[at] = angle;
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

  // print in the same type/ave/delta/count layout as fix shake, so fix ilves is
  // a drop-in for tools that read SHAKE statistics; the header also reports the
  // largest Newton iteration count used since the last output (ILVES-specific).

  if (comm->me == 0) {
    // when a labelmap is present, report each constrained type by its symbolic
    // label (matching the label-based input), with a numeric fallback for any
    // type that has no label; otherwise report numeric types as before.
    const bool uselabel = (atom->labelmapflag != 0);
    auto blabel = [&](int t) -> std::string {
      if (uselabel) {
        const std::string &s = atom->lmap->find_label(t, Atom::BOND);
        if (!s.empty()) return s;
      }
      return std::to_string(t);
    };
    auto alabel = [&](int t) -> std::string {
      if (uselabel) {
        const std::string &s = atom->lmap->find_label(t, Atom::ANGLE);
        if (!s.empty()) return s;
      }
      return std::to_string(t);
    };

    // type-column width: longest label (or type number) actually printed,
    // shared by the Bond: and Angle: rows so the columns stay aligned
    int width;
    if (uselabel) {
      width = 1;
      for (int t = 1; t < nb; ++t)
        if (bcount_all[t]) {
          const int w = (int) blabel(t).size();
          if (w > width) width = w;
        }
      for (int t = 1; t < na; ++t)
        if (acount_all[t]) {
          const int w = (int) alabel(t).size();
          if (w > width) width = w;
        }
    } else {
      int maxt = (nb > na) ? nb : na;
      if (maxt < 1) maxt = 1;
      width = (int) log10((double) maxt) + 2;
    }

    auto mesg = fmt::format("ILVES stats (type/ave/delta/count) on step {} ({}{} Newton "
                            "iterations)\n",
                            update->ntimestep, (fixed_iter_flag ? "" : "up to "), niter_max);
    for (int t = 1; t < nb; ++t) {
      if (bcount_all[t] == 0) continue;
      if (uselabel) {
        mesg += fmt::format("Bond:  {:<{}}   {:<9.6} {:<11.6} {:>8d}\n", blabel(t), width,
                            bsum_all[t] / (double) bcount_all[t], bmax_all[t] - bmin_all[t],
                            bcount_all[t]);
      } else {
        mesg += fmt::format("Bond:  {:>{}}   {:<9.6} {:<11.6} {:>8d}\n", t, width,
                            bsum_all[t] / (double) bcount_all[t], bmax_all[t] - bmin_all[t],
                            bcount_all[t]);
      }
    }
    for (int t = 1; t < na; ++t) {
      if (acount_all[t] == 0) continue;
      if (uselabel) {
        mesg += fmt::format("Angle: {:<{}}   {:<9.6} {:<11.6} {:>8d}\n", alabel(t), width,
                            asum_all[t] / (double) acount_all[t], amax_all[t] - amin_all[t],
                            acount_all[t]);
      } else {
        mesg += fmt::format("Angle: {:>{}}   {:<9.6} {:<11.6} {:>8d}\n", t, width,
                            asum_all[t] / (double) acount_all[t], amax_all[t] - amin_all[t],
                            acount_all[t]);
      }
    }

    utils::logmesg(lmp, mesg);
  }

  // reset the iteration counter for the next reporting interval
  niter_max = 0;
}
