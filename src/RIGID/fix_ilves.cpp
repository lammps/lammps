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

#include "atom.h"
#include "bond.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "ilves_asym.h"
#include "ilves_sym.h"
#include "memory.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

namespace {
enum { ILVES_FAST, ILVES_FULL };
constexpr double MASSDELTA = 0.1;

// keywords that terminate the b/a/t/m selector lists
int is_keyword(const char *s)
{
  return ((strcmp(s, "variant") == 0) || (strcmp(s, "store") == 0) ||
          (strcmp(s, "kbond") == 0) || (strcmp(s, "linearangle") == 0));
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
    variant(ILVES_FAST), molecular(0), types_negated(0), nconstraints(0), ilves_solver(nullptr),
    xpred(nullptr), maxatom(0), dtv(0.0), dtfsq(0.0), inv_dtfsq(0.0), x(nullptr), v(nullptr),
    f(nullptr), mass(nullptr), rmass(nullptr), type(nullptr), mask(nullptr), nlocal(0)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_ilves);

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

  // angle constraints are not yet implemented in this phase

  for (int i = 1; i <= atom->nangletypes; ++i)
    if (angle_flag[i])
      error->all(FLERR, "Fix ilves angle constraints are not yet implemented");

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
    } else {
      error->all(FLERR, "Unknown fix ilves keyword: {}", arg[iarg]);
    }
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
  // restore the bond types we negated so the bonded styles act on them again

  if (types_negated && atom->bond_type) negate_bond_types(1);

  delete ilves_solver;
  memory->destroy(xpred);
}

/* ---------------------------------------------------------------------- */

int FixIlves::setmask()
{
  int mask = 0;
  mask |= PRE_NEIGHBOR;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixIlves::init()
{
  if (!force->bond)
    error->all(FLERR, "Fix ilves requires a bond style to define equilibrium bond lengths");

  // equilibrium bond lengths per bond type, as in fix shake

  bond_distance.assign(atom->nbondtypes + 1, 0.0);
  for (int i = 1; i <= atom->nbondtypes; ++i)
    bond_distance[i] = force->bond->equilibrium_distance(i);

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
    types_negated = 1;
  }

  build_constraint_list();
}

/* ---------------------------------------------------------------------- */

void FixIlves::setup(int vflag)
{
  bigint nc = nconstraints, nctot = 0;
  MPI_Allreduce(&nc, &nctot, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (comm->me == 0)
    utils::logmesg(lmp, "Fix ilves: constraining {} bond(s)\n", nctot);

  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixIlves::pre_neighbor()
{
  build_constraint_list();
}

/* ---------------------------------------------------------------------- */

void FixIlves::post_force(int /*vflag*/)
{
  if (!ilves_solver) return;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  // (re)allocate the predicted-position buffer

  if (nlocal > maxatom) {
    memory->destroy(xpred);
    maxatom = atom->nmax;
    memory->create(xpred, maxatom, 3, "ilves:xpred");
  }

  // predict the unconstrained positions, as in fix shake unconstrained_update

  for (int i = 0; i < nlocal; ++i) {
    const double m = rmass ? rmass[i] : mass[type[i]];
    if (m <= 0.0) {
      xpred[i][0] = x[i][0];
      xpred[i][1] = x[i][1];
      xpred[i][2] = x[i][2];
      continue;
    }
    const double dtfm = dtfsq / m;
    xpred[i][0] = x[i][0] + dtv * v[i][0] + dtfm * f[i][0];
    xpred[i][1] = x[i][1] + dtv * v[i][1] + dtfm * f[i][1];
    xpred[i][2] = x[i][2] + dtv * v[i][2] + dtfm * f[i][2];
  }

  // solve the constraints (iterates xpred) and add the constraint forces to f.
  // when the solver takes no iterations the predicted positions already satisfy
  // the constraints, so no force is needed and the (stale) Lagrange multipliers
  // must not be applied -- mirror the numit > 0 guard of the reference code.

  auto result = ilves_solver->solve(x, xpred, tolerance, max_iter);
  if (result.second > 0) ilves_solver->add_constraint_forces(f, inv_dtfsq);
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

  if (num_bond) {
    for (int i = 0; i < nlocal; ++i) {
      if (!(mask[i] & groupbit)) continue;
      for (int m = 0; m < num_bond[i]; ++m) {
        const int btype = abs(bond_type[i][m]);
        if (btype == 0) continue;
        const int j = atom->map(bond_atom[i][m]);
        if (j < 0)
          error->one(FLERR, "Fix ilves bond atom missing on this processor");
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

  nconstraints = (int) clist_a.size();

  // per-atom inverse mass (1/m) for the constrained atoms, handed to the solver

  invmass.assign(nlocal, 0.0);
  for (int k = 0; k < nconstraints; ++k) {
    const int a = clist_a[k];
    const int b = clist_b[k];
    invmass[a] = rmass ? 1.0 / rmass[a] : 1.0 / mass[type[a]];
    invmass[b] = rmass ? 1.0 / rmass[b] : 1.0 / mass[type[b]];
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
