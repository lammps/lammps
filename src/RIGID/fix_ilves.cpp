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
   ILVES bond/angle constraint solver

   Reference:
     P. Garcia-Risueno et al., J. Chem. Theory Comput. (2026)
     DOI: 10.1021/acs.jctc.5c01376
     reference GROMACS implementation:
     https://github.com/LorienLV/_PAPER_ILVES

   Follows the user-interface of fix shake using the same b/a/t/m selectors
   (constrain by bond type, angle type, atom type, atom mass) and the same
   in-group requirement: a bond is constrained only when both atoms are
   in the fix group; an angle (i.e. the A-C "virtual" distance derived
   from A-B and B-C bonds) is constrained only when all three atoms are
   in the fix group AND the two flanking bonds are themselves selected by
   one of the above criteria.

   Unlike SHAKE, the constraint topology is not limited to small clusters
   (1+1, 1+2, 1+3 atoms or one 3-atom angle).  Constraints are stored as a
   flat list of (a, b) pairs with target distance d_k, so connected
   constraint networks of arbitrary size (e.g. all C-C backbone bonds in a
   protein) are admissible.

------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
   Created using Claude Code with Opus 4.7
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
#include "label_map.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "respa.h"
#include "update.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

namespace {
enum { ILVES_FULL, ILVES_FAST };
constexpr double BIG = 1.0e20;
constexpr double MASSDELTA = 0.1;
constexpr int DELTA_CONSTR = 256;

const char cite_fix_ilves[] =
    "fix ilves command: https://doi.org/10.1021/acs.jctc.5c01376\n\n"
    "@Article{LopezVillellas2025,\n"
    "author = {L{\'o}pez-Villellas,  Lori{\'e}n and Mikkelsen,  Carl Christian Kjelgaard "
    "and Galano-Frutos,  Juan Jos{\'e} and Marco-Sola,  Santiago and "
    "Alastruey-Bened{\'e},  Jes{\'u}s and Ib{\'a}{\~n}ez,  Pablo and "
    "Echenique,  Pablo and Moret{\'o},  Miquel and {De Rosa},  Maria Cristina and "
    "Garc{\'i}a-Risue{\~n}o,  Pablo},\n"
    "title = {ILVES: Accurate and Efficient Bond Length and Angle Constraints in Molecular "
    "Dynamics},\n"
    "volume = {21},\n"
    "ISSN = {1549-9626},\n"
    "url = {http://dx.doi.org/10.1021/acs.jctc.5c01376},\n"
    "DOI = {10.1021/acs.jctc.5c01376},\n"
    "number = {18},\n"
    "journal = {Journal of Chemical Theory and Computation},\n"
    "publisher = {American Chemical Society (ACS)},\n"
    "year = {2025},\n"
    "month = sep,\n"
    "pages = {8711-8719}\n"
    "}\n\n";
}    // namespace

/* ---------------------------------------------------------------------- */

FixIlves::FixIlves(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), bond_flag(nullptr), angle_flag(nullptr), type_flag(nullptr),
    mass_list(nullptr), bond_distance(nullptr), angle_distance(nullptr), angle_r1(nullptr),
    angle_r2(nullptr), angle_linear(nullptr), fstore(nullptr),
    ilves_flag(nullptr), xshake(nullptr), c_atom1(nullptr), c_atom2(nullptr),
    c_type(nullptr), c_dist(nullptr), c_lambda(nullptr), c_cluster(nullptr), c_rx(nullptr),
    c_ry(nullptr), c_rz(nullptr), c_rsq(nullptr), c_invma(nullptr), c_invmb(nullptr),
    cluster_offset(nullptr), c_perm(nullptr), c_slot(nullptr), lu_A(nullptr),
    lu_b(nullptr), lu_pivot(nullptr), cl_sx(nullptr), cl_sy(nullptr), cl_sz(nullptr),
    chol_pool(nullptr), chol_pool_offset(nullptr), cluster_bw(nullptr), cluster_cached(nullptr),
    x(nullptr), v(nullptr), f(nullptr), mass(nullptr), rmass(nullptr), type(nullptr),
    b_count(nullptr), b_count_all(nullptr), b_ave(nullptr), b_max(nullptr), b_min(nullptr),
    b_ave_all(nullptr), b_max_all(nullptr), b_min_all(nullptr), a_count(nullptr),
    a_count_all(nullptr), a_ave(nullptr), a_max(nullptr), a_min(nullptr), a_ave_all(nullptr),
    a_max_all(nullptr), a_min_all(nullptr)
{
  lu_alloc = 0;
  largest_cluster = 0;
  comm_mode = 0;
  global_topology_ready = false;
  chol_pool_alloc = 0;
  chol_offset_alloc = 0;
  cluster_bw_alloc = 0;
  cluster_cached_alloc = 0;
  largest_bw = 0;
  energy_global_flag = energy_peratom_flag = 1;
  virial_global_flag = virial_peratom_flag = 1;
  thermo_energy = thermo_virial = 1;
  create_attribute = 1;
  dof_flag = 1;
  scalar_flag = 1;
  extscalar = 1;
  stores_ids = 1;
  next_output = -1;

  vflag_post_force = 0;
  eflag_pre_reverse = 0;
  ebond = 0.0;
  has_angle = false;
  variant = ILVES_FAST;
  chol_calls = 0;
  chol_fallbacks = 0;
  linear_threshold = 165.0;

  store_flag = peratom_flag = 0;
  maxstore = -1;

  n_constr = 0;
  max_constr = 0;
  n_clusters = 0;

  respa = 0;
  nlevels_respa = 0;
  loop_respa = nullptr;
  step_respa = nullptr;
  fix_respa = nullptr;
  dtf_inner = 0.0;

  if (lmp->citeme) lmp->citeme->add(cite_fix_ilves);

  // error checks

  if (atom->molecular != Atom::MOLECULAR)
    error->all(FLERR, Error::COMMAND,
               "Fix ilves requires a molecular system (atom_style with bonds)");

  // perform initial allocation of atom-based per-atom arrays
  FixIlves::grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);

  // forward-comm payload: 3 doubles (xshake)
  comm_forward = 3;

  // parse arguments, same as for fix shake

  if (narg < 8) utils::missing_cmd_args(FLERR, "fix ilves", error);

  // clang-format off

  tolerance    = utils::numeric(FLERR,  arg[3], false, lmp);
  max_iter     = utils::inumeric(FLERR, arg[4], false, lmp);
  output_every = utils::inumeric(FLERR, arg[5], false, lmp);

  if (tolerance <= 0.0) error->all(FLERR, 3, "Fix ilves tolerance must be > 0");
  if (max_iter <= 0)    error->all(FLERR, 4, "Fix ilves max_iter must be > 0");
  if (output_every < 0) error->all(FLERR, 5, "Fix ilves output_every must be >= 0");

  // check typelabel/keyword conflicts (b/a/t/m vs typelabels)
  bool allow_typelabels = (atom->labelmapflag != 0);
  if (allow_typelabels) {
    for (int i = Atom::ATOM; i < Atom::DIHEDRAL; ++i) {
      if ((atom->lmap->find_type("b", i) >= 0) ||
          (atom->lmap->find_type("a", i) >= 0) ||
          (atom->lmap->find_type("t", i) >= 0) ||
          (atom->lmap->find_type("m", i) >= 0)) allow_typelabels = false;
    }
    if (!allow_typelabels && (comm->me == 0))
      error->warning(FLERR, "At least one typelabel conflicts with a fix ilves option: "
                     "support for typelabels is disabled.");
  }

  bond_flag  = new int[atom->nbondtypes  + 1]{};
  angle_flag = new int[atom->nangletypes + 1]{};
  type_flag  = new int[atom->ntypes      + 1]{};
  mass_list  = new double[atom->ntypes]{};
  nmass = 0;

  char mode = '\0';
  int next = 6;
  while (next < narg) {
    int i = -1;
    if      (strcmp(arg[next], "b") == 0) mode = 'b';
    else if (strcmp(arg[next], "a") == 0) mode = 'a';
    else if (strcmp(arg[next], "t") == 0) mode = 't';
    else if (strcmp(arg[next], "m") == 0) { mode = 'm'; atom->check_mass(FLERR); }

    // break on known optional keyword
    else if ((strcmp(arg[next], "kbond")       == 0) ||
             (strcmp(arg[next], "store")       == 0) ||
             (strcmp(arg[next], "variant")     == 0) ||
             (strcmp(arg[next], "linearangle") == 0)) {
      break;

    } else if (mode == 'b') {
      if (allow_typelabels) i = utils::expand_type_int(FLERR, arg[next], Atom::BOND, lmp);
      else                  i = utils::inumeric(FLERR, arg[next], false, lmp);
      if (i < 1 || i > atom->nbondtypes)
        error->all(FLERR, next, "Invalid bond type {} for fix ilves", arg[next]);
      bond_flag[i] = 1;

    } else if (mode == 'a') {
      if (allow_typelabels) i = utils::expand_type_int(FLERR, arg[next], Atom::ANGLE, lmp);
      else                  i = utils::inumeric(FLERR, arg[next], false, lmp);
      if (i < 1 || i > atom->nangletypes)
        error->all(FLERR, next, "Invalid angle type {} for fix ilves", arg[next]);
      angle_flag[i] = 1;
      has_angle = true;

    } else if (mode == 't') {
      if (allow_typelabels) i = utils::expand_type_int(FLERR, arg[next], Atom::ATOM, lmp);
      else                  i = utils::inumeric(FLERR, arg[next], false, lmp);
      if (i < 1 || i > atom->ntypes)
        error->all(FLERR, next, "Invalid atom type {} for fix ilves", arg[next]);
      type_flag[i] = 1;

    } else if (mode == 'm') {
      double m1 = utils::numeric(FLERR, arg[next], false, lmp);
      if (m1 <= 0.0) error->all(FLERR, next, "Invalid atom mass {} for fix ilves", arg[next]);
      if (nmass == atom->ntypes)
        error->all(FLERR, "Too many mass entries for fix ilves");
      mass_list[nmass++] = m1;

    } else {
      error->all(FLERR, "Unknown fix ilves command option: {}", arg[next]);
    }
    next++;
  }

  // optional args: kbond <v>, store <yes|no>, variant <full|fast>,
  //                linearangle <theta_deg> [Lmin]
  kbond = 1.0e9 * force->boltz;

  while (next < narg) {
    if (strcmp(arg[next], "kbond") == 0) {
      if (next + 2 > narg) utils::missing_cmd_args(FLERR, "fix ilves kbond", error);
      kbond = utils::numeric(FLERR, arg[next + 1], false, lmp);
      if (kbond <= 0.0) error->all(FLERR, next + 1, "Fix ilves kbond must be > 0");
      next += 2;

    } else if (strcmp(arg[next], "store") == 0) {
      if (next + 2 > narg) utils::missing_cmd_args(FLERR, "fix ilves store", error);
      store_flag = utils::logical(FLERR, arg[next + 1], false, lmp);
      if (store_flag) {
        peratom_flag = 1;
        size_peratom_cols = 3;
        peratom_freq = 1;
      }
      next += 2;

    } else if (strcmp(arg[next], "variant") == 0) {
      if (next + 2 > narg) utils::missing_cmd_args(FLERR, "fix ilves variant", error);
      if      (strcmp(arg[next + 1], "full") == 0) variant = ILVES_FULL;
      else if (strcmp(arg[next + 1], "fast") == 0) variant = ILVES_FAST;
      else error->all(FLERR, next + 1, "Unknown fix ilves variant {}", arg[next + 1]);
      next += 2;

    } else if (strcmp(arg[next], "linearangle") == 0) {
      if (next + 2 > narg) utils::missing_cmd_args(FLERR, "fix ilves linearangle", error);
      linear_threshold = utils::numeric(FLERR, arg[next + 1], false, lmp);
      if ((linear_threshold < 150.0) || (linear_threshold > 180.0))
        error->all(FLERR, next + 1,
                   "Fix ilves linearangle must be in the [150, 180] degrees range, got {}",
                   linear_threshold);
      next += 2;

    } else {
      error->all(FLERR, next, "Unknown fix ilves command option: {}", arg[next]);
    }
  }

  // allocate distance arrays (filled by init())
  bond_distance  = new double[atom->nbondtypes  + 1]{};
  angle_distance = new double[atom->nangletypes + 1]{};
  angle_r1       = new double[atom->nangletypes + 1]{};
  angle_r2       = new double[atom->nangletypes + 1]{};
  // Near-linear angle flag (filled by init_topology after
  // equilibrium_angle/_distance are available).  When set, the AC virtual
  // constraint is NOT added for that angle type; the standard angle
  // force-field term handles the angle.
  angle_linear   = new int[atom->nangletypes + 1]{};

  if (output_every) {
    const int nb = atom->nbondtypes  + 1;
    const int na = atom->nangletypes + 1;
    b_count = new bigint[nb]{};  b_count_all = new bigint[nb]{};
    b_ave   = new double[nb]{};  b_ave_all   = new double[nb]{};
    b_max   = new double[nb]{};  b_max_all   = new double[nb]{};
    b_min   = new double[nb];    b_min_all   = new double[nb];
    for (int i = 0; i < nb; ++i) b_min[i] = b_min_all[i] = BIG;
    a_count = new bigint[na]{};  a_count_all = new bigint[na]{};
    a_ave   = new double[na]{};  a_ave_all   = new double[na]{};
    a_max   = new double[na]{};  a_max_all   = new double[na]{};
    a_min   = new double[na];   a_min_all   = new double[na];
    for (int i = 0; i < na; ++i) a_min[i] = a_min_all[i] = BIG;
  }
}

/* ---------------------------------------------------------------------- */

FixIlves::~FixIlves()
{
  if (modify->get_fix_by_id(id)) atom->delete_callback(id, Atom::GROW);

  delete[] bond_flag;
  delete[] angle_flag;
  delete[] type_flag;
  delete[] mass_list;
  delete[] bond_distance;
  delete[] angle_distance;
  delete[] angle_r1;
  delete[] angle_r2;
  delete[] angle_linear;

  delete[] b_count;     delete[] b_count_all;
  delete[] b_ave;       delete[] b_ave_all;
  delete[] b_max;       delete[] b_max_all;
  delete[] b_min;       delete[] b_min_all;
  delete[] a_count;     delete[] a_count_all;
  delete[] a_ave;       delete[] a_ave_all;
  delete[] a_max;       delete[] a_max_all;
  delete[] a_min;       delete[] a_min_all;

  memory->destroy(ilves_flag);
  memory->destroy(xshake);
  memory->destroy(fstore);

  memory->destroy(c_atom1);
  memory->destroy(c_atom2);
  memory->destroy(c_type);
  memory->destroy(c_dist);
  memory->destroy(c_lambda);
  memory->destroy(c_cluster);
  memory->destroy(c_rx);
  memory->destroy(c_ry);
  memory->destroy(c_rz);
  memory->destroy(c_rsq);
  memory->destroy(c_invma);
  memory->destroy(c_invmb);
  memory->destroy(cluster_offset);
  memory->destroy(c_perm);
  memory->destroy(c_slot);
  memory->destroy(lu_A);
  memory->destroy(lu_b);
  memory->destroy(lu_pivot);
  memory->destroy(cl_sx);
  memory->destroy(cl_sy);
  memory->destroy(cl_sz);
  memory->destroy(chol_pool);
  memory->destroy(chol_pool_offset);
  memory->destroy(cluster_bw);
  memory->destroy(cluster_cached);
}

/* ---------------------------------------------------------------------- */

int FixIlves::setmask()
{
  int mask = 0;
  mask |= PRE_NEIGHBOR;
  mask |= POST_NEIGHBOR;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= END_OF_STEP;
  mask |= MIN_PRE_REVERSE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ----------------------------------------------------------------------
   set bond/angle equilibrium distances from the active bond/angle styles
------------------------------------------------------------------------- */

void FixIlves::init()
{
  // only one fix ilves instance allowed
  if (modify->get_fix_by_style("^ilves").size() > 1)
    error->all(FLERR, Error::NOLASTLINE, "More than one fix ilves instance");

  // detect rRESPA and stash the per-level step/loop arrays and the
  // FixRespa pointer (auto-created by run_style respa to hold f_level).
  // No special handling beyond that: the multi-level path is engaged via
  // post_force_respa() / unconstrained_update_respa().
  respa = 0;
  fix_respa = nullptr;
  if (utils::strmatch(update->integrate_style, "^respa")) {
    if (update->whichflag > 0) {
      auto fixes = modify->get_fix_by_style("^RESPA");
      if (fixes.size() > 0) fix_respa = dynamic_cast<FixRespa *>(fixes.front());
      else error->all(FLERR, Error::NOLASTLINE,
                      "Run style respa did not create fix RESPA");
    }
    auto *respa_ptr = dynamic_cast<Respa *>(update->integrate);
    if (!respa_ptr)
      error->all(FLERR, Error::NOLASTLINE,
                 "Failure to access Respa style {}", update->integrate_style);
    respa = 1;
    nlevels_respa = respa_ptr->nlevels;
    loop_respa = respa_ptr->loop;
    step_respa = respa_ptr->step;
  }

  // forbid box-changing fixes between fix ilves and integration?  we don't
  // care about the order strictly, but fix shake does this check and we
  // mirror it for consistency.
  bool boxflag = false;
  for (const auto &ifix : modify->get_fix_list()) {
    if (boxflag && utils::strmatch(ifix->style, "^ilves"))
      error->all(FLERR, Error::NOLASTLINE, "Fix ilves must come before any box-changing fix");
    if (ifix->box_change) boxflag = true;
  }

  // need a bond style
  if (force->bond == nullptr)
    error->all(FLERR, Error::NOLASTLINE, "Bond style must be defined for fix ilves");

  for (int i = 1; i <= atom->nbondtypes; ++i)
    bond_distance[i] = force->bond->equilibrium_distance(i);

  // sanity check: if angle constraints are requested, an angle_style must
  // be defined (the variant-specific init_topology() below uses
  // force->angle->equilibrium_angle to compute angle_distance[]).
  if (has_angle && force->angle == nullptr)
    error->all(FLERR, Error::NOLASTLINE,
               "Angle style must be defined for fix ilves with angle constraints");

  // Variant-specific topology setup: the local variant validates that
  // every constraint cluster fits within this rank's local-plus-ghost
  // reach; the global variant does an MPI_Allgatherv over all bonds and
  // angles.  Both fill angle_distance[] for the angle types selected.
  init_topology();

  // Identify near-linear angle types using the equilibrium geometry now
  // available in angle_r1/angle_r2/angle_distance.  The angle theta_0
  // is recovered via the law of cosines from r1, r2, d:
  //     cos(theta) = (r1^2 + r2^2 - d^2) / (2 r1 r2)
  // If the result corresponds to theta >= linear_threshold the angle's
  // AC virtual constraint is ill-conditioned (triangle inequality
  // saturates) and would risk excessive constraint forces.  Mark the
  // type as near-linear; build_constraint_list will skip the AC entry
  // and leave the angle_type slot unnegated so the standard angle
  // force-field term continues to act on it.  If equilibrium_angle was
  // not provided by the angle style, theta_0 defaults to 0 and the
  // derived theta is 0 -- correctly leaving the angle type non-linear.
  if (has_angle && linear_threshold < 180.0) {
    const double cos_threshold = cos(linear_threshold * MY_PI / 180.0);
    int n_linear = 0;
    for (int at = 1; at <= atom->nangletypes; ++at) {
      angle_linear[at] = 0;
      if (!angle_flag[at]) continue;
      const double r1 = angle_r1[at];
      const double r2 = angle_r2[at];
      const double d  = angle_distance[at];
      if (r1 <= 0.0 || r2 <= 0.0 || d <= 0.0) continue;
      double cos_theta = (r1*r1 + r2*r2 - d*d) / (2.0 * r1 * r2);
      if (cos_theta < -1.0) cos_theta = -1.0;
      if (cos_theta >  1.0) cos_theta =  1.0;
      if (cos_theta > cos_threshold) continue;     // theta < threshold
      angle_linear[at] = 1;
      ++n_linear;
    }
    if (n_linear && comm->me == 0)
      utils::logmesg(lmp,
                     "Fix ilves: skipping AC virtual constraint for {} near-linear angle "
                     "type(s) (theta_0 >= {} deg); the angle force-field term remains active "
                     "for these angles\n", n_linear, linear_threshold);
  } else {
    for (int at = 1; at <= atom->nangletypes; ++at) angle_linear[at] = 0;
  }

  // warn (don't error) on minimization use, since we will substitute strong
  // harmonic restraints in min_post_force.
  if ((comm->me == 0) && (update->whichflag == 2))
    error->warning(FLERR,
                   "Using fix ilves with minimization: substituting constraints with harmonic "
                   "restraint forces (kbond = {:.6g})", kbond);

  // negate bond_type / angle_type for every interaction that we are going
  // to constrain.  This makes the corresponding bond/angle style skip the
  // potential and force contribution for those interactions (it tests for
  // bondtype <= 0 / angletype <= 0).  Without this, the bonded forces and
  // our constraint forces would double up, injecting energy each step.
  negate_constrained_topology();
}

/* ----------------------------------------------------------------------
   For every bond and angle that our selectors will eventually constrain,
   flip the stored bond_type / angle_type to its negative.  Bond / angle
   styles in LAMMPS skip interactions with non-positive type, which is
   exactly the behaviour we want for constrained interactions.

   Idempotent: only flips a positive type, never a negative one.

   This is called once from init() before any time stepping, mirroring
   fix shake's find_clusters / bondtype_findset pattern.
------------------------------------------------------------------------- */

void FixIlves::negate_constrained_topology()
{
  int nlocal_now    = atom->nlocal;
  int **nb_type     = atom->bond_type;
  tagint **nb_atom  = atom->bond_atom;
  int *num_bond     = atom->num_bond;
  int **na_type     = atom->angle_type;
  tagint **na_atom1 = atom->angle_atom1;
  tagint **na_atom2 = atom->angle_atom2;
  tagint **na_atom3 = atom->angle_atom3;
  int *num_angle    = atom->num_angle;
  int *mask         = atom->mask;
  tagint *tag       = atom->tag;
  int *type         = atom->type;
  double *mass      = atom->mass;
  double *rmass     = atom->rmass;

  // bond pass
  for (int i = 0; i < nlocal_now; ++i) {
    if (!(mask[i] & groupbit)) continue;
    for (int b = 0; b < num_bond[i]; ++b) {
      int bt = nb_type[i][b];
      if (bt <= 0) continue;       // already negated (or disabled)
      if (bt > atom->nbondtypes) continue;
      tagint tj = nb_atom[i][b];
      int j = atom->map(tj);
      if (j < 0) continue;
      if (!(mask[j] & groupbit)) continue;

      bool selected = false;
      if (bond_flag[bt]) selected = true;
      else if (type_flag[type[i]] || type_flag[type[j]]) selected = true;
      else if (nmass) {
        double mi = rmass ? rmass[i] : mass[type[i]];
        double mj = rmass ? rmass[j] : mass[type[j]];
        if (masscheck(mi) || masscheck(mj)) selected = true;
      }
      if (!selected) continue;

      // negate the storage on this rank (the partner rank handles its own
      // copy when newton_bond=0; with newton_bond=on only one rank holds
      // the bond, and the negate happens there).
      nb_type[i][b] = -bt;
    }
  }

  // angle pass: negate the angle type when both flanking bonds are
  // themselves constrained (matching the build_constraint_list filter)
  if (has_angle) {
    for (int i = 0; i < nlocal_now; ++i) {
      if (!(mask[i] & groupbit)) continue;
      for (int m = 0; m < num_angle[i]; ++m) {
        int at = na_type[i][m];
        if (at <= 0) continue;
        if (at > atom->nangletypes) continue;
        if (!angle_flag[at]) continue;
        if (angle_linear[at]) continue;        // no AC constraint -- keep
                                               // the angle force-field term
        if (na_atom2[i][m] != tag[i]) continue;

        tagint t1 = na_atom1[i][m];
        tagint t3 = na_atom3[i][m];
        int ia1 = atom->map(t1);
        int ia3 = atom->map(t3);
        if (ia1 < 0 || ia3 < 0) continue;
        if (!(mask[ia1] & groupbit) || !(mask[ia3] & groupbit)) continue;

        // both flanking bonds (i-ia1 and i-ia3) must themselves be
        // constrained.  The local variant of bond_is_constrained()
        // walks local bond storage; the global variant consults the
        // replicated bond table -- so this works regardless of which
        // rank holds the bond under newton_bond=on.
        if (!bond_is_constrained(tag[i], t1)) continue;
        if (!bond_is_constrained(tag[i], t3)) continue;

        na_type[i][m] = -at;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixIlves::setup(int vflag)
{
  // refresh pointers, build initial constraint list
  pre_neighbor();
  post_neighbor();

  int n_info[2] = {n_constr, n_clusters};
  int n_info_all[2];
  int n_max_all;
  MPI_Reduce(n_info, n_info_all, 2, MPI_INT, MPI_SUM, 0, world);
  n_info_all[0] /= comm->nprocs;
  n_info_all[1] /= comm->nprocs;
  MPI_Reduce(&largest_cluster, &n_max_all, 1, MPI_INT, MPI_MAX, 0, world);

  if (comm->me == 0)
    utils::logmesg(lmp, "Fix ilves info:\n"
                   "  {} algorithm\n"
                   "  {} constraints/proc\n"
                   "  {} clusters/proc\n"
                   "  {} max. constraints/cluster\n",
                   (variant == ILVES_FULL) ? "ILVES (full, LU)" : "ILVES (fast, banded Cholesky)",
                   n_info_all[0], n_info_all[1], n_max_all);

  // At setup we need a different dtfsq convention than during normal stepping:
  //
  //   normal post_force(t): force change propagates through TWO half-steps
  //     (final_integrate at t  +  initial_integrate at t+dt), so position
  //     correction = dt^2/m * f_c -> dtfsq = dt^2.
  //
  //   setup post_force: force change propagates through only ONE half-step
  //     (the very first initial_integrate), so position correction is
  //     0.5*dt^2/m * f_c -> dtfsq = 0.5*dt^2.
  //
  // We mirror fix shake's correct_coordinates / correct_velocities /
  // shake_end_of_step pattern: project x onto the constraint surface
  // directly, remove velocity components along constraints, and add
  // constraint forces for the first integration step.
  //
  // For rRESPA, dtv = step_respa[0] (the outer step) and the constraint
  // forces must be applied at every level so that the per-level forces
  // saved in fix_respa->f_level reflect the constrained dynamics.

  if (!respa) {
    dtv = update->dt;
    dtfsq = 0.5 * update->dt * update->dt * force->ftm2v;
    correct_coordinates(vflag);
    post_force(vflag);
    dtfsq = update->dt * update->dt * force->ftm2v;
  } else {
    dtv = step_respa[0];
    dtf_inner = 0.5 * step_respa[0] * force->ftm2v;
    dtfsq = 0.5 * step_respa[0] * step_respa[0] * force->ftm2v;
    correct_coordinates(vflag);
    // apply constraint forces to every rRESPA level, mirroring
    // FixShake::shake_end_of_step.  copy_flevel_f swaps atom->f with
    // f_level[ilevel] so post_force_respa sees the per-level forces;
    // copy_f_flevel swaps them back, capturing our added constraint
    // contributions into the per-level store.
    auto *respa_ptr = dynamic_cast<Respa *>(update->integrate);
    for (int ilevel = 0; ilevel < nlevels_respa; ++ilevel) {
      respa_ptr->copy_flevel_f(ilevel);
      post_force_respa(vflag, ilevel, loop_respa[ilevel] - 1);
      respa_ptr->copy_f_flevel(ilevel);
    }
    dtf_inner = step_respa[0] * force->ftm2v;
  }

  if (output_every) {
    bigint nt = update->ntimestep;
    next_output = nt + output_every;
    if (nt % output_every != 0)
      next_output = (nt/output_every) * output_every + output_every;
    stats();
  } else next_output = -1;
}

/* ---------------------------------------------------------------------- */

void FixIlves::min_setup(int vflag)
{
  pre_neighbor();
  post_neighbor();
  min_post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixIlves::pre_neighbor()
{
  // refresh atom-class pointers; needed because reneighbor can change them
  x = atom->x;
  v = atom->v;
  f = atom->f;
  mass  = atom->mass;
  rmass = atom->rmass;
  type  = atom->type;
  nlocal = atom->nlocal;

  // Catch silent topology changes (fix bond/create, fix bond/break,
  // delete_atoms, create_atoms, set type ...) that would invalidate the
  // gathered global tables.  Skip on the first call (before init_topology
  // has run; global_topology_ready is the flag set at end of init).
  if (global_topology_ready) check_topology_unchanged();
}

/* ---------------------------------------------------------------------- */

void FixIlves::post_neighbor()
{
  // rebuild the flat constraint list after each reneighbor
  build_constraint_list();
}

/* ----------------------------------------------------------------------
   Helper: is a bond between atoms with tags ta and tb selected by the
   current b/t/m criteria?  Uses local atom data; both atoms must be at
   least locally accessible (own copy or ghost).

   Reads atom->type / atom->mass / atom->rmass directly rather than the
   cached member pointers because this is also called from init() (via
   negate_constrained_topology -> bond_is_constrained) where the member
   pointers have not yet been refreshed.
------------------------------------------------------------------------- */


bool FixIlves::bond_selected_for_atoms(int ia, int ib, int bt)
{
  if (bond_flag[bt]) return true;
  int *atype     = atom->type;
  double *amass  = atom->mass;
  double *armass = atom->rmass;
  if (type_flag[atype[ia]] || type_flag[atype[ib]]) return true;
  if (nmass) {
    double mia = armass ? armass[ia] : amass[atype[ia]];
    double mib = armass ? armass[ib] : amass[atype[ib]];
    if (masscheck(mia) || masscheck(mib)) return true;
  }
  return false;
}


/* ----------------------------------------------------------------------
   Union-find over the constraint graph, then sort constraints so that
   constraints belonging to the same cluster are contiguous in c_perm/c_slot.
   For cluster-ID purposes we canonicalize ghost periodic images via
   atom->map(atom->tag[idx]) so that the same physical atom always maps
   to the same index in the union-find.
------------------------------------------------------------------------- */

void FixIlves::group_by_cluster()
{
  n_clusters = 0;

  if (n_constr == 0) {
    if (!cluster_offset) memory->create(cluster_offset, 1, "ilves:cluster_offset");
    cluster_offset[0] = 0;
    return;
  }

  const int nmax = atom->nmax;
  int *parent = nullptr;
  memory->create(parent, nmax, "ilves:uf_parent");
  for (int i = 0; i < nmax; ++i) parent[i] = i;

  auto find_root = [&](int a) {
    while (parent[a] != a) { parent[a] = parent[parent[a]]; a = parent[a]; }
    return a;
  };

  tagint *tag = atom->tag;

  // canonicalize: collapse periodic-image ghosts of locally-owned atoms to
  // their local index, so that union-find sees the same physical atom as
  // one node.  True MPI ghosts (atoms not owned locally at all) keep their
  // ghost index.
  for (int k = 0; k < n_constr; ++k) {
    int a = c_atom1[k];                          // typically local
    int b_raw = c_atom2[k];
    int b_canon = atom->map(tag[b_raw]);
    if (b_canon < 0) b_canon = b_raw;            // shouldn't happen normally
    int ra = find_root(a);
    int rb = find_root(b_canon);
    if (ra != rb) parent[ra] = rb;
  }

  int *atomlabel = nullptr;
  memory->create(atomlabel, nmax, "ilves:uf_label");
  for (int i = 0; i < nmax; ++i) atomlabel[i] = -1;

  for (int k = 0; k < n_constr; ++k) {
    int r = find_root(c_atom1[k]);
    if (atomlabel[r] < 0) atomlabel[r] = n_clusters++;
    c_cluster[k] = atomlabel[r];
  }

  // build cluster_offset[] via counting sort on c_cluster
  memory->grow(cluster_offset, n_clusters + 1, "ilves:cluster_offset");
  for (int c = 0; c <= n_clusters; ++c) cluster_offset[c] = 0;
  for (int k = 0; k < n_constr; ++k) cluster_offset[c_cluster[k] + 1]++;
  for (int c = 0; c < n_clusters; ++c) cluster_offset[c + 1] += cluster_offset[c];

  // place each constraint into its slot
  int *cursor = nullptr;
  memory->create(cursor, n_clusters, "ilves:cursor");
  for (int c = 0; c < n_clusters; ++c) cursor[c] = cluster_offset[c];
  for (int k = 0; k < n_constr; ++k) {
    int c = c_cluster[k];
    int s = cursor[c]++;
    c_perm[s] = k;
    c_slot[k] = s;
  }

  // track largest cluster size for workspace sizing
  largest_cluster = 0;
  for (int c = 0; c < n_clusters; ++c) {
    int sz = cluster_offset[c + 1] - cluster_offset[c];
    if (sz > largest_cluster) largest_cluster = sz;
  }

  memory->destroy(cursor);
  memory->destroy(parent);
  memory->destroy(atomlabel);
}

/* ----------------------------------------------------------------------
   Precompute per-constraint reference vector r_k = x[a]-x[b] (closest
   periodic image), |r_k|^2, and inverse masses.  These are constant
   throughout the constraint solve (they use x at the start of the step).
------------------------------------------------------------------------- */

void FixIlves::precompute_constraint_data()
{
  for (int k = 0; k < n_constr; ++k) {
    int a = c_atom1[k];
    int b = c_atom2[k];
    // Apply minimum_image to the raw vector r_k = x[a]-x[b].  closest_image
    // can return inconsistent ghost copies for atoms that span PBC with
    // different image flags (e.g. neighbouring atoms in the same molecule
    // stored in different periodic cells), giving a wrong "bond vector".
    // minimum_image directly wraps the displacement to the closest periodic
    // image, independent of which storage index was returned.
    double rx = atom->x[a][0] - atom->x[b][0];
    double ry = atom->x[a][1] - atom->x[b][1];
    double rz = atom->x[a][2] - atom->x[b][2];
    domain->minimum_image(FLERR, rx, ry, rz);
    c_rx[k] = rx;
    c_ry[k] = ry;
    c_rz[k] = rz;
    c_rsq[k] = rx*rx + ry*ry + rz*rz;
    c_invma[k] = rmass ? 1.0 / rmass[a] : 1.0 / mass[type[a]];
    c_invmb[k] = rmass ? 1.0 / rmass[b] : 1.0 / mass[type[b]];
  }

  // Per-cluster RCM reorder + bandwidth compute, then size the banded
  // Cholesky factor cache.  RCM must run before grow_factor_cache because
  // the pool size depends on cluster_bw[].
  rcm_reorder_clusters();
  grow_factor_cache();
}

/* ----------------------------------------------------------------------
   Per-cluster reverse Cuthill-McKee on the constraint-adjacency graph
   (two constraints share a graph edge iff they share an atom).  Applies
   the resulting permutation in-place to c_perm and c_slot, then computes
   each cluster's bandwidth (max over edges (k,l) with k<l of slot_l-slot_k
   in the new ordering) and stores it in cluster_bw[].

   The natural cluster ordering (insertion order from group_by_cluster)
   often has bandwidth close to n_c for large connected clusters such as
   protein backbones.  RCM typically reduces it to O(sqrt(n_c)) or
   smaller, which is what makes banded Cholesky competitive with full
   sparse direct factorization for these particular matrices.

   Cost: O(sum_c (n_c + nnz_c)) per call, where nnz_c is the number of
   shared-atom adjacencies in cluster c.  Runs only at reneighbor steps.
------------------------------------------------------------------------- */

void FixIlves::rcm_reorder_clusters()
{
  if (n_clusters == 0) {
    largest_bw = 0;
    return;
  }

  if (n_clusters > cluster_bw_alloc) {
    cluster_bw_alloc = n_clusters;
    memory->grow(cluster_bw, cluster_bw_alloc, "ilves:cluster_bw");
  }

  // adjacency list in slot space (per cluster, slot indices are 0..n_c-1).
  // Built by grouping constraints under their two atom tags: every pair of
  // constraints sharing an atom yields a (slot_a, slot_b) adjacency.
  std::unordered_map<tagint, std::vector<int>> atom_slots;
  std::vector<std::vector<int>> adj;
  std::vector<int> rcm_order;
  std::vector<int> deg;
  std::vector<int> bfs_queue;
  std::vector<char> visited;
  std::vector<int> level_start;
  std::vector<int> old_to_new;
  std::vector<int> old_perm;

  largest_bw = 0;

  for (int c = 0; c < n_clusters; ++c) {
    const int beg = cluster_offset[c];
    const int end = cluster_offset[c + 1];
    const int n_c = end - beg;
    if (n_c <= 1) {
      cluster_bw[c] = 0;
      continue;
    }

    // collect (atom tag) -> list of slots (relative to beg)
    atom_slots.clear();
    for (int s = 0; s < n_c; ++s) {
      int k = c_perm[beg + s];
      atom_slots[atom->tag[c_atom1[k]]].push_back(s);
      atom_slots[atom->tag[c_atom2[k]]].push_back(s);
    }

    // build adjacency list: for each atom shared by m constraints, the m
    // slots are pairwise adjacent.  Symmetric, deduplicated below.
    adj.assign(n_c, {});
    for (auto &kv : atom_slots) {
      const auto &vec = kv.second;
      for (size_t i = 0; i + 1 < vec.size(); ++i)
        for (size_t j = i + 1; j < vec.size(); ++j) {
          adj[vec[i]].push_back(vec[j]);
          adj[vec[j]].push_back(vec[i]);
        }
    }
    for (int s = 0; s < n_c; ++s) {
      std::sort(adj[s].begin(), adj[s].end());
      adj[s].erase(std::unique(adj[s].begin(), adj[s].end()), adj[s].end());
    }

    // node degrees
    deg.assign(n_c, 0);
    for (int s = 0; s < n_c; ++s) deg[s] = (int) adj[s].size();

    // RCM: BFS starting from a low-degree vertex.  For each level, sort
    // newly-visited neighbours by increasing degree (the standard CM
    // tie-breaker).  Reverse the final order at the end.
    visited.assign(n_c, 0);
    rcm_order.clear();
    rcm_order.reserve(n_c);

    int start = 0;
    for (int s = 1; s < n_c; ++s) if (deg[s] < deg[start]) start = s;

    // A cluster is by definition a connected component, but RCM handles
    // disjoint pieces too via the outer "while unvisited remain" loop.
    while ((int) rcm_order.size() < n_c) {
      if (visited[start]) {
        // find the next unvisited low-degree vertex
        start = -1;
        for (int s = 0; s < n_c; ++s) {
          if (!visited[s] && (start < 0 || deg[s] < deg[start])) start = s;
        }
        if (start < 0) break;
      }
      bfs_queue.clear();
      bfs_queue.push_back(start);
      visited[start] = 1;
      size_t head = 0;
      while (head < bfs_queue.size()) {
        int v = bfs_queue[head++];
        rcm_order.push_back(v);
        // enqueue unvisited neighbours, sorted by degree
        int nbeg = (int) bfs_queue.size();
        for (int u : adj[v]) {
          if (!visited[u]) {
            visited[u] = 1;
            bfs_queue.push_back(u);
          }
        }
        std::sort(bfs_queue.begin() + nbeg, bfs_queue.end(),
                  [&](int a, int b) { return deg[a] < deg[b]; });
      }
    }
    // reverse
    std::reverse(rcm_order.begin(), rcm_order.end());

    // old_to_new[old_slot] = new_slot
    old_to_new.assign(n_c, -1);
    for (int s = 0; s < n_c; ++s) old_to_new[rcm_order[s]] = s;

    // apply permutation to c_perm[beg..end) and update c_slot[]
    old_perm.assign(n_c, 0);
    for (int s = 0; s < n_c; ++s) old_perm[s] = c_perm[beg + s];
    for (int s = 0; s < n_c; ++s) {
      int new_slot_for_old_s = old_to_new[s];
      c_perm[beg + new_slot_for_old_s] = old_perm[s];
    }
    for (int s = 0; s < n_c; ++s) c_slot[c_perm[beg + s]] = beg + s;

    // bandwidth in the new ordering
    int bw_c = 0;
    for (int s_old = 0; s_old < n_c; ++s_old) {
      int s_new = old_to_new[s_old];
      for (int u_old : adj[s_old]) {
        int u_new = old_to_new[u_old];
        int d = u_new - s_new;
        if (d < 0) d = -d;
        if (d > bw_c) bw_c = d;
      }
    }
    cluster_bw[c] = bw_c;
    if (bw_c > largest_bw) largest_bw = bw_c;
  }
}

/* ---------------------------------------------------------------------- */

void FixIlves::add_constraint(int a, int b, int btype, double dist)
{
  if (n_constr == max_constr) grow_constraint_list(max_constr + DELTA_CONSTR);
  c_atom1[n_constr]   = a;
  c_atom2[n_constr]   = b;
  c_type[n_constr]    = btype;
  c_dist[n_constr]    = dist;
  c_lambda[n_constr]  = 0.0;
  c_cluster[n_constr] = -1;
  n_constr++;
}

/* ---------------------------------------------------------------------- */

void FixIlves::grow_constraint_list(int newcap)
{
  max_constr = newcap;
  memory->grow(c_atom1,   max_constr, "ilves:c_atom1");
  memory->grow(c_atom2,   max_constr, "ilves:c_atom2");
  memory->grow(c_type,    max_constr, "ilves:c_type");
  memory->grow(c_dist,    max_constr, "ilves:c_dist");
  memory->grow(c_lambda,  max_constr, "ilves:c_lambda");
  memory->grow(c_cluster, max_constr, "ilves:c_cluster");
  memory->grow(c_rx,      max_constr, "ilves:c_rx");
  memory->grow(c_ry,      max_constr, "ilves:c_ry");
  memory->grow(c_rz,      max_constr, "ilves:c_rz");
  memory->grow(c_rsq,     max_constr, "ilves:c_rsq");
  memory->grow(c_invma,   max_constr, "ilves:c_invma");
  memory->grow(c_invmb,   max_constr, "ilves:c_invmb");
  memory->grow(c_perm,    max_constr, "ilves:c_perm");
  memory->grow(c_slot,    max_constr, "ilves:c_slot");
}

/* ---------------------------------------------------------------------- */

int FixIlves::masscheck(double m1)
{
  for (int i = 0; i < nmass; ++i)
    if (fabs(m1 - mass_list[i]) <= MASSDELTA) return 1;
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixIlves::unconstrained_update()
{
  // do this for all atoms (not just constrained ones) so that ghost
  // forward_comm of xshake delivers correct values even for atoms that
  // are not directly constrained on this rank but are partners of
  // constraints owned by neighbouring ranks.
  if (rmass) {
    for (int i = 0; i < nlocal; ++i) {
      double dtfm = dtfsq / rmass[i];
      xshake[i][0] = x[i][0] + dtv*v[i][0] + dtfm*f[i][0];
      xshake[i][1] = x[i][1] + dtv*v[i][1] + dtfm*f[i][1];
      xshake[i][2] = x[i][2] + dtv*v[i][2] + dtfm*f[i][2];
    }
  } else {
    for (int i = 0; i < nlocal; ++i) {
      double dtfm = dtfsq / mass[type[i]];
      xshake[i][0] = x[i][0] + dtv*v[i][0] + dtfm*f[i][0];
      xshake[i][1] = x[i][1] + dtv*v[i][1] + dtfm*f[i][1];
      xshake[i][2] = x[i][2] + dtv*v[i][2] + dtfm*f[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixIlves::post_force(int vflag)
{
  if (update->ntimestep == next_output) stats();

  // reset Cholesky stats
  chol_calls = 0;
  chol_fallbacks = 0;

  // refresh atom pointers in case they moved
  x = atom->x;
  v = atom->v;
  f = atom->f;
  mass  = atom->mass;
  rmass = atom->rmass;
  type  = atom->type;
  nlocal = atom->nlocal;

  // recompute reference vectors r_k = x[a]-x[b] (closest periodic image)
  for (int k = 0; k < n_constr; ++k) {
    int a = c_atom1[k];
    int b = c_atom2[k];
    double rx = x[a][0] - x[b][0];
    double ry = x[a][1] - x[b][1];
    double rz = x[a][2] - x[b][2];
    domain->minimum_image(FLERR, rx, ry, rz);
    c_rx[k]  = rx;
    c_ry[k]  = ry;
    c_rz[k]  = rz;
    c_rsq[k] = rx*rx + ry*ry + rz*rz;
  }

  // predict unconstrained move xshake, then share with ghosts
  unconstrained_update();
  comm->forward_comm(this);

  ev_init(eflag_pre_reverse, vflag);
  ebond = 0.0;

  // run Newton iteration; converged xshake satisfies all distance constraints
  solve_constraints();

  // convert accumulated Lagrange multipliers into atom forces (so that the
  // next initial_integrate produces the constrained positions)
  apply_constraint_forces(vflag);

  vflag_post_force = vflag;
}

/* ----------------------------------------------------------------------
   rRESPA variant of unconstrained_update.  The unconstrained position
   prediction at innermost-loop endpoint of level ilevel is

     xshake = x + dt0*v + (dtN/m)*f
                        + sum_{j<ilevel} 0.5*dt_j/m * f_level[i][j]

   where dtN is the time step at level ilevel and f is the current force
   (also at level ilevel since respa swaps atom->f with f_level on each
   level entry).  dtfsq is set to dt0 * dtN so the solver's downstream
   scaling produces a per-level Lagrange-force magnitude compatible with
   that level's propagation.  See FixShake::unconstrained_update_respa
   for the original derivation.
------------------------------------------------------------------------- */

void FixIlves::unconstrained_update_respa(int ilevel)
{
  const double dtf_innerhalf = 0.5 * step_respa[0] * force->ftm2v;
  double ***f_level = fix_respa->f_level;
  dtfsq = dtf_inner * step_respa[ilevel];

  if (rmass) {
    for (int i = 0; i < nlocal; ++i) {
      const double invmass = 1.0 / rmass[i];
      double dtfmsq = dtfsq * invmass;
      xshake[i][0] = x[i][0] + dtv*v[i][0] + dtfmsq*f[i][0];
      xshake[i][1] = x[i][1] + dtv*v[i][1] + dtfmsq*f[i][1];
      xshake[i][2] = x[i][2] + dtv*v[i][2] + dtfmsq*f[i][2];
      for (int jlevel = 0; jlevel < ilevel; ++jlevel) {
        dtfmsq = dtf_innerhalf * step_respa[jlevel] * invmass;
        xshake[i][0] += dtfmsq * f_level[i][jlevel][0];
        xshake[i][1] += dtfmsq * f_level[i][jlevel][1];
        xshake[i][2] += dtfmsq * f_level[i][jlevel][2];
      }
    }
  } else {
    for (int i = 0; i < nlocal; ++i) {
      const double invmass = 1.0 / mass[type[i]];
      double dtfmsq = dtfsq * invmass;
      xshake[i][0] = x[i][0] + dtv*v[i][0] + dtfmsq*f[i][0];
      xshake[i][1] = x[i][1] + dtv*v[i][1] + dtfmsq*f[i][1];
      xshake[i][2] = x[i][2] + dtv*v[i][2] + dtfmsq*f[i][2];
      for (int jlevel = 0; jlevel < ilevel; ++jlevel) {
        dtfmsq = dtf_innerhalf * step_respa[jlevel] * invmass;
        xshake[i][0] += dtfmsq * f_level[i][jlevel][0];
        xshake[i][1] += dtfmsq * f_level[i][jlevel][1];
        xshake[i][2] += dtfmsq * f_level[i][jlevel][2];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Enforce ILVES constraints from rRESPA.  Called once per (ilevel,iloop)
   pair by the respa integrator.  The structural skeleton mirrors
   FixShake::post_force_respa: stats() only fires at the outermost level;
   v_init / evflag are activated on the last iteration of each level for
   correct virial accumulation.  Coordinate prediction and constraint
   solving use the per-level dtfsq computed inside
   unconstrained_update_respa().
------------------------------------------------------------------------- */

void FixIlves::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa - 1 && update->ntimestep == next_output) stats();

  chol_calls = 0;
  chol_fallbacks = 0;

  x = atom->x;
  v = atom->v;
  f = atom->f;
  mass  = atom->mass;
  rmass = atom->rmass;
  type  = atom->type;
  nlocal = atom->nlocal;

  // recompute reference vectors r_k = x[a]-x[b] (same as post_force)
  for (int k = 0; k < n_constr; ++k) {
    int a = c_atom1[k];
    int b = c_atom2[k];
    double rx = x[a][0] - x[b][0];
    double ry = x[a][1] - x[b][1];
    double rz = x[a][2] - x[b][2];
    domain->minimum_image(FLERR, rx, ry, rz);
    c_rx[k]  = rx;
    c_ry[k]  = ry;
    c_rz[k]  = rz;
    c_rsq[k] = rx*rx + ry*ry + rz*rz;
  }

  // predict unconstrained move (sets dtfsq for the current level), then
  // share xshake with ghosts.
  unconstrained_update_respa(ilevel);
  comm->forward_comm(this);

  // virial init once per outer step at the innermost level's last loop;
  // evflag = 1 on the last loop of each level (so constraint virial is
  // accumulated once per level into the running totals; using ev_init
  // here would zero the per-atom virial accumulated at other levels).
  if (ilevel == 0 && iloop == loop_respa[ilevel] - 1 && vflag) v_init(vflag);
  if (iloop == loop_respa[ilevel] - 1) evflag = 1;
  else evflag = 0;
  ebond = 0.0;

  solve_constraints();
  apply_constraint_forces(vflag);

  vflag_post_force = vflag;
}

/* ----------------------------------------------------------------------
   Grow the per-cluster dense-LU workspace to handle clusters of size n
------------------------------------------------------------------------- */

void FixIlves::grow_lu_workspace(int n)
{
  if (n <= lu_alloc) return;
  lu_alloc = n;
  memory->grow(lu_A,     lu_alloc * lu_alloc, "ilves:lu_A");
  memory->grow(lu_b,     lu_alloc,            "ilves:lu_b");
  memory->grow(lu_pivot, lu_alloc,            "ilves:lu_pivot");
  memory->grow(cl_sx,    lu_alloc,            "ilves:cl_sx");
  memory->grow(cl_sy,    lu_alloc,            "ilves:cl_sy");
  memory->grow(cl_sz,    lu_alloc,            "ilves:cl_sz");
}

/* ----------------------------------------------------------------------
   In-place LU factorization with partial pivoting + forward/back solve.
   lu_A is row-major n x n; lu_b is the rhs in/out.  Returns 0 on success,
   nonzero if the pivot is below tolerance (matrix is singular).
------------------------------------------------------------------------- */

int FixIlves::lu_factor_solve(int n)
{
  constexpr double TINY = 1.0e-30;

  // factorization with row pivoting
  for (int i = 0; i < n; ++i) lu_pivot[i] = i;

  for (int k = 0; k < n; ++k) {
    // find pivot in column k from row k onward
    double pmax = fabs(lu_A[k*n + k]);
    int pidx = k;
    for (int i = k + 1; i < n; ++i) {
      double v = fabs(lu_A[i*n + k]);
      if (v > pmax) { pmax = v; pidx = i; }
    }
    if (pmax < TINY) return 1;     // singular

    if (pidx != k) {
      // swap rows k and pidx (full rows of A, and pivot record)
      for (int j = 0; j < n; ++j) {
        double t = lu_A[k*n + j];
        lu_A[k*n + j] = lu_A[pidx*n + j];
        lu_A[pidx*n + j] = t;
      }
      int t = lu_pivot[k]; lu_pivot[k] = lu_pivot[pidx]; lu_pivot[pidx] = t;
    }

    // eliminate
    double diag_inv = 1.0 / lu_A[k*n + k];
    for (int i = k + 1; i < n; ++i) {
      double factor = lu_A[i*n + k] * diag_inv;
      lu_A[i*n + k] = factor;       // store L in lower triangle
      for (int j = k + 1; j < n; ++j) {
        lu_A[i*n + j] -= factor * lu_A[k*n + j];
      }
    }
  }

  // permute b according to pivot
  // lu_pivot[k] = original row index now sitting in slot k
  // apply permutation P to b: b_new[k] = b_old[ lu_pivot[k] ]
  // (using cl_sx as scratch since it's at least lu_alloc long)
  for (int k = 0; k < n; ++k) cl_sx[k] = lu_b[lu_pivot[k]];
  for (int k = 0; k < n; ++k) lu_b[k] = cl_sx[k];

  // forward substitution L y = b (L has unit diagonal, stored below diag)
  for (int i = 1; i < n; ++i) {
    double s = lu_b[i];
    for (int j = 0; j < i; ++j) s -= lu_A[i*n + j] * lu_b[j];
    lu_b[i] = s;
  }
  // back substitution U x = y
  for (int i = n - 1; i >= 0; --i) {
    double s = lu_b[i];
    for (int j = i + 1; j < n; ++j) s -= lu_A[i*n + j] * lu_b[j];
    lu_b[i] = s / lu_A[i*n + i];
  }
  return 0;
}

/* ----------------------------------------------------------------------
   In-place Cholesky factorization A = L L^T (lower triangular) followed
   by forward + back substitution to solve A x = b.  Requires A to be
   symmetric positive definite; reads only the lower triangle (i >= j) of
   the assembled matrix and overwrites it with L.  The upper triangle
   (i < j) is left untouched.

   Returns 0 on success, 1 if the running diagonal becomes non-positive
   (matrix is not SPD).  In that case the caller should fall back to
   lu_factor_solve, which can handle indefinite or asymmetric A.

   This is the workhorse for the symmetric ILVES_FAST position Jacobian
   and for the (always-symmetric) velocity-projection matrix.  Asymptotic
   cost is ~n^3/6 flops vs ~n^3/3 for the LU above.
------------------------------------------------------------------------- */

int FixIlves::chol_factor_solve(int n)
{
  int info = chol_factor(n, lu_A);
  if (info) return info;
  chol_solve(n, lu_A, lu_b);
  return 0;
}

/* ----------------------------------------------------------------------
   In-place Cholesky factor of symmetric positive-definite A (row-major,
   n x n).  Reads only the lower triangle (i >= j); the upper triangle is
   ignored.  Writes L into the lower triangle (including the diagonal).
   Returns 0 on success or 1 if a running diagonal becomes non-positive
   (matrix is not SPD).  Cost: ~n^3/6 flops.
------------------------------------------------------------------------- */

int FixIlves::chol_factor(int n, double *A)
{
  constexpr double TINY = 1.0e-30;

  for (int k = 0; k < n; ++k) {
    double diag = A[k*n + k];
    for (int j = 0; j < k; ++j) {
      double Lkj = A[k*n + j];
      diag -= Lkj * Lkj;
    }
    if (diag < TINY) return 1;     // not SPD
    double Lkk = sqrt(diag);
    A[k*n + k] = Lkk;
    double inv_Lkk = 1.0 / Lkk;
    for (int i = k + 1; i < n; ++i) {
      double t = A[i*n + k];
      for (int j = 0; j < k; ++j) t -= A[i*n + j] * A[k*n + j];
      A[i*n + k] = t * inv_Lkk;
    }
  }
  return 0;
}

/* ----------------------------------------------------------------------
   Forward + back substitution to solve L L^T x = b given a factor L that
   was produced by chol_factor.  Operates in-place on b.  L is read from
   the lower triangle of the n x n matrix at L[].  Cost: ~n^2 flops --
   much cheaper than chol_factor's ~n^3/6, which is why caching the
   factor across Newton iterations is worth the extra storage.
------------------------------------------------------------------------- */

void FixIlves::chol_solve(int n, const double *L, double *b)
{
  // Forward substitution: L y = b  (overwrite b with y).
  for (int i = 0; i < n; ++i) {
    double s = b[i];
    for (int j = 0; j < i; ++j) s -= L[i*n + j] * b[j];
    b[i] = s / L[i*n + i];
  }
  // Back substitution: L^T x = y.  L^T entry at (j, i) is L[j*n + i]
  // for j > i (lower-triangle storage of L).
  for (int i = n - 1; i >= 0; --i) {
    double s = b[i];
    for (int j = i + 1; j < n; ++j) s -= L[j*n + i] * b[j];
    b[i] = s / L[i*n + i];
  }
}

/* ----------------------------------------------------------------------
   Banded Cholesky factor of symmetric positive-definite A stored in
   lower-band packed row-major form: AB[i*(bw+1) + (i-j)] holds A[i][j]
   for j in [max(0, i-bw), i].  AB[i*(bw+1) + 0] is the diagonal, AB[i*
   (bw+1) + d] is the d'th sub-diagonal entry of row i.

   Entries (i, j) with j > i are not stored (upper triangle implicit by
   symmetry); the band Cholesky factor inherits the same bandwidth as A
   (no fill outside the band).  Returns 0 on success, 1 if any pivot is
   non-positive.  Cost: ~n*(bw+1)^2/2 flops -- vs ~n^3/6 for dense.

   For bw = 0 the routine reduces to scaling the diagonal by 1/sqrt(.).
------------------------------------------------------------------------- */

int FixIlves::band_chol_factor(int n, int bw, double *AB)
{
  constexpr double TINY = 1.0e-30;
  const int row_stride = bw + 1;

  for (int k = 0; k < n; ++k) {
    double diag = AB[k*row_stride + 0];
    // diag -= sum_{j=max(0,k-bw)}^{k-1} L[k][j]^2
    const int jlo = (k - bw > 0) ? (k - bw) : 0;
    for (int j = jlo; j < k; ++j) {
      double Lkj = AB[k*row_stride + (k - j)];
      diag -= Lkj * Lkj;
    }
    if (diag < TINY) return 1;
    double Lkk = sqrt(diag);
    AB[k*row_stride + 0] = Lkk;
    double inv_Lkk = 1.0 / Lkk;

    // update L[i][k] for i in (k, min(n-1, k+bw)]
    int ihi = (k + bw < n - 1) ? (k + bw) : (n - 1);
    for (int i = k + 1; i <= ihi; ++i) {
      // d_ik = i - k in [1, bw], so L[i][k] is at AB[i*stride + (i-k)]
      double t = AB[i*row_stride + (i - k)];
      // subtract sum_j L[i][j] * L[k][j] where both rows have a non-zero
      // entry at column j.  L[i][j] is zero (and not stored) for j < i-bw,
      // so j must satisfy j >= i-bw.  L[k][j] is also zero for j < k-bw,
      // but i > k so i-bw > k-bw -- the binding constraint is j >= i-bw.
      const int jlo2 = (i - bw > 0) ? (i - bw) : 0;
      for (int j = jlo2; j < k; ++j) {
        t -= AB[i*row_stride + (i - j)] * AB[k*row_stride + (k - j)];
      }
      AB[i*row_stride + (i - k)] = t * inv_Lkk;
    }
  }
  return 0;
}

/* ----------------------------------------------------------------------
   Forward + back substitution to solve L L^T x = b given a lower-band
   packed L produced by band_chol_factor.  Operates in-place on b.
   Cost: ~n*(bw+1) flops per pass, so ~2*n*(bw+1) total.
------------------------------------------------------------------------- */

void FixIlves::band_chol_solve(int n, int bw, const double *AB, double *b)
{
  const int row_stride = bw + 1;

  // Forward substitution: L y = b
  for (int i = 0; i < n; ++i) {
    double s = b[i];
    const int jlo = (i - bw > 0) ? (i - bw) : 0;
    for (int j = jlo; j < i; ++j) {
      s -= AB[i*row_stride + (i - j)] * b[j];
    }
    b[i] = s / AB[i*row_stride + 0];
  }
  // Back substitution: L^T x = y.  L[j][i] for j > i is at row j, sub-
  // diagonal offset (j - i).
  for (int i = n - 1; i >= 0; --i) {
    double s = b[i];
    int jhi = (i + bw < n - 1) ? (i + bw) : (n - 1);
    for (int j = i + 1; j <= jhi; ++j) {
      s -= AB[j*row_stride + (j - i)] * b[j];
    }
    b[i] = s / AB[i*row_stride + 0];
  }
}

/* ----------------------------------------------------------------------
   Allocate / grow the per-cluster Cholesky factor cache.  Called from
   precompute_constraint_data() after cluster_offset[] is final.  Pool
   total size is sum_c n_c^2 doubles.  For variant != ILVES_FAST we do
   not allocate the pool at all (the solver doesn't reuse the factor).
------------------------------------------------------------------------- */

void FixIlves::grow_factor_cache()
{
  if (variant != ILVES_FAST) return;

  if (n_clusters + 1 > chol_offset_alloc) {
    chol_offset_alloc = n_clusters + 1;
    memory->grow(chol_pool_offset, chol_offset_alloc, "ilves:chol_pool_offset");
  }
  if (n_clusters > cluster_cached_alloc) {
    cluster_cached_alloc = n_clusters;
    memory->grow(cluster_cached, cluster_cached_alloc, "ilves:cluster_cached");
  }

  // CSR-style offsets in doubles.  Each cluster c stores its lower band-
  // packed L matrix as n_c rows of (cluster_bw[c] + 1) doubles each.
  chol_pool_offset[0] = 0;
  for (int c = 0; c < n_clusters; ++c) {
    int n_c = cluster_offset[c + 1] - cluster_offset[c];
    int bw_c = cluster_bw ? cluster_bw[c] : 0;
    chol_pool_offset[c + 1] = chol_pool_offset[c]
                            + (bigint) n_c * (bw_c + 1);
  }
  bigint need = chol_pool_offset[n_clusters];
  if (need > chol_pool_alloc) {
    bigint newcap = need + need / 4 + 16;
    memory->grow(chol_pool, newcap, "ilves:chol_pool");
    chol_pool_alloc = newcap;
  }
}

/* ----------------------------------------------------------------------
   Solve all constraints by ILVES Newton iteration.  For each cluster c:
     - assemble the local n_c x n_c Jacobian A_c (structurally symmetric,
       numerically asymmetric -- equals the exact Newton Jacobian)
     - assemble the rhs as the current constraint residuals g_k
     - LU-solve A_c * d-lambda = -g_c
     - update xshake using d-lambda
     - accumulate c_lambda[k] += d-lambda[k]
   Iterate until max relative bond-length violation falls below tolerance
   or max_iter is reached.

   Returns true if converged.  Phase 2 implementation is serial-only; for
   MPI the inter-cluster reduction will be added in Phase 4.
------------------------------------------------------------------------- */

bool FixIlves::solve_constraints()
{
  if (n_constr == 0) return true;

  grow_lu_workspace(largest_cluster);

  // zero accumulated Lagrange multipliers at start of step
  for (int k = 0; k < n_constr; ++k) c_lambda[k] = 0.0;

  // invalidate the per-cluster Cholesky factor cache: c_rx/c_ry/c_rz and
  // c_rsq were just recomputed in post_force, so every factor from the
  // previous step is stale.  Cache is repopulated lazily at iter == 0 of
  // each cluster (see below).  Only meaningful for ILVES_FAST.
  if (variant == ILVES_FAST && cluster_cached)
    for (int c = 0; c < n_clusters; ++c) cluster_cached[c] = 0;

  // tolerance is on relative bond-length violation: |s_k|/d_k - 1.
  // g_k = 0.5*(|s_k|^2 - d_k^2); for small violation,
  // g_k / d_k^2 ~ (|s_k|/d_k - 1).  Use that as the residual measure.
  const double tol_g = tolerance;

  bool converged = false;

  // 2-atom k = (a, b): gradient +r_k at a, -r_k at b; diagonal weight
  // (sum w_p^2/m_p) = 1/m_a + 1/m_b.
  auto diag_factor = [&](int k) -> double {
    return c_invma[k] + c_invmb[k];
  };

  // Off-diagonal A_kl contribution: sum over shared atoms p of
  //   w_k_p * w_l_p * (1/m_p) * dot
  // For two 2-atom constraints there is at most one shared atom, so this
  // is an early-return chain.  'dot' is the caller-supplied inner product
  // r_k.r_l (ILVES_FAST / velocity solve) or s_k.r_l (ILVES_FULL).
  auto offdiag = [&](int k, int l, double dot) -> double {
    tagint tag_ak = atom->tag[c_atom1[k]];
    tagint tag_bk = atom->tag[c_atom2[k]];
    tagint tag_al = atom->tag[c_atom1[l]];
    tagint tag_bl = atom->tag[c_atom2[l]];
    double invma_k = c_invma[k];
    double invmb_k = c_invmb[k];
    if (tag_ak == tag_al) return  invma_k * dot;
    if (tag_ak == tag_bl) return -invma_k * dot;
    if (tag_bk == tag_al) return -invmb_k * dot;
    if (tag_bk == tag_bl) return  invmb_k * dot;
    return 0.0;
  };

  for (int iter = 0; iter < max_iter; ++iter) {
    double max_relres = 0.0;

    for (int c = 0; c < n_clusters; ++c) {
      const int beg = cluster_offset[c];
      const int end = cluster_offset[c + 1];
      const int n_c = end - beg;

      // For ILVES_FAST the matrix entries depend only on c_invma/c_invmb/
      // c_rsq/c_r{x,y,z} -- all step-constant.  After we factor a cluster
      // once at iter == 0 we keep the L matrix in chol_pool[c] and only
      // re-do the cheap O(n_c^2) forward+back substitution on subsequent
      // iterations.  ILVES_FULL has an s_k-dependent Jacobian that must be
      // reassembled every Newton iter, so it doesn't benefit from caching.
      const bool can_reuse = (variant == ILVES_FAST) && cluster_cached
                          && cluster_cached[c];

      // s_k = (current iter's "bond" vector) under PBC minimum-image, rhs
      // g_k, and the residual indicator are always recomputed -- they
      // depend on xshake which changes between iterations.
      for (int s = 0; s < n_c; ++s) {
        int k = c_perm[beg + s];
        int a = c_atom1[k];
        int b = c_atom2[k];
        double sx = xshake[a][0] - xshake[b][0];
        double sy = xshake[a][1] - xshake[b][1];
        double sz = xshake[a][2] - xshake[b][2];
        domain->minimum_image(FLERR, sx, sy, sz);
        cl_sx[s] = sx;
        cl_sy[s] = sy;
        cl_sz[s] = sz;

        // residual (same for both variants): g_k = 0.5*(|s_k|^2 - d_k^2)
        double ssq = sx*sx + sy*sy + sz*sz;
        double gk = 0.5 * (ssq - c_dist[k]*c_dist[k]);
        lu_b[s] = -gk;

        double relres = fabs(gk) / (c_dist[k]*c_dist[k]);
        if (relres > max_relres) max_relres = relres;
      }

      int info = 0;

      if (can_reuse) {
        // Reuse the cached L (band format): only the cheap forward+back
        // substitution -- O(n_c*bw) work instead of O(n_c*bw^2) for a
        // re-factor or O(n_c^3) for dense.
        const int bw_c = cluster_bw[c];
        double *Lptr = chol_pool + chol_pool_offset[c];
        band_chol_solve(n_c, bw_c, Lptr, lu_b);
      } else if (variant == ILVES_FAST) {
        // ILVES_FAST: assemble symmetric A into the cluster's banded slot
        // in chol_pool, factor in-place (band Cholesky), and solve.  On
        // success the cached factor is reused on subsequent Newton iters.
        const int bw_c = cluster_bw[c];
        const int row_stride = bw_c + 1;
        double *Aptr = chol_pool + chol_pool_offset[c];

        // zero the cluster's band storage (n_c rows of (bw_c+1) doubles)
        for (int i = 0; i < n_c * row_stride; ++i) Aptr[i] = 0.0;

        // diagonal A_kk = (sum w_p^2/m_p) * |r_k|^2 -- diag_factor()
        // collapses to (1/m_a + 1/m_b) for 2-atom and
        // 1/m_B + 1/(4 m_A) + 1/(4 m_C) for 3-atom B-M.
        for (int s = 0; s < n_c; ++s) {
          int k = c_perm[beg + s];
          Aptr[s*row_stride + 0] = diag_factor(k) * c_rsq[k];
        }

        // off-diagonals: scan all pairs (s, t) with s < t in the cluster,
        // add A[t][s] = sum_shared_atom w_k_p * w_l_p * (r_k.r_l) / m_p.
        // After RCM the pairs that actually share an atom satisfy
        // |t-s| <= bw_c, so the band slot Aptr[t*stride + (t-s)] always
        // exists.  Pairs further apart silently contribute zero by
        // construction (no shared atom -> no contribution).
        for (int s = 0; s < n_c; ++s) {
          int k = c_perm[beg + s];
          int t_end = s + bw_c + 1;
          if (t_end > n_c) t_end = n_c;
          for (int t = s + 1; t < t_end; ++t) {
            int l = c_perm[beg + t];
            double rkrl = c_rx[k]*c_rx[l] + c_ry[k]*c_ry[l] + c_rz[k]*c_rz[l];
            double val = offdiag(k, l, rkrl);
            if (val != 0.0) Aptr[t*row_stride + (t - s)] = val;
          }
        }

        ++chol_calls;
        info = band_chol_factor(n_c, bw_c, Aptr);
        if (info == 0) {
          cluster_cached[c] = 1;
          band_chol_solve(n_c, bw_c, Aptr, lu_b);
        } else {
          ++chol_fallbacks;
          cluster_cached[c] = 0;
          // band-Cholesky failed (matrix not SPD).  Re-assemble densely
          // into lu_A and run LU.  rhs lu_b was not touched by the failed
          // factor so we don't need to rebuild it.
          for (int i = 0; i < n_c * n_c; ++i) lu_A[i] = 0.0;
          for (int s2 = 0; s2 < n_c; ++s2) {
            int k = c_perm[beg + s2];
            lu_A[s2*n_c + s2] = diag_factor(k) * c_rsq[k];
          }
          for (int s2 = 0; s2 < n_c; ++s2) {
            int k = c_perm[beg + s2];
            for (int t = s2 + 1; t < n_c; ++t) {
              int l = c_perm[beg + t];
              double rkrl = c_rx[k]*c_rx[l] + c_ry[k]*c_ry[l] + c_rz[k]*c_rz[l];
              double val = offdiag(k, l, rkrl);
              if (val != 0.0) {
                lu_A[s2*n_c + t] = val;
                lu_A[t*n_c + s2] = val;
              }
            }
          }
          info = lu_factor_solve(n_c);
        }
      } else {
        // ILVES_FULL: assemble the exact-Jacobian asymmetric matrix into
        // dense lu_A (cl_s* was already filled above with current s_k) and
        // solve with LU.  No caching: the matrix depends on s_k and so
        // changes every Newton iteration.  Diagonal uses s_k.r_k; off-
        // diagonal uses s_k.r_l (asymmetric: A_kl != A_lk in general).
        for (int i = 0; i < n_c * n_c; ++i) lu_A[i] = 0.0;
        for (int s = 0; s < n_c; ++s) {
          int k = c_perm[beg + s];
          double sr = cl_sx[s]*c_rx[k] + cl_sy[s]*c_ry[k] + cl_sz[s]*c_rz[k];
          lu_A[s*n_c + s] = diag_factor(k) * sr;
        }
        for (int s = 0; s < n_c; ++s) {
          int k = c_perm[beg + s];
          for (int t = 0; t < n_c; ++t) {
            if (t == s) continue;
            int l = c_perm[beg + t];
            double srl = cl_sx[s]*c_rx[l] + cl_sy[s]*c_ry[l] + cl_sz[s]*c_rz[l];
            double val = offdiag(k, l, srl);
            if (val != 0.0) lu_A[s*n_c + t] += val;
          }
        }
        info = lu_factor_solve(n_c);
      }

      if (info) {
        error->one(FLERR, Error::NOLASTLINE,
                   "Fix ilves: singular Jacobian in cluster {} (size {}, iter {}). "
                   "This usually indicates a degenerate or overdetermined "
                   "constraint topology.", c, n_c, iter);
      }

      // apply d-lambda to xshake and accumulate c_lambda.  Route partner
      // updates to the LOCAL owner via atom->map(tag).  For a true MPI
      // ghost (owner on another rank), atom->map returns the ghost index
      // -- skip the partner update locally; the partner rank holds the
      // same constraint and will update its local atom, and forward_comm
      // at the end of the iteration brings the corrected ghost position
      // back to this rank.
      for (int s = 0; s < n_c; ++s) {
        int k = c_perm[beg + s];
        double dl = lu_b[s];
        c_lambda[k] += dl;
        const int a = c_atom1[k];
        const int b = c_atom2[k];

        // Per-atom delta = w_p * dl * (1/m_p) * r_k: +r_k/m_a at a, -r_k/m_b at b.
        const double rx = c_rx[k];
        const double ry = c_ry[k];
        const double rz = c_rz[k];
        const double dl_a = dl * c_invma[k];
        const double dl_b = dl * c_invmb[k];

        // Route each atom's update to its local owner (atom->map(tag)
        // collapses PBC ghosts to local indices).  Updates to xshake[ghost]
        // would be overwritten by the end-of-iter forward_comm anyway, so
        // skip those.  PBC self-images route to their local owner via
        // atom->map so the forward_comm picks them up.
        int a_local = atom->map(atom->tag[a]);
        int b_local = atom->map(atom->tag[b]);
        if (a_local >= 0 && a_local < nlocal) {
          xshake[a_local][0] += dl_a * rx;
          xshake[a_local][1] += dl_a * ry;
          xshake[a_local][2] += dl_a * rz;
        }
        if (b_local >= 0 && b_local < nlocal) {
          xshake[b_local][0] -= dl_b * rx;
          xshake[b_local][1] -= dl_b * ry;
          xshake[b_local][2] -= dl_b * rz;
        }
      }
    }

    // refresh ghost xshake: keep redundantly-solving ranks in sync.
    comm->forward_comm(this);

    // synchronize convergence across MPI ranks: without this, ranks
    // with different residuals would exit the Newton loop at different
    // iterations and then deadlock on the next forward_comm call.
    double global_relres = max_relres;
    if (comm->nprocs > 1)
      MPI_Allreduce(&max_relres, &global_relres, 1, MPI_DOUBLE, MPI_MAX, world);
    if (global_relres < tol_g) { converged = true; break; }
  }

  if (!converged && (comm->me == 0)) {
    error->warning(FLERR,
                   "Fix ilves: max_iter ({}) reached without reaching tol={} at step {}",
                   max_iter, tolerance, update->ntimestep);
  }
  return converged;
}

/* ----------------------------------------------------------------------
   Convert the accumulated Lagrange multipliers into constraint forces
   added to atom->f, so that the next initial_integrate puts the atoms at
   the constrained xshake positions.

   For each constraint k = (a, b):
     f[a] += (lambda_k / dtfsq) * r_k
     f[b] -= (lambda_k / dtfsq) * r_k
   Forces are applied only to atoms that are local (< nlocal).  For ghost
   partners (MPI cross-rank bonds), the partner's owner rank applies its
   own copy of the force; Phase 4 will add reverse_comm for newton_bond=1
   cases where only one rank holds the constraint.
------------------------------------------------------------------------- */

void FixIlves::apply_constraint_forces(int vflag)
{
  if (n_constr == 0) return;
  const double inv_dtfsq = 1.0 / dtfsq;
  int *mask = atom->mask;
  double bond_v[6];

  if (store_flag) {
    if (maxstore < atom->nmax) {
      maxstore = atom->nmax;
      memory->destroy(fstore);
      memory->create(fstore, maxstore, 3, "ilves:fstore");
    }
    for (int i = 0; i < maxstore; ++i)
      fstore[i][0] = fstore[i][1] = fstore[i][2] = 0.0;
    array_atom = fstore;
  }

  for (int k = 0; k < n_constr; ++k) {
    double scale = c_lambda[k] * inv_dtfsq;
    double fx = scale * c_rx[k];
    double fy = scale * c_ry[k];
    double fz = scale * c_rz[k];
    int a = c_atom1[k];
    int b = c_atom2[k];

    // Route each endpoint's force to its local owner (atom->map(tag)
    // collapses PBC ghosts to local indices).  Any c_atom* may be a
    // ghost for cluster-completion constraints in MPI; that rank does
    // not apply the force locally -- the owning rank handles its own
    // atoms.
    int a_local = atom->map(atom->tag[a]);
    int b_local = atom->map(atom->tag[b]);

    if (a_local >= 0 && a_local < nlocal) {
      f[a_local][0] += fx; f[a_local][1] += fy; f[a_local][2] += fz;
      if (store_flag) {
        fstore[a_local][0] += fx; fstore[a_local][1] += fy; fstore[a_local][2] += fz;
      }
    }
    if (b_local >= 0 && b_local < nlocal) {
      f[b_local][0] -= fx; f[b_local][1] -= fy; f[b_local][2] -= fz;
      if (store_flag) {
        fstore[b_local][0] -= fx; fstore[b_local][1] -= fy; fstore[b_local][2] -= fz;
      }
    }

    if (evflag) {
      int atomlist[2];
      int count = 0;
      if (a_local >= 0 && a_local < nlocal) atomlist[count++] = a_local;
      if (b_local >= 0 && b_local < nlocal) atomlist[count++] = b_local;
      bond_v[0] = scale * c_rx[k]*c_rx[k];
      bond_v[1] = scale * c_ry[k]*c_ry[k];
      bond_v[2] = scale * c_rz[k]*c_rz[k];
      bond_v[3] = scale * c_rx[k]*c_ry[k];
      bond_v[4] = scale * c_rx[k]*c_rz[k];
      bond_v[5] = scale * c_ry[k]*c_rz[k];
      double fpair[1] = {scale};
      double dellist[1][3] = {{c_rx[k], c_ry[k], c_rz[k]}};
      int pairlist[1][2] = {{a_local, b_local}};
      v_tally(count, atomlist, 2.0, bond_v, nlocal, 1, pairlist, fpair, dellist);
    }
  }
  (void) mask;
  (void) vflag;
}

/* ----------------------------------------------------------------------
   Project atom positions onto the constraint surface at setup time.
   The trick (borrowed from fix shake's correct_coordinates): temporarily
   zero atom->f and atom->v, run the Newton solver to find constraint
   forces, then apply those forces directly to atom->x as position
   corrections (x -> x + dtfsq/m * f) and restore the saved f and v.

   This must be called with dtfsq = 0.5*dt^2*ftm2v.
------------------------------------------------------------------------- */

void FixIlves::correct_coordinates(int vflag)
{
  if (n_constr == 0) return;

  // save current f and v, then zero them so unconstrained_update reduces
  // to xshake = x (the current positions)
  double **f_save = nullptr;
  double **v_save = nullptr;
  memory->create(f_save, nlocal, 3, "ilves:fsave");
  memory->create(v_save, nlocal, 3, "ilves:vsave");
  for (int i = 0; i < nlocal; ++i) {
    f_save[i][0] = f[i][0]; f_save[i][1] = f[i][1]; f_save[i][2] = f[i][2];
    v_save[i][0] = v[i][0]; v_save[i][1] = v[i][1]; v_save[i][2] = v[i][2];
    f[i][0] = f[i][1] = f[i][2] = 0.0;
    v[i][0] = v[i][1] = v[i][2] = 0.0;
  }

  // run the Newton solver as in post_force.  xshake will be initialized
  // to x by unconstrained_update (since v = f = 0), and the solve finds
  // the constraint forces (added to atom->f) needed to project x to the
  // constraint surface via a single integration half-step.
  unconstrained_update();
  comm->forward_comm(this);
  solve_constraints();
  apply_constraint_forces(vflag);

  // project: x_new = x + dtfsq/m * f_c (where f is currently only f_c
  // because we zeroed f before solving)
  if (rmass) {
    for (int i = 0; i < nlocal; ++i) {
      double dfm = dtfsq / rmass[i];
      x[i][0] += dfm * f[i][0];
      x[i][1] += dfm * f[i][1];
      x[i][2] += dfm * f[i][2];
    }
  } else {
    for (int i = 0; i < nlocal; ++i) {
      double dfm = dtfsq / mass[type[i]];
      x[i][0] += dfm * f[i][0];
      x[i][1] += dfm * f[i][1];
      x[i][2] += dfm * f[i][2];
    }
  }

  // restore f and v
  for (int i = 0; i < nlocal; ++i) {
    f[i][0] = f_save[i][0]; f[i][1] = f_save[i][1]; f[i][2] = f_save[i][2];
    v[i][0] = v_save[i][0]; v[i][1] = v_save[i][1]; v[i][2] = v_save[i][2];
  }
  memory->destroy(f_save);
  memory->destroy(v_save);

  // updated x must be propagated to ghost atoms for the velocity correction
  // and subsequent force computation; xshake forward-comm is a convenient
  // way to do this since we already have the pack/unpack for xshake.
  // (positions stored in atom->x are also communicated by LAMMPS's normal
  //  flow, so this isn't strictly required here, but matches fix shake.)
  double **xtmp = xshake;
  xshake = x;
  comm->forward_comm(this);
  xshake = xtmp;

  // re-precompute reference vectors r_k = x[a]-x[b] now that x has been
  // updated -- subsequent solves (e.g. the post_force call following this
  // function) must use the corrected positions
  precompute_constraint_data();
}

/* ----------------------------------------------------------------------
   Remove velocity components along constrained directions so that
   d/dt(|x[a]-x[b]|^2) = 0 (i.e. (v[a]-v[b]) . r_k = 0) for every
   constraint k.  This is the velocity-RATTLE step.  At setup, the
   initial velocities from the data file (or velocity command) are
   generally not constraint-consistent; this routine projects them
   onto the constrained-velocity manifold while conserving linear
   momentum within each cluster.

   For each cluster:
     - assemble symmetric positive-definite matrix
         A_kl = sum_p sigma_k^p sigma_l^p (r_k . r_l) / m_p
       (the same shared-atom pattern as the position solve, but with
        s_k replaced by r_k -- the matrix is genuinely symmetric here)
     - rhs g_k = (v[a_k] - v[b_k]) . r_k
     - solve A * mu = g (note: rhs in lu_b will be +g, dlambda is -mu)
     - subtract mu_k from velocities:
         v[a] -= mu_k / m_a * r_k
         v[b] += mu_k / m_b * r_k

   Linear system -- a single solve is exact (no iteration needed).
------------------------------------------------------------------------- */

void FixIlves::correct_velocities()
{
  if (n_constr == 0) return;
  grow_lu_workspace(largest_cluster);

  // 2-atom helpers (mirror solve_constraints).
  auto diag_factor = [&](int k) -> double {
    return c_invma[k] + c_invmb[k];
  };
  auto offdiag = [&](int k, int l, double dot) -> double {
    tagint tag_ak = atom->tag[c_atom1[k]];
    tagint tag_bk = atom->tag[c_atom2[k]];
    tagint tag_al = atom->tag[c_atom1[l]];
    tagint tag_bl = atom->tag[c_atom2[l]];
    double invma_k = c_invma[k];
    double invmb_k = c_invmb[k];
    if (tag_ak == tag_al) return  invma_k * dot;
    if (tag_ak == tag_bl) return -invma_k * dot;
    if (tag_bk == tag_al) return -invmb_k * dot;
    if (tag_bk == tag_bl) return  invmb_k * dot;
    return 0.0;
  };

  // rhs for constraint k: g_k = (v[a] - v[b]) . r_k (time derivative of
  // the bond-length constraint).
  auto velocity_rhs = [&](int k) -> double {
    int a = c_atom1[k];
    int b = c_atom2[k];
    double vxd = v[a][0] - v[b][0];
    double vyd = v[a][1] - v[b][1];
    double vzd = v[a][2] - v[b][2];
    return vxd*c_rx[k] + vyd*c_ry[k] + vzd*c_rz[k];
  };

  for (int c = 0; c < n_clusters; ++c) {
    const int beg = cluster_offset[c];
    const int end = cluster_offset[c + 1];
    const int n_c = end - beg;

    // assemble symmetric A and rhs g
    for (int i = 0; i < n_c * n_c; ++i) lu_A[i] = 0.0;

    for (int s = 0; s < n_c; ++s) {
      int k = c_perm[beg + s];
      lu_A[s*n_c + s] = diag_factor(k) * c_rsq[k];
      lu_b[s] = velocity_rhs(k);
    }

    // off-diagonals: same shared-atom logic as the position solve, with
    // (r_k . r_l) -- symmetric in k,l so we only need s < t
    for (int s = 0; s < n_c; ++s) {
      int k = c_perm[beg + s];
      for (int t = s + 1; t < n_c; ++t) {
        int l = c_perm[beg + t];
        double rkrl = c_rx[k]*c_rx[l] + c_ry[k]*c_ry[l] + c_rz[k]*c_rz[l];
        double val = offdiag(k, l, rkrl);
        if (val != 0.0) {
          lu_A[s*n_c + t] = val;
          lu_A[t*n_c + s] = val;     // symmetric
        }
      }
    }

    // Velocity-projection matrix is always symmetric (r_k.r_l off-diagonals)
    // regardless of variant -- always SPD for non-degenerate clusters.
    ++chol_calls;
    int info = chol_factor_solve(n_c);
    if (info) {
      ++chol_fallbacks;
      // re-zero and re-assemble for the LU fallback (Cholesky overwrote L)
      for (int i = 0; i < n_c * n_c; ++i) lu_A[i] = 0.0;
      for (int s = 0; s < n_c; ++s) {
        int k = c_perm[beg + s];
        lu_A[s*n_c + s] = diag_factor(k) * c_rsq[k];
        lu_b[s] = velocity_rhs(k);
      }
      for (int s = 0; s < n_c; ++s) {
        int k = c_perm[beg + s];
        for (int t = s + 1; t < n_c; ++t) {
          int l = c_perm[beg + t];
          double rkrl = c_rx[k]*c_rx[l] + c_ry[k]*c_ry[l] + c_rz[k]*c_rz[l];
          double val = offdiag(k, l, rkrl);
          if (val != 0.0) {
            lu_A[s*n_c + t] = val;
            lu_A[t*n_c + s] = val;
          }
        }
      }
      info = lu_factor_solve(n_c);
    }
    if (info)
      error->one(FLERR,  Error::NOLASTLINE, "Fix ilves: singular velocity-correction matrix "
                 "in cluster {} (size {}).  Check for degenerate constraint topology.", c, n_c);

    // apply: v[p] += w_p * mu_k * (1/m_p) * (-r_k) for each participating
    // atom.  Note the OVERALL sign flips relative to the position solve:
    // the velocity rhs g_k uses +r_k for the v-dot-r convention, while
    // the constraint Jacobian is the same as for the position update.
    // Net effect for atom_p: dv = -w_p * mu_k / m_p * r_k.  For 2-atom
    // this reduces to v[a] -= mu * (1/m_a) * r_k and v[b] += mu *
    // (1/m_b) * r_k as in the original code.
    for (int s = 0; s < n_c; ++s) {
      int k = c_perm[beg + s];
      double mu = lu_b[s];
      int a = c_atom1[k];
      int b = c_atom2[k];
      int a_local = atom->map(atom->tag[a]);
      int b_local = atom->map(atom->tag[b]);
      double da = mu * c_invma[k];
      double db = mu * c_invmb[k];
      if (a_local >= 0 && a_local < nlocal) {
        v[a_local][0] -= da * c_rx[k];
        v[a_local][1] -= da * c_ry[k];
        v[a_local][2] -= da * c_rz[k];
      }
      if (b_local >= 0 && b_local < nlocal) {
        v[b_local][0] += db * c_rx[k];
        v[b_local][1] += db * c_ry[k];
        v[b_local][2] += db * c_rz[k];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   RATTLE-style velocity correction at the END_OF_STEP hook.  At this
   point the integrator has already advanced v to v(t+dt) using both
   half-steps -- we project v onto the constrained-velocity manifold so
   the dynamics stays on the constraint surface step-by-step instead of
   gradually redistributing KE/PE between bonded and non-bonded modes.

   c_rx/c_ry/c_rz/c_rsq are recomputed first because the bond vectors
   have rotated since post_force ran on this step.
------------------------------------------------------------------------- */

void FixIlves::end_of_step()
{
  x = atom->x;
  v = atom->v;
  f = atom->f;
  mass  = atom->mass;
  rmass = atom->rmass;
  type  = atom->type;
  nlocal = atom->nlocal;

  for (int k = 0; k < n_constr; ++k) {
    int a = c_atom1[k];
    int b = c_atom2[k];
    double rx = x[a][0] - x[b][0];
    double ry = x[a][1] - x[b][1];
    double rz = x[a][2] - x[b][2];
    domain->minimum_image(FLERR, rx, ry, rz);
    c_rx[k]  = rx;
    c_ry[k]  = ry;
    c_rz[k]  = rz;
    c_rsq[k] = rx*rx + ry*ry + rz*rz;
  }

  // forward-communicate atom->v so ghost velocities reflect the values
  // just-set on the owning rank's final_integrate.  Without this, ghost
  // v is stale (last value before exchange), and the cross-rank velocity
  // projection in correct_velocities computes inconsistent multipliers
  // on partner ranks, injecting energy step by step.
  if (comm->nprocs > 1) {
    comm_mode = 1;
    comm->forward_comm(this);
    comm_mode = 0;
  }

  correct_velocities();
}

void FixIlves::min_pre_reverse(int eflag, int /*vflag*/)
{
  eflag_pre_reverse = eflag;
}

/* ----------------------------------------------------------------------
   during minimization, the SHAKE/RATTLE machinery (which works through
   the integrator) is inert -- there are no velocities or time steps to
   step through.  We substitute the constraints with strong harmonic
   restraint potentials V_k = 0.5 * kbond * (|r_k| - d_k)^2, where kbond
   defaults to 1e9 * boltz (same as fix shake).  Energy and force are
   contributed to the minimizer just like a bond potential.
------------------------------------------------------------------------- */

void FixIlves::min_post_force(int vflag)
{
  // refresh atom pointers (allocation may have moved on reneighbor)
  x = atom->x;
  v = atom->v;
  f = atom->f;
  mass  = atom->mass;
  rmass = atom->rmass;
  type  = atom->type;
  nlocal = atom->nlocal;

  ev_init(eflag_pre_reverse, vflag);
  ebond = 0.0;

  // statistics accumulators
  const int nb = atom->nbondtypes  + 1;
  const int na = atom->nangletypes + 1;
  if (output_every) {
    for (int i = 0; i < nb; ++i) { b_count[i] = 0; b_ave[i] = b_max[i] = 0.0; b_min[i] = BIG; }
    for (int i = 0; i < na; ++i) { a_count[i] = 0; a_ave[i] = a_max[i] = 0.0; a_min[i] = BIG; }
  }

  // allocate / reset per-atom constraint-force storage
  if (store_flag) {
    if (maxstore < atom->nmax) {
      maxstore = atom->nmax;
      memory->destroy(fstore);
      memory->create(fstore, maxstore, 3, "ilves:fstore");
    }
    for (int i = 0; i < maxstore; ++i)
      fstore[i][0] = fstore[i][1] = fstore[i][2] = 0.0;
    array_atom = fstore;
  }

  tagint *tag = atom->tag;
  double bond_v[6];

  for (int k = 0; k < n_constr; ++k) {
    int a = c_atom1[k];
    int b = c_atom2[k];
    tagint ta = tag[a], tb = tag[b];

    // global dedup: only the rank locally owning the lowest-tag atom of
    // the constraint accounts for it (otherwise both ranks of a mirrored
    // cross-rank cluster would double-count this constraint's harmonic
    // restraint energy and the minimizer would see double its strength).
    tagint t_lower = (ta < tb) ? ta : tb;
    int lower_local = atom->map(t_lower);
    if (lower_local < 0 || lower_local >= nlocal) continue;

    int a_local = atom->map(ta);
    int b_local = atom->map(tb);

    double dx = x[a][0] - x[b][0];
    double dy = x[a][1] - x[b][1];
    double dz = x[a][2] - x[b][2];
    domain->minimum_image(FLERR, dx, dy, dz);
    double r = sqrt(dx*dx + dy*dy + dz*dz);
    double dr = r - c_dist[k];
    double rk = kbond * dr;
    double fbond = (r > 0.0) ? -2.0 * rk / r : 0.0;
    double eb = rk * dr;
    const double e_share = 0.5;

    int atomlist[2];
    int count = 0;
    if (a_local >= 0 && a_local < nlocal) {
      f[a_local][0] += dx * fbond;
      f[a_local][1] += dy * fbond;
      f[a_local][2] += dz * fbond;
      if (store_flag) {
        fstore[a_local][0] += dx * fbond;
        fstore[a_local][1] += dy * fbond;
        fstore[a_local][2] += dz * fbond;
      }
      atomlist[count++] = a_local;
      ebond += e_share * eb;
    }
    if (b_local >= 0 && b_local < nlocal) {
      f[b_local][0] -= dx * fbond;
      f[b_local][1] -= dy * fbond;
      f[b_local][2] -= dz * fbond;
      if (store_flag) {
        fstore[b_local][0] -= dx * fbond;
        fstore[b_local][1] -= dy * fbond;
        fstore[b_local][2] -= dz * fbond;
      }
      atomlist[count++] = b_local;
      ebond += e_share * eb;
    }
    if (evflag) {
      bond_v[0] = 0.5 * dx*dx * fbond;
      bond_v[1] = 0.5 * dy*dy * fbond;
      bond_v[2] = 0.5 * dz*dz * fbond;
      bond_v[3] = 0.5 * dx*dy * fbond;
      bond_v[4] = 0.5 * dx*dz * fbond;
      bond_v[5] = 0.5 * dy*dz * fbond;
      ev_tally(count, atomlist, 2.0, eb, bond_v);
    }

    if (output_every) {
      int t = c_type[k];
      if (t > 0) {
        b_count[t]++;
        b_ave[t] += r;
        if (r > b_max[t]) b_max[t] = r;
        if (r < b_min[t]) b_min[t] = r;
      } else {
        int at = -t;
        a_count[at]++;
        a_ave[at] += r;
        if (r > a_max[at]) a_max[at] = r;
        if (r < a_min[at]) a_min[at] = r;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   per-atom array lifecycle
------------------------------------------------------------------------- */

double FixIlves::memory_usage()
{
  double bytes = 0.0;
  bytes += (double) atom->nmax * sizeof(int);          // ilves_flag
  bytes += (double) atom->nmax * 3 * sizeof(double);   // xshake
  if (store_flag) bytes += (double) maxstore * 3 * sizeof(double);
  bytes += (double) max_constr * (5*sizeof(int) + 3*sizeof(double));
  // Replicated global topology: see gather_global_topology() for the
  // component breakdown.  These bytes are present on every MPI rank.
  bytes += 2.0 * (double) gb_a.size()    * sizeof(tagint);
  bytes +=       (double) gb_type.size() * sizeof(int);
  bytes += 3.0 * (double) ga1.size()     * sizeof(tagint);
  bytes +=       (double) ga_type.size() * sizeof(int);
  bytes += 48.0 * (double) tag_cluster.size();
  // Cholesky factor cache (only allocated for the "fast" variant).
  bytes += (double) chol_pool_alloc * sizeof(double);
  return bytes;
}

void FixIlves::grow_arrays(int nmax)
{
  memory->grow(ilves_flag, nmax, "ilves:ilves_flag");
  memory->destroy(xshake);
  memory->create(xshake, nmax, 3, "ilves:xshake");
}

void FixIlves::copy_arrays(int i, int j, int /*delflag*/)
{
  ilves_flag[j] = ilves_flag[i];
  // xshake is a transient working buffer rebuilt each step; no need to copy.
}

void FixIlves::set_arrays(int i)
{
  ilves_flag[i] = 0;
}

/* ---------------------------------------------------------------------- */

int FixIlves::pack_exchange(int i, double *buf)
{
  buf[0] = (double) ilves_flag[i];
  return 1;
}

int FixIlves::unpack_exchange(int nlocal_in, double *buf)
{
  ilves_flag[nlocal_in] = (int) buf[0];
  return 1;
}

/* ---------------------------------------------------------------------- */

int FixIlves::pack_forward_comm(int n, int *list, double *buf,
                                int pbc_flag, int *pbc)
{
  int m = 0;
  // comm_mode == 1 means pack atom->v (no PBC shift needed for velocities)
  if (comm_mode == 1) {
    for (int i = 0; i < n; ++i) {
      int j = list[i];
      buf[m++] = atom->v[j][0];
      buf[m++] = atom->v[j][1];
      buf[m++] = atom->v[j][2];
    }
    return m;
  }
  // default (comm_mode == 0): pack xshake with PBC shift applied
  if (pbc_flag == 0) {
    for (int i = 0; i < n; ++i) {
      int j = list[i];
      buf[m++] = xshake[j][0];
      buf[m++] = xshake[j][1];
      buf[m++] = xshake[j][2];
    }
  } else {
    double dx, dy, dz;
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    for (int i = 0; i < n; ++i) {
      int j = list[i];
      buf[m++] = xshake[j][0] + dx;
      buf[m++] = xshake[j][1] + dy;
      buf[m++] = xshake[j][2] + dz;
    }
  }
  return m;
}

void FixIlves::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  if (comm_mode == 1) {
    for (int i = first; i < last; ++i) {
      atom->v[i][0] = buf[m++];
      atom->v[i][1] = buf[m++];
      atom->v[i][2] = buf[m++];
    }
    return;
  }
  for (int i = first; i < last; ++i) {
    xshake[i][0] = buf[m++];
    xshake[i][1] = buf[m++];
    xshake[i][2] = buf[m++];
  }
}

/* ----------------------------------------------------------------------
   count DOFs removed: every constraint subtracts one (across all ranks)
   the global sum is what the integrator uses for temperature normalization
------------------------------------------------------------------------- */

bigint FixIlves::dof(int igroup)
{
  // each constraint k contributes 1 dof iff both atoms are in igroup.
  // With cross-rank cluster mirroring (newton on or off) each constraint
  // can appear on multiple ranks.  Dedupe by counting a constraint only
  // on the rank that locally owns its LOWER-TAG atom.  This way every
  // constraint contributes exactly one dof globally.
  int gbit = group->bitmask[igroup];
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  const int nlocal_now = atom->nlocal;
  bigint n = 0;
  for (int k = 0; k < n_constr; ++k) {
    int a = c_atom1[k];
    int b = c_atom2[k];
    tagint ta = tag[a], tb = tag[b];
    tagint t_lower = (ta < tb) ? ta : tb;
    int lower_local = atom->map(t_lower);
    if (lower_local < 0 || lower_local >= nlocal_now) continue;     // partner rank counts it

    // group filter: need both atoms in group.  Use whichever local
    // copies are available (ghost mask is valid too via comm shells).
    int a_local = atom->map(ta);
    int b_local = atom->map(tb);
    if (a_local < 0 || b_local < 0) continue;
    if (!(mask[a_local] & gbit)) continue;
    if (!(mask[b_local] & gbit)) continue;
    ++n;
  }
  bigint nall;
  MPI_Allreduce(&n, &nall, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  return nall;
}

/* ---------------------------------------------------------------------- */

void FixIlves::reset_dt()
{
  if (utils::strmatch(update->integrate_style, "^verlet")) {
    dtv   = update->dt;
    dtfsq = update->dt * update->dt * force->ftm2v;
    respa = 0;
  } else {
    auto *respa_ptr = dynamic_cast<Respa *>(update->integrate);
    if (!respa_ptr)
      error->all(FLERR, Error::NOLASTLINE,
                 "Failure to access Respa style {}", update->integrate_style);
    respa = 1;
    nlevels_respa = respa_ptr->nlevels;
    loop_respa = respa_ptr->loop;
    step_respa = respa_ptr->step;
    dtv = step_respa[0];
    dtf_inner = step_respa[0] * force->ftm2v;
    // dtfsq is recomputed each level inside unconstrained_update_respa
  }
}

/* ----------------------------------------------------------------------
   compute_scalar: total restraint energy (only nonzero during minimization;
   Phase 5 will fill this with sum_k 0.5*kbond*(|r_ab|-d_k)^2)
------------------------------------------------------------------------- */

double FixIlves::compute_scalar()
{
  double e_all;
  MPI_Allreduce(&ebond, &e_all, 1, MPI_DOUBLE, MPI_SUM, world);
  return e_all;
}

/* ----------------------------------------------------------------------
   bond/angle statistics (current bond length distribution per type)
------------------------------------------------------------------------- */

void FixIlves::stats()
{
  if (!output_every) return;

  const int nb = atom->nbondtypes  + 1;
  const int na = atom->nangletypes + 1;

  for (int i = 0; i < nb; ++i) { b_count[i] = 0; b_ave[i] = b_max[i] = 0.0; b_min[i] = BIG; }
  for (int i = 0; i < na; ++i) { a_count[i] = 0; a_ave[i] = a_max[i] = 0.0; a_min[i] = BIG; }

  const int nlocal_now = atom->nlocal;
  tagint *tag = atom->tag;
  for (int k = 0; k < n_constr; ++k) {
    int a = c_atom1[k];
    int b = c_atom2[k];
    tagint ta = tag[a], tb = tag[b];

    // dedup: only count this constraint on the rank that owns the lower-tag
    // endpoint locally.  Required because the same constraint can live on
    // multiple ranks when its cluster spans subdomains (we solve them
    // redundantly).
    {
      tagint t_lower = (ta < tb) ? ta : tb;
      int lower_local = atom->map(t_lower);
      if (lower_local < 0 || lower_local >= nlocal_now) continue;
    }

    double dx = atom->x[a][0] - atom->x[b][0];
    double dy = atom->x[a][1] - atom->x[b][1];
    double dz = atom->x[a][2] - atom->x[b][2];
    domain->minimum_image(FLERR, dx, dy, dz);
    double r = sqrt(dx*dx + dy*dy + dz*dz);
    int t = c_type[k];
    if (t > 0) {
      b_count[t]++;
      b_ave[t] += r;
      if (r > b_max[t]) b_max[t] = r;
      if (r < b_min[t]) b_min[t] = r;
    } else {
      int at = -t;
      a_count[at]++;
      a_ave[at] += r;
      if (r > a_max[at]) a_max[at] = r;
      if (r < a_min[at]) a_min[at] = r;
    }
  }

  MPI_Allreduce(b_count, b_count_all, nb, MPI_LMP_BIGINT, MPI_SUM, world);
  MPI_Allreduce(b_ave,   b_ave_all,   nb, MPI_DOUBLE,     MPI_SUM, world);
  MPI_Allreduce(b_max,   b_max_all,   nb, MPI_DOUBLE,     MPI_MAX, world);
  MPI_Allreduce(b_min,   b_min_all,   nb, MPI_DOUBLE,     MPI_MIN, world);
  MPI_Allreduce(a_count, a_count_all, na, MPI_LMP_BIGINT, MPI_SUM, world);
  MPI_Allreduce(a_ave,   a_ave_all,   na, MPI_DOUBLE,     MPI_SUM, world);
  MPI_Allreduce(a_max,   a_max_all,   na, MPI_DOUBLE,     MPI_MAX, world);
  MPI_Allreduce(a_min,   a_min_all,   na, MPI_DOUBLE,     MPI_MIN, world);

  if (comm->me == 0) {
    auto mesg = fmt::format("ILVES stats on step {} (type / ave / delta / count)\n"
                            "  bonds: length (Angstrom); angles: value (degrees)\n",
                            update->ntimestep);
    for (int i = 1; i < nb; ++i) {
      if (b_count_all[i])
        mesg += fmt::format("Bond:  {:>4d}   {:<9.6} {:<11.6} {:>8d}\n",
                            i, b_ave_all[i] / b_count_all[i],
                            b_max_all[i] - b_min_all[i], b_count_all[i]);
    }
    for (int i = 1; i < na; ++i) {
      if (a_count_all[i]) {
        double r1 = angle_r1[i], r2 = angle_r2[i];
        if (r1 > 0.0 && r2 > 0.0) {
          // convert end-to-end distances to angles (degrees) via law of cosines
          double inv2r1r2 = 1.0 / (2.0 * r1 * r2);
          double r1sq_r2sq = r1 * r1 + r2 * r2;
          auto dist_to_deg = [&](double r) {
            double c = (r1sq_r2sq - r * r) * inv2r1r2;
            if (c > 1.0) c = 1.0;
            if (c < -1.0) c = -1.0;
            return acos(c) * (180.0 / MY_PI);
          };
          double theta_ave = dist_to_deg(a_ave_all[i] / a_count_all[i]);
          double theta_min = dist_to_deg(a_min_all[i]);
          double theta_max = dist_to_deg(a_max_all[i]);
          mesg += fmt::format("Angle: {:>4d}   {:<9.4f} {:<11.6} {:>8d}\n",
                              i, theta_ave, theta_max - theta_min, a_count_all[i]);
        } else {
          mesg += fmt::format("Angle: {:>4d}   {:<9.6} {:<11.6} {:>8d}\n",
                              i, a_ave_all[i] / a_count_all[i],
                              a_max_all[i] - a_min_all[i], a_count_all[i]);
        }
      }
    }
    if (chol_calls > 0)
      mesg += fmt::format("Cholesky: {} calls, {} LU fallbacks ({:.4}%)\n",
                          chol_calls, chol_fallbacks,
                          100.0 * (double) chol_fallbacks / (double) chol_calls);
    utils::logmesg(lmp, mesg);
  }

  bigint nt = update->ntimestep;
  next_output = nt + output_every;
  if (nt % output_every != 0)
    next_output = (nt / output_every) * output_every + output_every;
}

/* ----------------------------------------------------------------------
   init_topology: gather the global bond/angle table, build the
   tag-to-cluster map, and fill angle_distance[] for selected angle types
   using the replicated angle table.  Runs once from FixIlves::init().
------------------------------------------------------------------------- */

void FixIlves::init_topology()
{
  gather_global_topology();

  // angle distances: for each constrained angle type, scan the global
  // angle table for the bond types of its flanking bonds.  Computing this
  // from the global table avoids problems where a given rank's local
  // angle is of a different angle type or where the relevant bond is
  // stored on a remote rank with newton_bond=on.
  //
  // Filtering: an angle is constrained only when all three atoms are in
  // the fix group AND both flanking bonds are themselves constrained
  // (see build_constraint_list).  Angles that fail either test are
  // irrelevant to the angle_distance computation and must not feed into
  // the mixed-bond check below -- otherwise an angle type shared by both
  // constrained and unconstrained instances (e.g. one with type-5 bonds
  // selected by `b 5` and another with type-1 bonds outside the group)
  // would trigger a spurious "mixed type bonds" error.
  if (has_angle) {
    const int na = (int) ga1.size();

    // Per-angle filter: 1 if this rank can verify the angle is constrained,
    // 0 otherwise.  At init() time ghosts are not yet built, so atoms on
    // remote ranks are not locally accessible (atom->map returns -1) and
    // mask/type lookups are unavailable for them; an MPI_Allreduce(MAX)
    // below combines per-rank determinations, so an angle is kept as long
    // as some rank can confirm it.
    std::vector<char> ang_constrain(na, 0);
    int *mask_now = atom->mask;
    for (int gi = 0; gi < na; ++gi) {
      int at = ga_type[gi];
      if (at <= 0 || at > atom->nangletypes) continue;
      if (!angle_flag[at]) continue;
      tagint t1 = ga1[gi];
      tagint t2 = ga2[gi];
      tagint t3 = ga3[gi];
      int i1 = atom->map(t1);
      int i2 = atom->map(t2);
      int i3 = atom->map(t3);
      if (i1 < 0 || i2 < 0 || i3 < 0) continue;
      if (!(mask_now[i1] & groupbit)) continue;
      if (!(mask_now[i2] & groupbit)) continue;
      if (!(mask_now[i3] & groupbit)) continue;
      if (!bond_is_constrained(t2, t1)) continue;
      if (!bond_is_constrained(t2, t3)) continue;
      ang_constrain[gi] = 1;
    }
    if ((comm->nprocs > 1) && (na > 0))
      MPI_Allreduce(MPI_IN_PLACE, ang_constrain.data(), na,
                    MPI_CHAR, MPI_MAX, world);

    auto lookup_bt = [&](tagint t1, tagint t2) -> int {
      tagint lo_t = MIN(t1, t2), hi_t = MAX(t1, t2);
      int lo = 0, hi = (int) gb_a.size();
      while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (gb_a[mid] < lo_t || (gb_a[mid] == lo_t && gb_b[mid] < hi_t)) lo = mid + 1;
        else hi = mid;
      }
      if (lo == (int) gb_a.size() || gb_a[lo] != lo_t || gb_b[lo] != hi_t) return 0;
      return gb_type[lo];
    };

    for (int at = 1; at <= atom->nangletypes; ++at) {
      if (!angle_flag[at]) { angle_distance[at] = 0.0; continue; }

      int b1 = 0, b2 = 0;
      int conflict = 0;
      for (int gi = 0; gi < na; ++gi) {
        if (ga_type[gi] != at) continue;
        if (!ang_constrain[gi]) continue;
        tagint tA = ga1[gi];
        tagint tB = ga2[gi];     // middle atom
        tagint tC = ga3[gi];
        int bt1 = lookup_bt(tA, tB);
        int bt2 = lookup_bt(tB, tC);
        if (bt1 == 0 || bt2 == 0) continue;
        int bmin = MIN(bt1, bt2);
        int bmax = MAX(bt1, bt2);
        if (b1 == 0) { b1 = bmin; b2 = bmax; }
        else if (bmin != b1 || bmax != b2) { conflict = 1; break; }
      }
      if (conflict)
        error->all(FLERR, Error::NOLASTLINE, "Fix ilves: angle type {} has mixed type bonds", at);
      if (b1 == 0) { angle_distance[at] = 0.0; continue; }
      const double theta0 = force->angle->equilibrium_angle(at);
      const double r1 = bond_distance[b1];
      const double r2 = bond_distance[b2];
      angle_r1[at] = r1;
      angle_r2[at] = r2;
      angle_distance[at] = sqrt(r1*r1 + r2*r2 - 2.0*r1*r2*cos(theta0));
    }
  }

  record_topology_baseline();
}

/* ----------------------------------------------------------------------
   Walk this rank's locally-stored bond and angle slots once, count those
   of constrained types, and MPI_Allreduce to globals.  The bond/angle
   counts are RAW per-slot sums -- not deduplicated for newton_bond -- so
   under newton off they are 2x the unique-bond count.  That's fine: we
   only ever compare these numbers against the at-init baseline, and the
   newton mode doesn't change within a run, so the absolute scale is
   irrelevant.  natoms is reported as-is (already a global value).
------------------------------------------------------------------------- */

void FixIlves::count_constrained_topology(bigint &natoms,
                                          bigint &nconstrbonds,
                                          bigint &nconstrangles) const
{
  int *num_bond_l    = atom->num_bond;
  int **bond_type_l  = atom->bond_type;
  int *num_angle_l   = atom->num_angle;
  int **angle_type_l = atom->angle_type;
  const int nlocal_now = atom->nlocal;
  const int nbtypes = atom->nbondtypes;
  const int natypes = atom->nangletypes;

  bigint local_bonds = 0;
  bigint local_angles = 0;

  if (num_bond_l && bond_type_l) {
    for (int i = 0; i < nlocal_now; ++i) {
      const int nb = num_bond_l[i];
      for (int b = 0; b < nb; ++b) {
        int bt = bond_type_l[i][b];
        if (bt == 0) continue;
        if (bt < 0) bt = -bt;
        if (bt <= nbtypes && bond_flag && bond_flag[bt]) ++local_bonds;
      }
    }
  }

  if (num_angle_l && angle_type_l) {
    for (int i = 0; i < nlocal_now; ++i) {
      const int na = num_angle_l[i];
      for (int m = 0; m < na; ++m) {
        int at = angle_type_l[i][m];
        if (at == 0) continue;
        if (at < 0) at = -at;
        if (at <= natypes && angle_flag && angle_flag[at]) ++local_angles;
      }
    }
  }

  bigint locals[2] = {local_bonds, local_angles};
  bigint globals[2];
  MPI_Allreduce(locals, globals, 2, MPI_LMP_BIGINT, MPI_SUM, world);

  natoms        = atom->natoms;
  nconstrbonds  = globals[0];
  nconstrangles = globals[1];
}

/* ---------------------------------------------------------------------- */

void FixIlves::record_topology_baseline()
{
  count_constrained_topology(natoms_at_init,
                             nconstrbonds_sum_at_init,
                             nconstrangles_sum_at_init);
}

/* ----------------------------------------------------------------------
   Compare the current constrained-topology counts against the baseline
   recorded at init_topology() time.  Any mismatch indicates a sibling
   fix (e.g. fix bond/create, fix bond/break, delete_atoms, create_atoms)
   has modified atoms or constrained bond/angle slots during the run --
   fix ilves does not support dynamic topology and would silently use
   stale gathered tables, so we abort with a clear error.
------------------------------------------------------------------------- */

void FixIlves::check_topology_unchanged()
{
  bigint natoms_now, nbonds_now, nangles_now;
  count_constrained_topology(natoms_now, nbonds_now, nangles_now);

  if (natoms_now != natoms_at_init) {
    error->all(FLERR, Error::NOLASTLINE,
               "Fix ilves: total atom count changed from {} (at init) to {} "
               "during the run.  Fix ilves does not support dynamic topology; "
               "rerun with a fixed atom set or invoke a fresh run command "
               "after the topology change so init() can re-gather.",
               natoms_at_init, natoms_now);
  }
  if (nbonds_now != nconstrbonds_sum_at_init) {
    error->all(FLERR, Error::NOLASTLINE,
               "Fix ilves: number of constrained bond slots changed from {} "
               "(at init) to {} during the run.  A sibling fix (fix bond/create, "
               "fix bond/break, set type, ...) appears to have modified the "
               "constrained topology.  Fix ilves does not support dynamic "
               "topology; rerun with a fixed topology.",
               nconstrbonds_sum_at_init, nbonds_now);
  }
  if (nangles_now != nconstrangles_sum_at_init) {
    error->all(FLERR, Error::NOLASTLINE,
               "Fix ilves: number of constrained angle slots changed from {} "
               "(at init) to {} during the run.  A sibling fix has modified "
               "the constrained topology.  Fix ilves does not support dynamic "
               "topology; rerun with a fixed topology.",
               nconstrangles_sum_at_init, nangles_now);
  }
}

/* ----------------------------------------------------------------------
   Gather global bond and angle topology onto every rank via MPI_Allgatherv.
   We store (lower_tag, higher_tag, type) for bonds and
   (atom1_tag, middle_tag, atom3_tag, type) for angles -- deduplicated.

   This is called once at init and the result lives in the gb_ and ga_
   member arrays.  The topology of bonds and angles is treated as fixed
   for the lifetime of the fix; commands that change topology mid-run
   (e.g. fix bond/create) are not currently supported by fix ilves.

   The global topology enables a uniform constraint-list build that does
   not depend on which rank stores a given bond (newton_bond on vs off).
   Every rank that owns at least one atom of a given bond or angle will
   add the corresponding constraint to its local list and apply the
   constraint force to its local atoms only.  Both ranks of a cross-rank
   constraint thus compute the same Lagrange multiplier from synchronized
   positions and apply equal-and-opposite forces to the two atoms.
------------------------------------------------------------------------- */

void FixIlves::gather_global_topology()
{
  tagint *tag       = atom->tag;
  int **nb_type     = atom->bond_type;
  tagint **nb_atom  = atom->bond_atom;
  int *num_bond     = atom->num_bond;
  int **na_type     = atom->angle_type;
  tagint **na_atom1 = atom->angle_atom1;
  tagint **na_atom2 = atom->angle_atom2;
  tagint **na_atom3 = atom->angle_atom3;
  int *num_angle    = atom->num_angle;
  const int nlocal_now = atom->nlocal;

  // -----------------------------------------------------------------
  // bonds: pack (min_tag, max_tag, |type|) for each locally-stored bond
  // -----------------------------------------------------------------
  std::vector<tagint> sendb;
  sendb.reserve(3 * nlocal_now);
  for (int i = 0; i < nlocal_now; ++i) {
    tagint ti = tag[i];
    for (int b = 0; b < num_bond[i]; ++b) {
      int bt = nb_type[i][b];
      if (bt == 0) continue;
      if (bt < 0) bt = -bt;
      tagint tj = nb_atom[i][b];
      sendb.push_back(MIN(ti, tj));
      sendb.push_back(MAX(ti, tj));
      sendb.push_back((tagint) bt);
    }
  }
  int my_b = (int) sendb.size();
  std::vector<int> cb(comm->nprocs), db(comm->nprocs);
  MPI_Allgather(&my_b, 1, MPI_INT, cb.data(), 1, MPI_INT, world);
  int totb = 0;
  for (int p = 0; p < comm->nprocs; ++p) { db[p] = totb; totb += cb[p]; }
  std::vector<tagint> recvb(totb);
  MPI_Allgatherv(sendb.data(), my_b, MPI_LMP_TAGINT,
                 recvb.data(), cb.data(), db.data(), MPI_LMP_TAGINT, world);

  // dedup bonds: sort by (min_tag, max_tag), keep unique
  struct BondEntry { tagint a, b; int type; };
  std::vector<BondEntry> bonds;
  bonds.reserve(totb / 3);
  for (int i = 0; i < totb; i += 3) {
    bonds.push_back({recvb[i], recvb[i+1], (int) recvb[i+2]});
  }
  std::sort(bonds.begin(), bonds.end(),
            [](const BondEntry &x, const BondEntry &y) {
              if (x.a != y.a) return x.a < y.a;
              if (x.b != y.b) return x.b < y.b;
              return x.type < y.type;
            });
  auto last_b = std::unique(bonds.begin(), bonds.end(),
                            [](const BondEntry &x, const BondEntry &y) {
                              return x.a == y.a && x.b == y.b;
                            });
  bonds.erase(last_b, bonds.end());

  gb_a.clear(); gb_b.clear(); gb_type.clear();
  gb_a.reserve(bonds.size()); gb_b.reserve(bonds.size()); gb_type.reserve(bonds.size());
  for (auto &e : bonds) { gb_a.push_back(e.a); gb_b.push_back(e.b); gb_type.push_back(e.type); }

  // -----------------------------------------------------------------
  // angles: pack (atom1_tag, atom2_tag, atom3_tag, |type|) for each
  // locally-stored angle.  atom2 = middle atom.  Pre-filter by angle_flag
  // to keep ga_ down to only the angle types we will constrain.
  // -----------------------------------------------------------------
  std::vector<tagint> senda;
  for (int i = 0; i < nlocal_now; ++i) {
    for (int m = 0; m < num_angle[i]; ++m) {
      int at = na_type[i][m];
      if (at == 0) continue;
      if (at < 0) at = -at;
      if (at > atom->nangletypes || !angle_flag[at]) continue;
      // dedupe: only middle atom owner packs the angle entry
      if (na_atom2[i][m] != tag[i]) continue;
      tagint o1 = na_atom1[i][m];
      tagint o3 = na_atom3[i][m];
      tagint t1 = MIN(o1, o3);
      tagint t3 = MAX(o1, o3);
      senda.push_back(t1);
      senda.push_back(na_atom2[i][m]);
      senda.push_back(t3);
      senda.push_back((tagint) at);
    }
  }
  int my_a = (int) senda.size();
  std::vector<int> ca(comm->nprocs), da(comm->nprocs);
  MPI_Allgather(&my_a, 1, MPI_INT, ca.data(), 1, MPI_INT, world);
  int tota = 0;
  for (int p = 0; p < comm->nprocs; ++p) { da[p] = tota; tota += ca[p]; }
  std::vector<tagint> recva(tota);
  MPI_Allgatherv(senda.data(), my_a, MPI_LMP_TAGINT,
                 recva.data(), ca.data(), da.data(), MPI_LMP_TAGINT, world);

  struct AngleEntry { tagint a, b, c; int type; };
  std::vector<AngleEntry> angles;
  angles.reserve(tota / 4);
  for (int i = 0; i < tota; i += 4) {
    angles.push_back({recva[i], recva[i+1], recva[i+2], (int) recva[i+3]});
  }
  std::sort(angles.begin(), angles.end(),
            [](const AngleEntry &x, const AngleEntry &y) {
              if (x.b != y.b) return x.b < y.b;
              if (x.a != y.a) return x.a < y.a;
              if (x.c != y.c) return x.c < y.c;
              return x.type < y.type;
            });
  auto last_a = std::unique(angles.begin(), angles.end(),
                            [](const AngleEntry &x, const AngleEntry &y) {
                              return x.a == y.a && x.b == y.b && x.c == y.c;
                            });
  angles.erase(last_a, angles.end());

  ga1.clear(); ga2.clear(); ga3.clear(); ga_type.clear();
  ga1.reserve(angles.size()); ga2.reserve(angles.size());
  ga3.reserve(angles.size()); ga_type.reserve(angles.size());
  for (auto &e : angles) {
    ga1.push_back(e.a); ga2.push_back(e.b); ga3.push_back(e.c);
    ga_type.push_back(e.type);
  }

  // -----------------------------------------------------------------
  // Build a global cluster-id table via union-find over all involved
  // tags (linked: bond a-b, and angle legs a-b, b-c).
  // -----------------------------------------------------------------
  tag_cluster.clear();
  tag_cluster.reserve(gb_a.size() * 2 + ga1.size() * 3);
  for (size_t i = 0; i < gb_a.size(); ++i) {
    tag_cluster.try_emplace(gb_a[i], gb_a[i]);
    tag_cluster.try_emplace(gb_b[i], gb_b[i]);
  }
  for (size_t i = 0; i < ga1.size(); ++i) {
    tag_cluster.try_emplace(ga1[i], ga1[i]);
    tag_cluster.try_emplace(ga2[i], ga2[i]);
    tag_cluster.try_emplace(ga3[i], ga3[i]);
  }

  auto find = [&](tagint t) {
    auto it = tag_cluster.find(t);
    while (it->second != t) {
      auto pit = tag_cluster.find(it->second);
      it->second = pit->second;
      t = it->second;
      it = pit;
    }
    return t;
  };
  auto unite = [&](tagint a, tagint b) {
    tagint ra = find(a), rb = find(b);
    if (ra != rb) tag_cluster[ra] = rb;
  };

  for (size_t i = 0; i < gb_a.size(); ++i) unite(gb_a[i], gb_b[i]);
  for (size_t i = 0; i < ga1.size(); ++i) {
    unite(ga1[i], ga2[i]);
    unite(ga2[i], ga3[i]);
  }
  for (auto &kv : tag_cluster) kv.second = find(kv.first);

  if (comm->me == 0) {
    // Estimate the bytes held by the replicated global topology tables.
    // These are duplicated on every MPI rank, so the total fix-ilves
    // memory at high rank counts can be (Nranks * this) -- the
    // dominant cost for very large systems.  Components:
    //   gb_a, gb_b   : 2 * nbonds * sizeof(tagint)
    //   gb_type      : nbonds * sizeof(int)
    //   ga1, ga2, ga3: 3 * nangles * sizeof(tagint)
    //   ga_type      : nangles * sizeof(int)
    //   tag_cluster  : ~48 bytes per involved tag (libstdc++ unordered_map
    //                  node + bucket; an over-estimate for libc++)
    const bigint sz_tagint = (bigint) sizeof(tagint);
    const bigint sz_int    = (bigint) sizeof(int);
    bigint bytes = 2 * (bigint) gb_a.size() * sz_tagint
                 +     (bigint) gb_type.size() * sz_int
                 + 3 * (bigint) ga1.size() * sz_tagint
                 +     (bigint) ga_type.size() * sz_int
                 + 48 * (bigint) tag_cluster.size();
    utils::logmesg(lmp,
                   "Fix ilves: gathered global topology info for {} bonds, {} angles, "
                   "and {} involved atoms\n"
                   "  replicated topology storage ~ {} bytes/rank ({:.2f} MB)\n",
                   gb_a.size(), ga1.size(), tag_cluster.size(),
                   bytes, (double) bytes / (1024.0 * 1024.0));
  }

  global_topology_ready = true;
}

/* ----------------------------------------------------------------------
   Helper: given two tags, return true if the (ta, tb) bond is selected
   for constraint.  Looks up the bond type via the global topology table.
   Both atoms must be either local or available as ghosts on this rank
   for the type/mass selectors to work.  Returns false otherwise.
------------------------------------------------------------------------- */

bool FixIlves::bond_is_constrained(tagint ta, tagint tb)
{
  tagint tmin = MIN(ta, tb);
  tagint tmax = MAX(ta, tb);
  int lo = 0, hi = (int) gb_a.size();
  while (lo < hi) {
    int mid = (lo + hi) / 2;
    if (gb_a[mid] < tmin || (gb_a[mid] == tmin && gb_b[mid] < tmax)) lo = mid + 1;
    else hi = mid;
  }
  if (lo == (int) gb_a.size() || gb_a[lo] != tmin || gb_b[lo] != tmax) return false;

  int ia = atom->map(ta);
  int ib = atom->map(tb);
  if (ia < 0 || ib < 0) return false;
  return bond_selected_for_atoms(ia, ib, gb_type[lo]);
}

/* ----------------------------------------------------------------------
   Build the flat constraint list from the replicated bond/angle tables.

   - bond constraints: pair (i, partner) for every bond where:
       both atoms are in fix group
       AND bond type is in bond_flag, OR
           either atom type is in type_flag, OR
           either atom mass matches mass_list
   - angle "virtual" A-C constraints: derived from angle entries where
       angle type is in angle_flag
       AND both flanking bonds (A-B, B-C) are also selected (per above)
       AND all three atoms are in the fix group
   - additionally, for near-linear angle types (theta_0 >= linear_threshold),
     a 3-atom B-M virtual constraint is emitted and the bond between B
     and the higher-tag endpoint of {A, C} is dropped.
------------------------------------------------------------------------- */

void FixIlves::build_constraint_list()
{
  const int nlocal_now = atom->nlocal;
  int *mask = atom->mask;
  tagint *tag = atom->tag;

  // mark per-atom flags as we go (zero first)
  for (int i = 0; i < atom->nmax; ++i) ilves_flag[i] = 0;
  n_constr = 0;

  // Lazy gather (in case init_topology was bypassed; normally already done)
  if (!global_topology_ready) gather_global_topology();

  // Identify clusters whose atoms intersect my local atoms.  Each
  // participating rank includes EVERY constraint in any intersecting
  // cluster, even constraints between two ghost atoms, so that all
  // ranks running the redundant cluster solve compute the same
  // Lagrange multipliers from synchronized positions.
  std::unordered_map<tagint, char> my_cluster;
  my_cluster.reserve(nlocal_now);
  for (int i = 0; i < nlocal_now; ++i) {
    tagint t = tag[i];
    auto it = tag_cluster.find(t);
    if (it != tag_cluster.end()) my_cluster[it->second] = 1;
  }

  auto in_my_cluster = [&](tagint t) -> bool {
    auto it = tag_cluster.find(t);
    if (it == tag_cluster.end()) return false;
    return my_cluster.count(it->second) > 0;
  };

  // Pick c_atom1 and c_atom2 indices for a (ta, tb) pair.  c_atom1 is
  // preferentially the locally-owned copy when one of the two is local;
  // ties go to lower tag.  b_out is closest-image to a_out for PBC
  // consistency.
  auto pick_atoms = [&](tagint ta, tagint tb, int &a_out, int &b_out) -> bool {
    int ja = atom->map(ta);
    int jb = atom->map(tb);
    if (ja < 0 || jb < 0) return false;
    bool ja_local = (ja < nlocal_now);
    bool jb_local = (jb < nlocal_now);
    if (ja_local && jb_local) {
      if (tag[ja] < tag[jb]) { a_out = ja; b_out = jb; }
      else                   { a_out = jb; b_out = ja; }
    } else if (ja_local) {
      a_out = ja; b_out = jb;
    } else if (jb_local) {
      a_out = jb; b_out = ja;
    } else {
      // both ghost: include for cluster completion; this rank does not
      // apply force to either atom (apply_constraint_forces tests < nlocal).
      if (tag[ja] < tag[jb]) { a_out = ja; b_out = jb; }
      else                   { a_out = jb; b_out = ja; }
    }
    b_out = domain->closest_image(a_out, b_out);
    return true;
  };

  // -----------------------------------------------------------------
  // Pre-phase: identify bond legs to skip for near-linear angles.
  // -----------------------------------------------------------------
  struct PairHash {
    size_t operator()(const std::pair<tagint,tagint> &p) const noexcept {
      return std::hash<tagint>()(p.first) ^ (std::hash<tagint>()(p.second) << 1);
    }
  };
  std::unordered_set<std::pair<tagint,tagint>, PairHash> dropped_bonds;
  if (has_angle) {
    const int na_global_pre = (int) ga1.size();
    for (int gi = 0; gi < na_global_pre; ++gi) {
      int at = ga_type[gi];
      if (at <= 0 || at > atom->nangletypes) continue;
      if (!angle_flag[at]) continue;
      if (!angle_linear[at]) continue;

      tagint t1 = ga1[gi];     // lower-tag endpoint
      tagint t2 = ga2[gi];     // middle (B)
      tagint t3 = ga3[gi];     // higher-tag endpoint

      // Both flanking bonds must be constrained for the angle to be
      // emitted -- otherwise Phase B wouldn't touch it and the drop set
      // would silently disable a wanted bond.
      if (!bond_is_constrained(t2, t1)) continue;
      if (!bond_is_constrained(t2, t3)) continue;

      // Drop the bond to the higher-tag endpoint (t3).  In gb_a/gb_b
      // canonical order, lo = min(t2, t3), hi = max(t2, t3).
      tagint lo = (t2 < t3) ? t2 : t3;
      tagint hi = (t2 < t3) ? t3 : t2;
      dropped_bonds.emplace(lo, hi);
    }
  }

  // -----------------------------------------------------------------
  // Phase A: bond constraints from clusters I'm in
  // -----------------------------------------------------------------
  const int nb_global = (int) gb_a.size();
  for (int gi = 0; gi < nb_global; ++gi) {
    int bt = gb_type[gi];
    if (bt <= 0 || bt > atom->nbondtypes) continue;
    tagint ta = gb_a[gi];
    tagint tb = gb_b[gi];

    if (!in_my_cluster(ta) && !in_my_cluster(tb)) continue;

    if (!dropped_bonds.empty() && dropped_bonds.count({ta, tb})) continue;

    int a_idx, b_idx;
    if (!pick_atoms(ta, tb, a_idx, b_idx)) continue;

    if (!(mask[a_idx] & groupbit) || !(mask[b_idx] & groupbit)) continue;
    if (!bond_selected_for_atoms(a_idx, b_idx, bt)) continue;

    add_constraint(a_idx, b_idx, bt, bond_distance[bt]);
    ilves_flag[a_idx] = 1;
    ilves_flag[b_idx] = 1;
  }

  // -----------------------------------------------------------------
  // Phase B: angle virtual A-C constraints.  For angle types flagged as
  // near-linear (theta_0 >= linear_threshold) the AC entry would be
  // rank-deficient against the two flanking bonds, so we silently skip
  // it.  negate_constrained_topology() will skip negating those angle
  // types' sign too, leaving the standard angle force-field term active.
  // -----------------------------------------------------------------
  if (has_angle) {
    const int na_global = (int) ga1.size();
    for (int gi = 0; gi < na_global; ++gi) {
      int at = ga_type[gi];
      if (at <= 0 || at > atom->nangletypes) continue;
      if (!angle_flag[at]) continue;
      if (angle_linear[at]) continue;        // angle term handles it

      tagint t1 = ga1[gi];
      tagint t2 = ga2[gi];
      tagint t3 = ga3[gi];

      if (!in_my_cluster(t1) && !in_my_cluster(t2) && !in_my_cluster(t3)) continue;

      int i1 = atom->map(t1);
      int i2 = atom->map(t2);
      int i3 = atom->map(t3);
      if (i1 < 0 || i2 < 0 || i3 < 0) continue;

      if (!(mask[i1] & groupbit)) continue;
      if (!(mask[i2] & groupbit)) continue;
      if (!(mask[i3] & groupbit)) continue;

      if (!bond_is_constrained(t2, t1)) continue;
      if (!bond_is_constrained(t2, t3)) continue;

      int a_idx, b_idx;
      if (!pick_atoms(t1, t3, a_idx, b_idx)) continue;

      add_constraint(a_idx, b_idx, -at, angle_distance[at]);
      ilves_flag[a_idx] = 1;
      ilves_flag[b_idx] = 1;
      ilves_flag[i2]    = 1;
    }
  }

  // connected-component labelling and cluster grouping
  group_by_cluster();
  precompute_constraint_data();
}
