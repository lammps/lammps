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
#include <cstdint>
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
    "author = {L{\\'o}pez-Villellas,  Lori{\\'e}n and Mikkelsen,  Carl Christian Kjelgaard "
    "and Galano-Frutos,  Juan Jos{\\'e} and Marco-Sola,  Santiago and "
    "Alastruey-Bened{\\'e},  Jes{\\'u}s and Ib{\\'a}{\\~n}ez,  Pablo and "
    "Echenique,  Pablo and Moret{\\'o},  Miquel and {De Rosa},  Maria Cristina and "
    "Garc{\\'i}a-Risue{\\~n}o,  Pablo},\n"
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
    angle_r2(nullptr), angle_linear(nullptr), fstore(nullptr), xshake(nullptr), c_atom1(nullptr),
    c_atom2(nullptr), c_type(nullptr), c_dist(nullptr), c_lambda(nullptr), c_cluster(nullptr),
    c_rx(nullptr), c_ry(nullptr), c_rz(nullptr), c_rsq(nullptr), c_invma(nullptr), c_invmb(nullptr),
    cluster_offset(nullptr), c_perm(nullptr), apply_order(nullptr), lu_A(nullptr), lu_b(nullptr),
    lu_pivot(nullptr), cl_sx(nullptr), cl_sy(nullptr), cl_sz(nullptr), chol_pool(nullptr),
    chol_pool_offset(nullptr), cluster_bw(nullptr), cluster_cached(nullptr), x(nullptr), v(nullptr),
    f(nullptr), mass(nullptr), rmass(nullptr), type(nullptr), b_count(nullptr),
    b_count_all(nullptr), b_ave(nullptr), b_max(nullptr), b_min(nullptr), b_ave_all(nullptr),
    b_max_all(nullptr), b_min_all(nullptr), a_count(nullptr), a_count_all(nullptr), a_ave(nullptr),
    a_max(nullptr), a_min(nullptr), a_ave_all(nullptr), a_max_all(nullptr), a_min_all(nullptr)
{
  lu_alloc = 0;
  largest_cluster = 0;
  comm_mode = 0;
  baseline_ready = false;
  natoms_at_init = nconstrbonds_sum_at_init = nconstrangles_sum_at_init = 0;
  angle_btype1 = angle_btype2 = nullptr;
  lang_fbuf = nullptr;
  lang_fbuf_alloc = 0;
  lang_vbuf = nullptr;
  lang_vbuf_alloc = 0;
  cluster_tag = nullptr;
  cluster_tag_alloc = 0;
  cluster_global_id = nullptr;
  cluster_owner_rank = nullptr;
  cluster_id_alloc = 0;
  owner_peers_off = nullptr;
  owner_peer_rank = nullptr;
  owner_peers_alloc_clusters = 0;
  owner_peer_alloc_entries = 0;
  n_owner_peer_entries = 0;
  owned_aug_constr_off = nullptr;
  owned_aug_constr_ta = nullptr;
  owned_aug_constr_tb = nullptr;
  owned_aug_constr_type = nullptr;
  owned_aug_constr_alloc = 0;
  ilv_num_bond = nullptr;
  ilv_bond_atom = nullptr;
  ilv_bond_type = nullptr;
  ilv_bond_per_atom = 0;
  ilv_nmax_alloc = 0;
  // comm_forward = 3 for forward_comm of xshake or v (see pack_forward_comm).
  // comm_reverse = 6 to accommodate the per-atom virial reverse_comm of
  // lang_vbuf (6 doubles per atom).  The lang_fbuf force reverse_comm
  // (3 doubles per atom) fits within the same buffer.  Only
  // active when linearangle K > 0 && newton bond on, but set unconditionally
  // so LAMMPS sizes the comm buffer correctly).
  comm_forward = 3;
  comm_reverse = 6;
  chol_pool_alloc = 0;
  chol_offset_alloc = 0;
  cluster_bw_alloc = 0;
  apply_order_alloc = 0;
  cluster_cached_alloc = 0;
  energy_global_flag = energy_peratom_flag = 1;
  virial_global_flag = virial_peratom_flag = 1;
  thermo_energy = thermo_virial = 1;
  dof_flag = 1;
  scalar_flag = 1;
  extscalar = 1;
  next_output = -1;

  eflag_pre_reverse = 0;
  ebond = 0.0;
  has_angle = false;
  variant = ILVES_FAST;
  chol_calls = 0;
  chol_fallbacks = 0;
  newton_iter_sum = 0;
  newton_iter_max = 0;
  newton_solve_count = 0;
  linear_threshold = 165.0;
  linear_angle_K = 0.0;

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

  // forward-comm payload size per atom.
  //   Mode 0 (xshake) / mode 1 (atom->v): 3 doubles.
  //   Mode 2 (ghost bond data refresh): 1 + 2*bond_per_atom values.
  // bond_per_atom is set by the atom_style + read_data before the fix is
  // created, so we know the size here.  Set comm_forward to the max of
  // the two so LAMMPS allocates the comm buffer correctly at first use.
  ilv_bond_per_atom = atom->bond_per_atom;
  const int bond_data_payload = 1 + 2 * ilv_bond_per_atom;
  comm_forward = (bond_data_payload > 3) ? bond_data_payload : 3;

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
      if ((i < 1) || (i > atom->nbondtypes))
        error->all(FLERR, next, "Invalid bond type {} for fix ilves", arg[next]);
      bond_flag[i] = 1;

    } else if (mode == 'a') {
      if (allow_typelabels) i = utils::expand_type_int(FLERR, arg[next], Atom::ANGLE, lmp);
      else                  i = utils::inumeric(FLERR, arg[next], false, lmp);
      if ((i < 1) || (i > atom->nangletypes))
        error->all(FLERR, next, "Invalid angle type {} for fix ilves", arg[next]);
      angle_flag[i] = 1;
      has_angle = true;

    } else if (mode == 't') {
      if (allow_typelabels) i = utils::expand_type_int(FLERR, arg[next], Atom::ATOM, lmp);
      else                  i = utils::inumeric(FLERR, arg[next], false, lmp);
      if ((i < 1) || (i > atom->ntypes))
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
      // optional K argument: when > 0, fix ilves applies a stiff angle
      // potential E = K*(1 + cos(theta)) to near-linear angles, in place
      // of the user's angle_style.  When 0 (default), near-linear angles
      // are left to the user's angle_style.
      if (next < narg) {
        char *endp = nullptr;
        double v = strtod(arg[next], &endp);
        if (endp != arg[next] && *endp == '\0') {
          if (v < 0.0)
            error->all(FLERR, next, "Fix ilves linearangle K must be >= 0.0");
          linear_angle_K = v;
          ++next;
        }
      }

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
  angle_btype1   = new int[atom->nangletypes + 1]{};
  angle_btype2   = new int[atom->nangletypes + 1]{};

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
  delete[] angle_btype1;
  delete[] angle_btype2;

  delete[] b_count;     delete[] b_count_all;
  delete[] b_ave;       delete[] b_ave_all;
  delete[] b_max;       delete[] b_max_all;
  delete[] b_min;       delete[] b_min_all;
  delete[] a_count;     delete[] a_count_all;
  delete[] a_ave;       delete[] a_ave_all;
  delete[] a_max;       delete[] a_max_all;
  delete[] a_min;       delete[] a_min_all;

  memory->destroy(xshake);
  memory->destroy(fstore);
  memory->destroy(lang_fbuf);
  memory->destroy(lang_vbuf);
  memory->destroy(cluster_tag);
  memory->destroy(cluster_global_id);
  memory->destroy(cluster_owner_rank);
  memory->destroy(owner_peers_off);
  memory->destroy(owner_peer_rank);
  memory->destroy(owned_aug_constr_off);
  memory->destroy(owned_aug_constr_ta);
  memory->destroy(owned_aug_constr_tb);
  memory->destroy(owned_aug_constr_type);
  memory->destroy(ilv_num_bond);
  memory->destroy(ilv_bond_atom);
  memory->destroy(ilv_bond_type);

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
  memory->destroy(apply_order);
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

  // Newton on or off for bond is supported by the constraint solver.
  // Under newton off bond, each bond is stored at both endpoints;
  // lookup_local_bond_type can read the type from either side.  Under
  // newton on bond, the bond is stored only at the lower-tag endpoint;
  // when both endpoints are ghosts on this rank, lookup_local_bond_type
  // returns 0 and the angle's flanking-bond check falls back to the
  // angle_btype1/btype2 cache built by init_topology (which uses an
  // MPI_Allreduce(MAX) and so works as long as ANY rank globally can
  // resolve each angle type's flanking types).
  //
  // Stiff angle force (linearangle K > 0) supports newton on via a
  // dedicated lang_fbuf force buffer + custom reverse_comm, and
  // lang_vbuf for per-atom virial when stress/atom is in use.  Set
  // comm_reverse = 6 so LAMMPS's comm allocator sizes the buffer for
  // the larger of the two payloads (lang_vbuf at 6 doubles/atom).
  if (linear_angle_K > 0.0 && force->newton_bond != 0) {
    comm_reverse = 6;
  }

  // ilv_bond_per_atom and comm_forward are set in the constructor; the
  // Schwarz overlap ghost-bond-data forward_comm uses comm_mode == 2.

  for (int i = 1; i <= atom->nbondtypes; ++i)
    bond_distance[i] = force->bond->equilibrium_distance(i);

  // sanity check: if angle constraints are requested, an angle_style must
  // be defined (the variant-specific init_topology() below uses
  // force->angle->equilibrium_angle to compute angle_distance[]).
  if (has_angle && force->angle == nullptr)
    error->all(FLERR, Error::NOLASTLINE,
               "Angle style must be defined for fix ilves with angle constraints");

  // Topology setup: compute equilibrium AC distance for every selected
  // angle type and consensus-check its flanking bond types across ranks.
  // No global bond/angle gather -- each rank walks its own local angle
  // storage and reduces by MPI_Allreduce(MAX) over the per-type pair.
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
    if (n_linear && (comm->me == 0)) {
      if (linear_angle_K > 0.0)
        utils::logmesg(lmp,
                       "Fix ilves: skipping AC virtual constraint for {} near-linear angle "
                       "type(s) (theta_0 >= {} deg); using a cosine/delta angle potential with "
                       "K_theta = {} instead\n", n_linear, linear_threshold, linear_angle_K);
      else
        utils::logmesg(lmp,
                       "Fix ilves: skipping AC virtual constraint for {} near-linear angle "
                       "type(s) (theta_0 >= {} deg); the original angle force-field term "
                       "remains active for these angles\n", n_linear, linear_threshold);
    }
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
        // Near-linear angle: only negate if the stiff angle force is
        // active (linearangle K > 0).  When K == 0 we leave the angle
        // positive so the user's angle_style still acts on it.
        if (angle_linear[at] && linear_angle_K <= 0.0) continue;
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

  // count how many local clusters this rank OWNS vs how many are
  // remotely owned (single-rank cluster ownership groundwork).
  int n_owned = 0, n_remote = 0;
  for (int c = 0; c < n_clusters; ++c) {
    if (cluster_owner_rank[c] == comm->me) ++n_owned;
    else ++n_remote;
  }
  int owned_info[2] = {n_owned, n_remote};
  int owned_info_all[2];
  MPI_Reduce(owned_info, owned_info_all, 2, MPI_INT, MPI_SUM, 0, world);

  // Comm-graph extent: number of (cluster, peer-rank) entries this rank
  // owns plus the augmented-constraint count contributed by peers.
  int my_peer_entries = n_owner_peer_entries;
  int my_aug_constr = (n_clusters > 0) ? owned_aug_constr_off[n_clusters] : 0;
  int peer_info[2] = {my_peer_entries, my_aug_constr};
  int peer_info_all[2];
  MPI_Reduce(peer_info, peer_info_all, 2, MPI_INT, MPI_SUM, 0, world);

  if (comm->me == 0)
    utils::logmesg(lmp, "Fix ilves info:\n"
                   "  {} algorithm\n"
                   "  {} constraints/proc\n"
                   "  {} clusters/proc\n"
                   "  {} max. constraints/cluster\n"
                   "  {} owned + {} remote clusters (global)\n"
                   "  {} owner-peer entries, {} augmented constraints (global)\n",
                   (variant == ILVES_FULL) ? "ILVES (full, LU)" : "ILVES (fast, banded Cholesky)",
                   n_info_all[0], n_info_all[1], n_max_all,
                   owned_info_all[0], owned_info_all[1],
                   peer_info_all[0], peer_info_all[1]);

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

/* ----------------------------------------------------------------------
   Re-run topology negation now that ghosts are populated.  The init()
   call to negate_constrained_topology() runs before comm->borders, so
   atom->map(partner) returns -1 for any cross-rank bond/angle partner
   and that constraint's storage stays positive -- which leaks past the
   PARTIAL bondlist filter and causes the user's bond/angle style to
   double-count force at setup.  setup_pre_neighbor() runs after
   comm->borders but before neighbor->build, so the bondlist that comes
   out of this setup pass correctly excludes every constrained bond.
   Idempotent: the bt <= 0 / at <= 0 guards inside the negate routine
   skip entries already flipped by init().
------------------------------------------------------------------------- */

void FixIlves::setup_pre_neighbor()
{
  negate_constrained_topology();
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
  // per-rank constraint list built at init.  Skip on the first call
  // before init_topology has run (baseline_ready is the gate).
  if (baseline_ready) check_topology_unchanged();
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
   constraints belonging to the same cluster are contiguous in c_perm.
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
   Refresh c_rx/c_ry/c_rz/c_rsq from the current positions in x[].  Uses
   minimum_image to wrap the displacement to the closest periodic image
   (closest_image can return inconsistent ghost copies for atoms whose
   periodic images differ, so we operate on the difference directly).
   Called at the top of every solver entry point (post_force,
   post_force_respa, end_of_step) and from precompute_constraint_data
   inside build_constraint_list.
------------------------------------------------------------------------- */

void FixIlves::refresh_constraint_geometry()
{
  for (int k = 0; k < n_constr; ++k) {
    const int a = c_atom1[k];
    const int b = c_atom2[k];
    double rx = x[a][0] - x[b][0];
    double ry = x[a][1] - x[b][1];
    double rz = x[a][2] - x[b][2];
    domain->minimum_image(FLERR, rx, ry, rz);
    c_rx[k] = rx;
    c_ry[k] = ry;
    c_rz[k] = rz;
    c_rsq[k] = rx*rx + ry*ry + rz*rz;
  }
}

/* ----------------------------------------------------------------------
   Precompute per-constraint reference vectors, |r_k|^2, and inverse
   masses.  Used by all subsequent solver entries until the next
   reneighbor.  RCM reorder + factor cache sizing also run here since
   the per-cluster bandwidths depend on the constraint list, which is
   just-finalized at this point in build_constraint_list.
------------------------------------------------------------------------- */

void FixIlves::precompute_constraint_data()
{
  refresh_constraint_geometry();
  for (int k = 0; k < n_constr; ++k) {
    const int a = c_atom1[k];
    const int b = c_atom2[k];
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
   Build the tag-pair-sorted apply order used by apply_constraint_forces.
   Sorted by (min_tag, max_tag, type) so different MPI ranks accumulate
   force / per-atom virial contributions in the same order on shared
   atoms -- eliminates ULP-level rank-dependence in the per-atom output
   (and the Newton iteration that follows from it).  The constraint list
   is stable until the next reneighbor, so this runs once per neighbor
   list build instead of every solve.
------------------------------------------------------------------------- */

void FixIlves::build_apply_order()
{
  if (n_constr > apply_order_alloc) {
    apply_order_alloc = n_constr;
    memory->grow(apply_order, apply_order_alloc, "ilves:apply_order");
  }
  for (int k = 0; k < n_constr; ++k) apply_order[k] = k;
  tagint *tag = atom->tag;
  std::sort(apply_order, apply_order + n_constr, [&](int x, int y) {
    tagint xa = tag[c_atom1[x]];
    tagint xb = tag[c_atom2[x]];
    tagint ya = tag[c_atom1[y]];
    tagint yb = tag[c_atom2[y]];
    if (xa > xb) std::swap(xa, xb);
    if (ya > yb) std::swap(ya, yb);
    if (xa != ya) return xa < ya;
    if (xb != yb) return xb < yb;
    return c_type[x] < c_type[y];
  });
}

/* ----------------------------------------------------------------------
   Per-cluster reverse Cuthill-McKee on the constraint-adjacency graph
   (two constraints share a graph edge iff they share an atom).  Applies
   the resulting permutation in-place to c_perm, then computes each
   cluster's bandwidth (max over edges (k,l) with k<l of slot_l-slot_k
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
  if (n_clusters == 0) return;

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

    // apply permutation to c_perm[beg..end)
    old_perm.assign(n_c, 0);
    for (int s = 0; s < n_c; ++s) old_perm[s] = c_perm[beg + s];
    for (int s = 0; s < n_c; ++s) {
      int new_slot_for_old_s = old_to_new[s];
      c_perm[beg + new_slot_for_old_s] = old_perm[s];
    }

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

  refresh_constraint_geometry();

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

  // optional stiff angle force for near-linear angles (when linearangle K > 0)
  apply_linear_angle_forces(vflag);
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

  refresh_constraint_geometry();

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
  apply_linear_angle_forces(vflag);
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
     - assemble the local n_c x n_c Jacobian A_c (variant=full uses the
       exact-Newton asymmetric form; variant=fast uses the symmetric
       quasi-Newton form whose factor is cached across iterations)
     - assemble the rhs as the current constraint residuals g_k
     - solve A_c * d-lambda = -g_c
     - update xshake using d-lambda
     - accumulate c_lambda[k] += d-lambda[k]
   Iterate until the globally-reduced max relative bond-length
   violation falls below tolerance, or max_iter is reached.  At
   multi-rank the per-iteration forward_comm(xshake) propagates
   inter-cluster updates between owner and peer ranks; the
   MPI_Allreduce(MAX) on the residual keeps every rank's loop in
   lockstep.  Returns true if converged.
------------------------------------------------------------------------- */

bool FixIlves::solve_constraints()
{
  // NOTE: do NOT early-return on n_constr == 0.  Under distributed
  // topology some ranks may have zero local constraints while others
  // have many.  The forward_comm and MPI_Allreduce in the Newton
  // iteration loop are collective, so all ranks must reach them in
  // lockstep -- skipping the loop on empty ranks deadlocks the others.

  grow_lu_workspace(largest_cluster);

  // Cold-start: zero c_lambda for every constraint.  Per-iter updates
  // inside the Newton loop accumulate into these.
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

  int iter = 0;
  for (; iter < max_iter; ++iter) {
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
      // ghost (owner on another rank), atom->map returns the ghost
      // index -- skip the partner update locally; the partner rank
      // holds the same constraint and will update its local atom, and
      // forward_comm at the end of the iteration brings the corrected
      // ghost position back to this rank.
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

  // Track Newton iteration counts for stats reporting.  iter is the loop
  // counter at exit; convert to 1-based count (iter+1 if we broke out on
  // convergence, max_iter if we hit the cap).
  const bigint iters_used = (bigint)(converged ? iter + 1 : max_iter);
  newton_iter_sum   += iters_used;
  if (iters_used > newton_iter_max) newton_iter_max = iters_used;
  newton_solve_count++;

  if (!converged && (comm->me == 0)) {
    error->warning(FLERR,
                   "Fix ilves: max_iter ({}) reached without reaching tol={} at step {}",
                   max_iter, tolerance, update->ntimestep);
  }

  return converged;
}

/* ----------------------------------------------------------------------
   Convert the accumulated Lagrange multipliers into constraint forces
   added to atom->f, so that the next initial_integrate puts the atoms
   at the constrained xshake positions.

   For each constraint k = (a, b):
     f[a] += (lambda_k / dtfsq) * r_k
     f[b] -= (lambda_k / dtfsq) * r_k
   Forces are applied only to atoms that are local (< nlocal).  Cross-
   rank constraints are owned by every participating rank after the
   owner-to-peers broadcast in build_constraint_list, so each rank's
   apply pass writes Newton-3rd-law-symmetric contributions to its own
   local atoms.
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

  // Iterate in canonical (tag-pair-sorted) apply_order built at the
  // last reneighbor.  Different ranks may have stored a given shared
  // constraint at different positions in c_atom1/c_atom2; iterating
  // raw 0..n_constr-1 mixes writes to the same atom's f / vatom in
  // different orders across ranks, producing ULP-level diffs that
  // propagate to ULP-level x_projected differences in
  // correct_coordinates, then through Newton iteration to several
  // percent c_lambda differences on small multipliers.  A stable
  // canonical order eliminates the rank-dependence and lets per-atom
  // virial match np=1 bit-for-bit on owned atoms.
  for (int idx = 0; idx < n_constr; ++idx) {
    const int k = apply_order[idx];
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
   Apply a stiff angle force to angle types flagged as near-linear, using
   the cosine form E = K*(1 + cos(theta)).  Force formula is identical to
   LAMMPS's angle_style cosine (see angle_cosine.cpp); the 1/sin(theta)
   factor that breaks angle_style harmonic at theta=180 deg cancels out
   here because the energy is a polynomial in cos(theta), not theta.

   Activated by `linearangle <theta_deg> <K>` with K > 0.  When active,
   negate_constrained_topology() also negates the near-linear angle types
   so the user's angle_style->compute skips them -- avoiding double
   counting between the user's angle force and this one.

   Iteration model under newton off bond: each local atom of the angle
   has the angle in its storage.  We iterate every local atom and, for
   each angle entry, compute the full force trio (f1, f2, f3) but apply
   ONLY the share for the local atom -- duplicating the computation
   three times globally but avoiding any reverse_comm.
------------------------------------------------------------------------- */

void FixIlves::apply_linear_angle_forces(int vflag)
{
  if (linear_angle_K <= 0.0) return;
  if (!has_angle) return;

  const int newton_bond = force->newton_bond;
  const int nlocal_now = atom->nlocal;
  const int nmax_now   = atom->nmax;
  tagint *tag         = atom->tag;
  int *mask           = atom->mask;
  int *num_angle      = atom->num_angle;
  int **na_type       = atom->angle_type;
  tagint **na_atom1   = atom->angle_atom1;
  tagint **na_atom2   = atom->angle_atom2;
  tagint **na_atom3   = atom->angle_atom3;
  double **xc         = atom->x;
  double **fc         = atom->f;
  const int natypes   = atom->nangletypes;

  // Under newton on, accumulate force into our own per-atom buffer (sized
  // nmax * 3) so we can reverse_comm ghost contributions to owners without
  // disturbing atom->f's ghost values (which other code may rely on).
  // When vflag_atom is active, also allocate and zero lang_vbuf for the
  // per-atom virial reverse_comm.
  const bool need_vatom = (vflag_atom && newton_bond);
  if (newton_bond) {
    grow_lang_fbuf(nmax_now);
    const bigint nb = 3 * (bigint) nmax_now;
    for (bigint k = 0; k < nb; ++k) lang_fbuf[k] = 0.0;
  }
  if (need_vatom) {
    grow_lang_vbuf(nmax_now);
    const bigint nbv = 6 * (bigint) nmax_now;
    for (bigint k = 0; k < nbv; ++k) lang_vbuf[k] = 0.0;
  }

  for (int i = 0; i < nlocal_now; ++i) {
    if (!(mask[i] & groupbit)) continue;
    const tagint ti = tag[i];
    const int nang = num_angle[i];
    for (int m = 0; m < nang; ++m) {
      int at = na_type[i][m];
      if (at == 0) continue;
      // The angle_type was negated for near-linear constrained angles
      // (see negate_constrained_topology); accept either sign for the
      // angle-type lookup.
      const int at_abs = (at < 0) ? -at : at;
      if (at_abs > natypes) continue;
      if (!angle_flag[at_abs]) continue;
      if (!angle_linear[at_abs]) continue;

      // Endpoints and middle.  Resolve all three to local indices in
      // this rank's halo.  Skip if any is not reachable.
      const tagint t1 = na_atom1[i][m];
      const tagint t2 = na_atom2[i][m];
      const tagint t3 = na_atom3[i][m];
      int i1 = atom->map(t1);
      int i2 = atom->map(t2);
      int i3 = atom->map(t3);
      if (i1 < 0 || i2 < 0 || i3 < 0) continue;

      // Cosine angle force (same math as LAMMPS angle_style cosine):
      //   E = K * (1 + cos(theta))
      //   c = cos(theta) = (r1 . r2) / (|r1| |r2|)
      //   force on atom 1 (= ti) = -dE/dx1 expanded via chain rule
      double dx1 = xc[i1][0] - xc[i2][0];
      double dy1 = xc[i1][1] - xc[i2][1];
      double dz1 = xc[i1][2] - xc[i2][2];
      domain->minimum_image(FLERR, dx1, dy1, dz1);
      const double rsq1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
      const double r1   = sqrt(rsq1);

      double dx2 = xc[i3][0] - xc[i2][0];
      double dy2 = xc[i3][1] - xc[i2][1];
      double dz2 = xc[i3][2] - xc[i2][2];
      domain->minimum_image(FLERR, dx2, dy2, dz2);
      const double rsq2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
      const double r2   = sqrt(rsq2);

      if (r1 <= 0.0 || r2 <= 0.0) continue;

      double c = dx1*dx2 + dy1*dy2 + dz1*dz2;
      c /= r1*r2;
      if (c >  1.0) c =  1.0;
      if (c < -1.0) c = -1.0;

      const double a    = linear_angle_K;
      const double a11  = a*c / rsq1;
      const double a12  = -a / (r1*r2);
      const double a22  = a*c / rsq2;

      double f1[3], f3[3];
      f1[0] = a11*dx1 + a12*dx2;
      f1[1] = a11*dy1 + a12*dy2;
      f1[2] = a11*dz1 + a12*dz2;
      f3[0] = a22*dx2 + a12*dx1;
      f3[1] = a22*dy2 + a12*dy1;
      f3[2] = a22*dz2 + a12*dz1;

      if (newton_bond) {
        // newton on bond: angle stored only at middle.  Only the
        // middle-owner processes; write force trio into lang_fbuf for
        // ALL three atoms (some may be ghosts).  reverse_comm later
        // sums ghost contributions into local owners.
        if (ti != t2) continue;
        lang_fbuf[3*i1+0] += f1[0]; lang_fbuf[3*i1+1] += f1[1]; lang_fbuf[3*i1+2] += f1[2];
        lang_fbuf[3*i2+0] -= f1[0] + f3[0];
        lang_fbuf[3*i2+1] -= f1[1] + f3[1];
        lang_fbuf[3*i2+2] -= f1[2] + f3[2];
        lang_fbuf[3*i3+0] += f3[0]; lang_fbuf[3*i3+1] += f3[1]; lang_fbuf[3*i3+2] += f3[2];
      } else {
        // newton off bond: angle stored at all three atoms.  Each rank
        // with the angle locally stored applies only its own atom's
        // share to atom->f directly (no reverse_comm needed).
        if (ti == t1) {
          fc[i1][0] += f1[0]; fc[i1][1] += f1[1]; fc[i1][2] += f1[2];
        } else if (ti == t2) {
          fc[i2][0] -= f1[0] + f3[0];
          fc[i2][1] -= f1[1] + f3[1];
          fc[i2][2] -= f1[2] + f3[2];
        } else if (ti == t3) {
          fc[i3][0] += f3[0]; fc[i3][1] += f3[1]; fc[i3][2] += f3[2];
        }
      }

      // Energy and global virial tally.  Pattern follows Angle::ev_tally:
      //   newton on bond:  only middle iterates (storage convention).
      //     Middle's rank adds FULL energy/virial in one shot.
      //   newton off bond: angle stored at all three atoms; each local
      //     atom owner adds its 1/3 share.  Three ranks each add 1/3,
      //     sum globally = full.
      const double eangle = a * (1.0 + c);
      double bond_v[6];
      bond_v[0] = dx1*f1[0] + dx2*f3[0];
      bond_v[1] = dy1*f1[1] + dy2*f3[1];
      bond_v[2] = dz1*f1[2] + dz2*f3[2];
      bond_v[3] = dx1*f1[1] + dx2*f3[1];
      bond_v[4] = dx1*f1[2] + dx2*f3[2];
      bond_v[5] = dy1*f1[2] + dy2*f3[2];
      if (newton_bond) {
        if (ti == t2) {
          ebond += eangle;
          if (evflag && vflag_global) {
            for (int kk = 0; kk < 6; ++kk) virial[kk] += bond_v[kk];
          }
          // Per-atom virial under newton on bond: the middle rank
          // sees the angle's full virial; mirror Angle::ev_tally and
          // write v/3 to each of i1, i2, i3 (some may be ghosts on
          // this rank; reverse_comm of lang_vbuf below delivers
          // ghost contributions to the owners).
          if (evflag && vflag_atom) {
            const double third = 1.0 / 3.0;
            const double v0 = third * bond_v[0];
            const double v1 = third * bond_v[1];
            const double v2 = third * bond_v[2];
            const double v3 = third * bond_v[3];
            const double v4 = third * bond_v[4];
            const double v5 = third * bond_v[5];
            lang_vbuf[6*i1+0] += v0; lang_vbuf[6*i1+1] += v1;
            lang_vbuf[6*i1+2] += v2; lang_vbuf[6*i1+3] += v3;
            lang_vbuf[6*i1+4] += v4; lang_vbuf[6*i1+5] += v5;
            lang_vbuf[6*i2+0] += v0; lang_vbuf[6*i2+1] += v1;
            lang_vbuf[6*i2+2] += v2; lang_vbuf[6*i2+3] += v3;
            lang_vbuf[6*i2+4] += v4; lang_vbuf[6*i2+5] += v5;
            lang_vbuf[6*i3+0] += v0; lang_vbuf[6*i3+1] += v1;
            lang_vbuf[6*i3+2] += v2; lang_vbuf[6*i3+3] += v3;
            lang_vbuf[6*i3+4] += v4; lang_vbuf[6*i3+5] += v5;
          }
        }
      } else {
        ebond += eangle * (1.0 / 3.0);
        if (evflag && vflag_global) {
          for (int kk = 0; kk < 6; ++kk) virial[kk] += bond_v[kk] * (1.0 / 3.0);
        }
        // Per-atom virial under newton off bond: each of the three
        // ranks that has one of t1/t2/t3 locally processes the angle
        // once (selected by `ti == tX` below) and writes v/3 to ITS
        // OWN local atom's vatom.  Sum across the three ranks = v.
        if (evflag && vflag_atom) {
          int local_idx = -1;
          if      (ti == t1) local_idx = i1;
          else if (ti == t2) local_idx = i2;
          else if (ti == t3) local_idx = i3;
          if (local_idx >= 0 && local_idx < nlocal_now) {
            const double third = 1.0 / 3.0;
            vatom[local_idx][0] += third * bond_v[0];
            vatom[local_idx][1] += third * bond_v[1];
            vatom[local_idx][2] += third * bond_v[2];
            vatom[local_idx][3] += third * bond_v[3];
            vatom[local_idx][4] += third * bond_v[4];
            vatom[local_idx][5] += third * bond_v[5];
          }
        }
      }
    }
  }

  if (newton_bond) {
    // reverse_comm: sum ghost lang_fbuf contributions into the local
    // entries on the owning ranks.  After this, lang_fbuf[i] for local
    // i contains the full angle force contribution for atom i.
    comm_mode = 3;
    comm->reverse_comm(this);
    comm_mode = 0;

    // Apply lang_fbuf to atom->f for local atoms.
    for (int i = 0; i < nlocal_now; ++i) {
      fc[i][0] += lang_fbuf[3*i+0];
      fc[i][1] += lang_fbuf[3*i+1];
      fc[i][2] += lang_fbuf[3*i+2];
    }

    // Per-atom virial: deliver ghost-slot lang_vbuf contributions to
    // their local owners and then accumulate into fix's vatom.  Each
    // atom's per-angle share was already divided by 3 above; we just
    // need to sum the ghost-side writes into the local atom.
    if (need_vatom) {
      comm_mode = 6;
      comm->reverse_comm(this);
      comm_mode = 0;
      for (int i = 0; i < nlocal_now; ++i) {
        vatom[i][0] += lang_vbuf[6*i+0];
        vatom[i][1] += lang_vbuf[6*i+1];
        vatom[i][2] += lang_vbuf[6*i+2];
        vatom[i][3] += lang_vbuf[6*i+3];
        vatom[i][4] += lang_vbuf[6*i+4];
        vatom[i][5] += lang_vbuf[6*i+5];
      }
    }
  }
  (void) vflag;
}

/* ---------------------------------------------------------------------- */

void FixIlves::grow_lang_fbuf(int nmax)
{
  if (nmax <= lang_fbuf_alloc) return;
  memory->grow(lang_fbuf, (bigint) nmax * 3, "ilves:lang_fbuf");
  lang_fbuf_alloc = nmax;
}

void FixIlves::grow_lang_vbuf(int nmax)
{
  if (nmax <= lang_vbuf_alloc) return;
  memory->grow(lang_vbuf, (bigint) nmax * 6, "ilves:lang_vbuf");
  lang_vbuf_alloc = nmax;
}

void FixIlves::grow_cluster_tag(int nmax)
{
  if (nmax <= cluster_tag_alloc) return;
  memory->grow(cluster_tag, nmax, "ilves:cluster_tag");
  cluster_tag_alloc = nmax;
}

void FixIlves::grow_cluster_id(int n)
{
  if (n <= cluster_id_alloc) return;
  memory->grow(cluster_global_id, n, "ilves:cluster_global_id");
  memory->grow(cluster_owner_rank, n, "ilves:cluster_owner_rank");
  cluster_id_alloc = n;
}

/* ----------------------------------------------------------------------
   Compute per-local-cluster global ID (= smallest tag of any atom in the
   cluster across ALL ranks) and the cluster owner (rank holding that
   smallest-tag atom locally).  Algorithm:

     1) cluster_tag[i] = tag[i] for all atoms.  For atoms in any local
        cluster, replace with the min tag seen in that cluster's local
        view.
     2) Iterate: forward_comm cluster_tag (comm_mode == 5) so ghost
        atoms carry the OWNER rank's cluster_tag.  Walk local clusters;
        for each, if any ghost atom in the cluster has lower
        cluster_tag, lower the entire cluster's cluster_tag.  Repeat
        until no rank made a change (Allreduce(MAX) on a dirty flag).
     3) cluster_global_id[c] = cluster_tag at any atom of cluster c.
        cluster_owner_rank[c] = comm->me if atom->map(global_id) is
        local on this rank; otherwise we send a probe to discover the
        owner.  (Initial impl: defer owner-rank lookup to a separate
        pass via Allgather of "I locally own these cluster ids".)

   This function is called once per reneighbor, from build_constraint_list
   after group_by_cluster.  The owner/peer assignment it produces drives
   build_comm_graph and the per-cluster Alltoallv exchanges that follow.
------------------------------------------------------------------------- */

void FixIlves::identify_cluster_owners()
{
  const int nmax = atom->nmax;
  const int nlocal_now = atom->nlocal;
  tagint *tag = atom->tag;

  grow_cluster_tag(nmax);

  // Initialize cluster_tag[i] = tag[i] (atoms not in any cluster keep
  // their own tag).  For ghost slots we set to 0; forward_comm will
  // refresh them.
  for (int i = 0; i < nlocal_now; ++i) cluster_tag[i] = tag[i];
  for (int i = nlocal_now; i < nmax; ++i) cluster_tag[i] = 0;

  // For each local cluster c, compute the min tag across the atoms it
  // touches (looking at c_atom1 and c_atom2 of every constraint).  Map
  // atom_index -> cluster_idx via union-find result already in c_cluster.
  if (n_clusters > 0) {
    grow_cluster_id(n_clusters);
    for (int c = 0; c < n_clusters; ++c) cluster_global_id[c] = MAXTAGINT;

    // First pass: gather min tag per local cluster.
    for (int k = 0; k < n_constr; ++k) {
      const int c = c_cluster[k];
      const tagint ta = tag[c_atom1[k]];
      const tagint tb = tag[c_atom2[k]];
      if (ta > 0 && ta < cluster_global_id[c]) cluster_global_id[c] = ta;
      if (tb > 0 && tb < cluster_global_id[c]) cluster_global_id[c] = tb;
    }

    // Stamp cluster_tag at every CANONICAL atom slot touched by any
    // constraint.  c_atom2[k] may be a NON-canonical PBC ghost slot
    // (closest_image picks the spatially closest periodic image), but
    // forward_comm propagates cluster_tag from the canonical local
    // (or owner-local) slot.  If we stamp only the non-canonical
    // ghost slot, the next forward_comm overwrites it with the
    // canonical's unstamped value (= the atom's own tag), and the
    // re-stamp loop below will keep oscillating.  Route both endpoint
    // writes through atom->map(tag).
    for (int k = 0; k < n_constr; ++k) {
      const int c = c_cluster[k];
      const tagint gid = cluster_global_id[c];
      const int a_can = atom->map(atom->tag[c_atom1[k]]);
      const int b_can = atom->map(atom->tag[c_atom2[k]]);
      if (a_can >= 0 && (cluster_tag[a_can] == 0 || gid < cluster_tag[a_can]))
        cluster_tag[a_can] = gid;
      if (b_can >= 0 && (cluster_tag[b_can] == 0 || gid < cluster_tag[b_can]))
        cluster_tag[b_can] = gid;
    }
  }

  // Iterate: forward_comm cluster_tag so ghost atoms carry the owner's
  // cluster_tag, then walk local clusters and reduce.
  constexpr int max_passes = 10;
  int pass = 0;
  for (; pass < max_passes; ++pass) {
    comm_mode = 5;
    comm->forward_comm(this);
    comm_mode = 0;

    int dirty_local = 0;

    if (n_clusters > 0) {
      // For each cluster, scan its constraint atoms.  Read cluster_tag
      // at the CANONICAL slot for each tag (atom->map) so we see
      // owner-propagated values, not stale non-canonical PBC ghost
      // slots that forward_comm may have left untouched.
      for (int k = 0; k < n_constr; ++k) {
        const int c = c_cluster[k];
        const int a_can = atom->map(atom->tag[c_atom1[k]]);
        const int b_can = atom->map(atom->tag[c_atom2[k]]);
        const tagint ta = (a_can >= 0) ? cluster_tag[a_can] : 0;
        const tagint tb = (b_can >= 0) ? cluster_tag[b_can] : 0;
        tagint best = cluster_global_id[c];
        if (ta > 0 && ta < best) best = ta;
        if (tb > 0 && tb < best) best = tb;
        if (best < cluster_global_id[c]) {
          cluster_global_id[c] = best;
          dirty_local = 1;
        }
      }
      // Re-stamp cluster_tag at the CANONICAL local slot of every
      // constraint endpoint.  Writing the canonical is what
      // forward_comm later propagates to all of that tag's ghosts.
      for (int k = 0; k < n_constr; ++k) {
        const int c = c_cluster[k];
        const tagint gid = cluster_global_id[c];
        const int a_can = atom->map(atom->tag[c_atom1[k]]);
        const int b_can = atom->map(atom->tag[c_atom2[k]]);
        if (a_can >= 0 && cluster_tag[a_can] > gid) {
          cluster_tag[a_can] = gid;
          dirty_local = 1;
        }
        if (b_can >= 0 && cluster_tag[b_can] > gid) {
          cluster_tag[b_can] = gid;
          dirty_local = 1;
        }
      }
    }

    int dirty_global = 0;
    if (comm->nprocs > 1)
      MPI_Allreduce(&dirty_local, &dirty_global, 1, MPI_INT, MPI_MAX, world);
    else
      dirty_global = dirty_local;
    if (!dirty_global) break;
  }
  if (pass == max_passes && comm->me == 0)
    error->warning(FLERR,
                   "Fix ilves: cluster-id propagation did not converge in {} passes",
                   max_passes);

  // Determine cluster owner.  For each local cluster c, the owner is
  // the rank whose local atom has tag == cluster_global_id[c].  On this
  // rank we own iff atom->map(global_id) returns a local index.
  for (int c = 0; c < n_clusters; ++c) {
    const tagint gid = cluster_global_id[c];
    const int i = atom->map(gid);
    if (i >= 0 && i < nlocal_now) cluster_owner_rank[c] = comm->me;
    else                          cluster_owner_rank[c] = -1;  // peer, owner TBD
  }

  // Owner discovery for peer clusters: each rank Allgather's the list of
  // global IDs it owns; every rank then knows the owner of every global
  // cluster it participates in.  O(nprocs * n_clusters_owned) comm.
  if (comm->nprocs > 1 && n_clusters > 0) {
    std::vector<tagint> my_owned;
    my_owned.reserve(n_clusters);
    for (int c = 0; c < n_clusters; ++c)
      if (cluster_owner_rank[c] == comm->me)
        my_owned.push_back(cluster_global_id[c]);
    int my_n = (int) my_owned.size();

    std::vector<int> counts(comm->nprocs, 0), displs(comm->nprocs, 0);
    MPI_Allgather(&my_n, 1, MPI_INT, counts.data(), 1, MPI_INT, world);
    int total = 0;
    for (int p = 0; p < comm->nprocs; ++p) {
      displs[p] = total;
      total += counts[p];
    }

    std::vector<tagint> all_ids(total > 0 ? total : 1);
    MPI_Allgatherv(my_owned.data(), my_n, MPI_LMP_TAGINT,
                   all_ids.data(), counts.data(), displs.data(),
                   MPI_LMP_TAGINT, world);

    // Build hash map global_id -> owner_rank from the gathered data.
    std::unordered_map<tagint, int> id_to_owner;
    id_to_owner.reserve((size_t) total);
    int p = 0;
    for (int i = 0; i < total; ++i) {
      while (p + 1 < comm->nprocs && i >= displs[p + 1]) ++p;
      id_to_owner[all_ids[i]] = p;
    }

    // Fill in cluster_owner_rank for our peer clusters.
    for (int c = 0; c < n_clusters; ++c) {
      if (cluster_owner_rank[c] == -1) {
        auto it = id_to_owner.find(cluster_global_id[c]);
        if (it != id_to_owner.end()) cluster_owner_rank[c] = it->second;
      }
    }
  }

  // cluster_tag is only used inside this routine (the per-iteration
  // forward_comm(comm_mode == 5) reads/writes it).  Release the
  // nmax-sized buffer between reneighbors -- it gets grown back at the
  // top of the next call.
  memory->destroy(cluster_tag);
  cluster_tag = nullptr;
  cluster_tag_alloc = 0;
}

/* ----------------------------------------------------------------------
   Build the per-cluster comm graph: for each OWNED cluster, the list
   of peer ranks that hold atoms locally in this cluster (stored in
   owner_peers_off / owner_peer_rank), and the augmented constraint
   list contributed by peers (owned_aug_constr_*).

   Built via point-to-point exchange (MPI_Alltoallv): each peer rank
   sends to each owner rank a packed buffer of (global_id, n_tags,
   tags..., n_constr, (ta, tb, type) triples) entries -- one per peer
   cluster.  Owner unpacks, looks up the local cluster index via the
   global_id, records the peer rank, and dedups the constraint triples
   into the augmented pool.  Peer-reported atom tags are skipped over
   (the owner no longer stores them).

   Called from build_constraint_list after identify_cluster_owners.
------------------------------------------------------------------------- */

void FixIlves::build_comm_graph()
{
  // Free any prior allocations -- comm graph is rebuilt fresh at each
  // reneighbor and the counts can change arbitrarily.
  free_comm_graph();
  if (n_clusters == 0) return;

  const int me = comm->me;
  const int nlocal_now = atom->nlocal;
  tagint *tag = atom->tag;

  // ---------- Per-cluster local-atom CSR ----------
  // Built only to drive the PEER->OWNER send pack below; not retained.
  // Each constraint contributes c_atom1 and c_atom2; we route through
  // atom->map to canonicalize PBC self-images and skip pure ghosts.

  std::vector<int> mark(atom->nmax, -1);  // -1 = not yet seen in any cluster
  std::vector<int> cluster_lat_off(n_clusters + 1, 0);

  // First pass: count unique local atoms per cluster.
  for (int k = 0; k < n_constr; ++k) {
    const int c = c_cluster[k];
    for (int side = 0; side < 2; ++side) {
      const int raw = (side == 0) ? c_atom1[k] : c_atom2[k];
      const int local_idx = atom->map(tag[raw]);
      if (local_idx < 0 || local_idx >= nlocal_now) continue;
      if (mark[local_idx] == c) continue;   // already counted for this cluster
      mark[local_idx] = c;
      ++cluster_lat_off[c + 1];
    }
  }

  int total_local = 0;
  for (int c = 0; c < n_clusters; ++c) {
    int n = cluster_lat_off[c + 1];
    cluster_lat_off[c + 1] = total_local + n;
    total_local += n;
  }

  std::vector<int> cluster_lat_idx(total_local > 0 ? total_local : 1);

  // Reset marker; second pass fills cluster_lat_idx.
  std::fill(mark.begin(), mark.end(), -1);
  std::vector<int> cursor(n_clusters, 0);
  for (int c = 0; c < n_clusters; ++c) cursor[c] = cluster_lat_off[c];

  for (int k = 0; k < n_constr; ++k) {
    const int c = c_cluster[k];
    for (int side = 0; side < 2; ++side) {
      const int raw = (side == 0) ? c_atom1[k] : c_atom2[k];
      const int local_idx = atom->map(tag[raw]);
      if (local_idx < 0 || local_idx >= nlocal_now) continue;
      if (mark[local_idx] == c) continue;
      mark[local_idx] = c;
      cluster_lat_idx[cursor[c]++] = local_idx;
    }
  }

  // ---------- owner_peers via MPI_Alltoallv exchange ----------
  // Each PEER cluster contributes a message to its owner with this
  // layout (per cluster, packed back-to-back into one buffer per
  // destination rank):
  //
  //   [gid, n_atoms, tag_0..tag_{n-1},
  //         n_constr, ca_0, cb_0, type_0, ..., ca_{m-1}, cb_{m-1}, type_{m-1}]
  //
  // tag_* are atom tags (peer's local atoms in this cluster).
  // ca_i / cb_i are the constraint endpoints' tags; type_i is the
  // constraint type (positive for bonds, negative for angle virtual
  // constraints -- same encoding as c_type[]).  All values are
  // tagint-sized; we cast constraint types to tagint for transport.
  //
  // The owner uses the atom list to know what peer atoms it must
  // gather positions for each step, and uses the constraint list to
  // augment its OWN visible constraint list so the per-cluster Newton
  // matrix is complete (covers atoms beyond the owner's halo).

  const int nprocs = comm->nprocs;
  std::vector<std::vector<tagint> > send_bufs(nprocs);
  for (int c = 0; c < n_clusters; ++c) {
    if (cluster_owner_rank[c] == me) continue;       // owned cluster, no send
    if (cluster_owner_rank[c] < 0) continue;         // no known owner (shouldn't happen)
    const int dst = cluster_owner_rank[c];
    const tagint gid = cluster_global_id[c];
    const int lo = cluster_lat_off[c];
    const int hi = cluster_lat_off[c + 1];
    const int nlat = hi - lo;
    std::vector<tagint> &buf = send_bufs[dst];
    buf.push_back(gid);
    buf.push_back((tagint) nlat);
    for (int s = lo; s < hi; ++s)
      buf.push_back(tag[cluster_lat_idx[s]]);

    // Constraint defs for this peer cluster: walk cluster_offset[c]..
    // cluster_offset[c+1] in c_perm and emit tag pairs + types.
    const int kbeg = cluster_offset[c];
    const int kend = cluster_offset[c + 1];
    buf.push_back((tagint) (kend - kbeg));
    for (int s = kbeg; s < kend; ++s) {
      const int k = c_perm[s];
      buf.push_back(tag[c_atom1[k]]);
      buf.push_back(tag[c_atom2[k]]);
      buf.push_back((tagint) c_type[k]);
    }
  }

  std::vector<int> send_counts(nprocs, 0), send_displs(nprocs, 0);
  for (int p = 0; p < nprocs; ++p) send_counts[p] = (int) send_bufs[p].size();

  std::vector<int> recv_counts(nprocs, 0), recv_displs(nprocs, 0);
  MPI_Alltoall(send_counts.data(), 1, MPI_INT,
               recv_counts.data(), 1, MPI_INT, world);

  // Prefix-sum displacements and flatten send buffer.
  int send_total = 0, recv_total = 0;
  for (int p = 0; p < nprocs; ++p) {
    send_displs[p] = send_total;
    send_total += send_counts[p];
    recv_displs[p] = recv_total;
    recv_total += recv_counts[p];
  }

  std::vector<tagint> send_flat(send_total > 0 ? send_total : 1);
  for (int p = 0; p < nprocs; ++p) {
    if (send_counts[p] == 0) continue;
    std::copy(send_bufs[p].begin(), send_bufs[p].end(),
              send_flat.begin() + send_displs[p]);
  }
  std::vector<tagint> recv_flat(recv_total > 0 ? recv_total : 1);

  MPI_Alltoallv(send_flat.data(), send_counts.data(), send_displs.data(),
                MPI_LMP_TAGINT,
                recv_flat.data(), recv_counts.data(), recv_displs.data(),
                MPI_LMP_TAGINT, world);

  // ---------- Parse received messages and build owner_peer_* ----------
  // Build global_id -> local-cluster-index map for our OWNED clusters
  // so we can dispatch incoming entries.
  std::unordered_map<tagint, int> owned_gid_to_c;
  owned_gid_to_c.reserve((size_t) n_clusters);
  for (int c = 0; c < n_clusters; ++c)
    if (cluster_owner_rank[c] == me)
      owned_gid_to_c[cluster_global_id[c]] = c;

  // First pass: count peer entries per owned cluster.
  memory->grow(owner_peers_off, n_clusters + 1, "ilves:owner_peers_off");
  owner_peers_alloc_clusters = n_clusters;
  for (int c = 0; c <= n_clusters; ++c) owner_peers_off[c] = 0;

  // We'll buffer the parsed entries first so we can build CSR neatly.
  // Each entry remembers where its constraints are in recv_flat[] so
  // the actual copy is deferred to the second pass.  Peer-reported atom
  // tags (the second wire-format field) are skipped over: we no longer
  // store them on the owner side.
  struct PendingEntry {
    int c;            // owned cluster index
    int peer;         // source rank
    int constr_start; // offset in recv_flat[] for the (ta, tb, type) triples
    int n_constr;
  };
  std::vector<PendingEntry> pending;
  pending.reserve((size_t) recv_total / 7);  // rough estimate

  for (int p = 0; p < nprocs; ++p) {
    int pos = recv_displs[p];
    const int end = pos + recv_counts[p];
    while (pos < end) {
      const tagint gid = recv_flat[pos++];
      const int nlat = (int) recv_flat[pos++];
      pos += nlat;                       // skip peer's atom-tag list
      const int nconstr = (int) recv_flat[pos++];
      const int constr_pos = pos;
      pos += 3 * nconstr;
      auto it = owned_gid_to_c.find(gid);
      if (it != owned_gid_to_c.end()) {
        PendingEntry e;
        e.c = it->second;
        e.peer = p;
        e.constr_start = constr_pos;
        e.n_constr = nconstr;
        pending.push_back(e);
        ++owner_peers_off[e.c + 1];
      }
    }
  }

  // Prefix-sum owner_peers_off.
  int total_peer_entries = 0;
  for (int c = 0; c < n_clusters; ++c) {
    int n = owner_peers_off[c + 1];
    owner_peers_off[c + 1] = total_peer_entries + n;
    total_peer_entries += n;
  }
  owner_peers_off[0] = 0;
  n_owner_peer_entries = total_peer_entries;

  memory->grow(owner_peer_rank,
               total_peer_entries > 0 ? total_peer_entries : 1,
               "ilves:owner_peer_rank");
  owner_peer_alloc_entries = total_peer_entries;

  // Second pass: record peer rank per entry.
  std::vector<int> ec_cursor(n_clusters, 0);
  for (int c = 0; c < n_clusters; ++c) ec_cursor[c] = owner_peers_off[c];

  for (const auto &e : pending) {
    const int slot = ec_cursor[e.c]++;
    owner_peer_rank[slot] = e.peer;
  }

  // ---------- Build owned_aug_constr_* ----------
  // Per OWNED cluster, the augmented constraint list contains peer-
  // reported constraints (by tag pair) that the owner does not already
  // have in c_atom1/c_atom2 of its OWN constraints in the same cluster.

  memory->grow(owned_aug_constr_off, n_clusters + 1, "ilves:owned_aug_constr_off");
  for (int c = 0; c <= n_clusters; ++c) owned_aug_constr_off[c] = 0;

  // First pass: count augmented constraints per cluster (dedup against
  // owner's already-present (ta, tb) pairs in the same cluster).

  // Group pending entries by cluster (already CSR-ordered).
  std::vector<std::vector<int> > pending_by_c(n_clusters);
  for (int i = 0; i < (int) pending.size(); ++i)
    pending_by_c[pending[i].c].push_back(i);

  auto pair_key = [](tagint a, tagint b) -> uint64_t {
    uint64_t lo = (uint64_t)((a < b) ? a : b);
    uint64_t hi = (uint64_t)((a < b) ? b : a);
    return (lo << 32) ^ hi;
  };

  int aug_constr_total = 0;
  for (int c = 0; c < n_clusters; ++c) {
    if (cluster_owner_rank[c] != me) continue;
    if (pending_by_c[c].empty()) continue;

    std::unordered_set<uint64_t> seen_constr;
    const int kbeg = cluster_offset[c];
    const int kend = cluster_offset[c + 1];
    for (int s = kbeg; s < kend; ++s) {
      const int k = c_perm[s];
      seen_constr.insert(pair_key(tag[c_atom1[k]], tag[c_atom2[k]]));
    }

    int ncn = 0;
    for (int idx : pending_by_c[c]) {
      const auto &e = pending[idx];
      for (int t = 0; t < e.n_constr; ++t) {
        const tagint ta = recv_flat[e.constr_start + 3*t + 0];
        const tagint tb = recv_flat[e.constr_start + 3*t + 1];
        const uint64_t key = pair_key(ta, tb);
        if (seen_constr.count(key)) continue;
        seen_constr.insert(key);
        ++ncn;
      }
    }
    owned_aug_constr_off[c + 1] = ncn;
    aug_constr_total += ncn;
  }

  // Prefix-sum.
  for (int c = 0; c < n_clusters; ++c)
    owned_aug_constr_off[c + 1] += owned_aug_constr_off[c];

  memory->grow(owned_aug_constr_ta,
               aug_constr_total > 0 ? aug_constr_total : 1,
               "ilves:owned_aug_constr_ta");
  memory->grow(owned_aug_constr_tb,
               aug_constr_total > 0 ? aug_constr_total : 1,
               "ilves:owned_aug_constr_tb");
  memory->grow(owned_aug_constr_type,
               aug_constr_total > 0 ? aug_constr_total : 1,
               "ilves:owned_aug_constr_type");
  owned_aug_constr_alloc = aug_constr_total;

  // Second pass: fill the augmented constraint pool.
  std::vector<int> aug_constr_cursor(n_clusters, 0);
  for (int c = 0; c < n_clusters; ++c)
    aug_constr_cursor[c] = owned_aug_constr_off[c];

  for (int c = 0; c < n_clusters; ++c) {
    if (cluster_owner_rank[c] != me) continue;
    if (pending_by_c[c].empty()) continue;

    std::unordered_set<uint64_t> seen_constr;
    const int kbeg = cluster_offset[c];
    const int kend = cluster_offset[c + 1];
    for (int s = kbeg; s < kend; ++s) {
      const int k = c_perm[s];
      seen_constr.insert(pair_key(tag[c_atom1[k]], tag[c_atom2[k]]));
    }

    for (int idx : pending_by_c[c]) {
      const auto &e = pending[idx];
      for (int t = 0; t < e.n_constr; ++t) {
        const tagint ta = recv_flat[e.constr_start + 3*t + 0];
        const tagint tb = recv_flat[e.constr_start + 3*t + 1];
        const int   ty = (int) recv_flat[e.constr_start + 3*t + 2];
        const uint64_t key = pair_key(ta, tb);
        if (seen_constr.count(key)) continue;
        seen_constr.insert(key);
        const int slot = aug_constr_cursor[c]++;
        owned_aug_constr_ta[slot] = ta;
        owned_aug_constr_tb[slot] = tb;
        owned_aug_constr_type[slot] = ty;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Owner-to-peers broadcast of full cluster constraint set.

   For every OWNED cluster c (cluster_owner_rank[c] == comm->me) we
   pack the FULL constraint set (owner's local constraints from
   cluster_offset[c]..cluster_offset[c+1] PLUS the augmented entries
   collected in build_comm_graph from peer reports) and ship it to
   each peer rank for that cluster.  Each peer parses the messages,
   looks up its corresponding local cluster by global ID, and ADDS any
   constraint it doesn't already see (matched by tag pair).  After
   the exchange, every participating rank's view of every shared
   cluster has the same constraint set -- assuming both endpoint
   atoms are visible in the peer's halo.  Newton then produces
   identical c_lambda on every rank, and apply_constraint_forces
   conserves momentum across ranks.

   Constraints involving atoms NOT in the peer's halo cannot be added
   (no local index to attach to); they remain a small Schwarz residual
   but are typically rare given the augmented-count evidence.

   Returns the number of constraints added on this rank.  When
   non-zero, build_constraint_list re-runs group_by_cluster() so
   cluster_offset / c_perm reflect the augmented list.
------------------------------------------------------------------------- */

int FixIlves::broadcast_full_clusters_to_peers()
{
  const int me = comm->me;
  const int nprocs = comm->nprocs;
  if (nprocs <= 1) return 0;
  if (n_clusters == 0) return 0;

  tagint *tag = atom->tag;

  // ---------- Owner packs full cluster constraint sets ----------
  // Per peer, one message per owned cluster:
  //   [gid, n_constr, ta0, tb0, type0, ta1, tb1, type1, ...]
  std::vector<std::vector<tagint> > send_bufs(nprocs);
  for (int c = 0; c < n_clusters; ++c) {
    if (cluster_owner_rank[c] != me) continue;
    const int n_peers_c = owner_peers_off[c + 1] - owner_peers_off[c];
    if (n_peers_c == 0) continue;

    const tagint gid = cluster_global_id[c];
    const int own_n  = cluster_offset[c + 1] - cluster_offset[c];
    const int aug_n  = owned_aug_constr_off[c + 1] - owned_aug_constr_off[c];
    const int total_n = own_n + aug_n;

    // Pack once, then copy into each peer's buffer.
    std::vector<tagint> msg;
    msg.reserve(2 + 3 * total_n);
    msg.push_back(gid);
    msg.push_back((tagint) total_n);
    for (int s = cluster_offset[c]; s < cluster_offset[c + 1]; ++s) {
      const int k = c_perm[s];
      msg.push_back(tag[c_atom1[k]]);
      msg.push_back(tag[c_atom2[k]]);
      msg.push_back((tagint) c_type[k]);
    }
    for (int t = owned_aug_constr_off[c]; t < owned_aug_constr_off[c + 1]; ++t) {
      msg.push_back(owned_aug_constr_ta[t]);
      msg.push_back(owned_aug_constr_tb[t]);
      msg.push_back((tagint) owned_aug_constr_type[t]);
    }

    for (int p = owner_peers_off[c]; p < owner_peers_off[c + 1]; ++p) {
      const int dst = owner_peer_rank[p];
      std::vector<tagint> &dst_buf = send_bufs[dst];
      dst_buf.insert(dst_buf.end(), msg.begin(), msg.end());
    }
  }

  // ---------- MPI_Alltoallv ----------
  std::vector<int> send_counts(nprocs, 0), send_displs(nprocs, 0);
  for (int p = 0; p < nprocs; ++p) send_counts[p] = (int) send_bufs[p].size();
  std::vector<int> recv_counts(nprocs, 0), recv_displs(nprocs, 0);
  MPI_Alltoall(send_counts.data(), 1, MPI_INT,
               recv_counts.data(), 1, MPI_INT, world);

  int send_total = 0, recv_total = 0;
  for (int p = 0; p < nprocs; ++p) {
    send_displs[p] = send_total;
    send_total += send_counts[p];
    recv_displs[p] = recv_total;
    recv_total += recv_counts[p];
  }

  std::vector<tagint> send_flat(send_total > 0 ? send_total : 1);
  for (int p = 0; p < nprocs; ++p) {
    if (send_counts[p] == 0) continue;
    std::copy(send_bufs[p].begin(), send_bufs[p].end(),
              send_flat.begin() + send_displs[p]);
  }
  std::vector<tagint> recv_flat(recv_total > 0 ? recv_total : 1);

  MPI_Alltoallv(send_flat.data(), send_counts.data(), send_displs.data(),
                MPI_LMP_TAGINT,
                recv_flat.data(), recv_counts.data(), recv_displs.data(),
                MPI_LMP_TAGINT, world);

  // ---------- Peer parses + augments local constraint list ----------
  // Build peer's gid -> local cluster c map for PEER clusters.
  std::unordered_map<tagint, int> peer_gid_to_c;
  peer_gid_to_c.reserve((size_t) n_clusters);
  for (int c = 0; c < n_clusters; ++c)
    if (cluster_owner_rank[c] != me)
      peer_gid_to_c[cluster_global_id[c]] = c;

  auto pair_key = [](tagint a, tagint b) -> uint64_t {
    uint64_t lo = (uint64_t)((a < b) ? a : b);
    uint64_t hi = (uint64_t)((a < b) ? b : a);
    return (lo << 32) ^ hi;
  };

  // Build a SINGLE global seen set across ALL constraints already on
  // this rank (any cluster, owned or peer).  Per-cluster dedup is
  // unsafe because peer can have multiple disjoint local clusters
  // that all belong to the same global cluster -- a constraint
  // present in peer-local-cluster X must not be added again to
  // peer-local-cluster Y when both map to the same owner-side gid.
  std::unordered_set<uint64_t> seen_global;
  seen_global.reserve((size_t) n_constr);
  for (int k = 0; k < n_constr; ++k) {
    seen_global.insert(pair_key(tag[c_atom1[k]], tag[c_atom2[k]]));
  }

  int n_added = 0;
  int n_skipped = 0;
  for (int p = 0; p < nprocs; ++p) {
    int pos = recv_displs[p];
    const int end = pos + recv_counts[p];
    while (pos < end) {
      const tagint gid = recv_flat[pos++];
      const int nc = (int) recv_flat[pos++];
      auto it = peer_gid_to_c.find(gid);
      if (it == peer_gid_to_c.end()) {
        pos += 3 * nc;  // skip; we don't have this cluster
        continue;
      }
      const int c = it->second;

      for (int t = 0; t < nc; ++t) {
        const tagint ta = recv_flat[pos++];
        const tagint tb = recv_flat[pos++];
        const int   ty = (int) recv_flat[pos++];
        const uint64_t key = pair_key(ta, tb);
        if (seen_global.count(key)) continue;

        // Atoms must be in our halo to add the constraint.  If they
        // aren't, this rank's view of the cluster will be incomplete
        // and the cross-rank c_lambda agreement will break.  Count
        // these and abort below if any occurred globally.
        const int ia = atom->map(ta);
        const int ib = atom->map(tb);
        if (ia < 0 || ib < 0) { ++n_skipped; continue; }

        // Build canonical (lower-tag, higher-tag) entry like the bond
        // loop in build_constraint_list, then PBC-image c_atom2.
        int a_idx, b_idx;
        if (ta < tb) { a_idx = ia; b_idx = ib; }
        else         { a_idx = ib; b_idx = ia; }
        b_idx = domain->closest_image(a_idx, b_idx);

        const double dist = (ty > 0) ? bond_distance[ty] : angle_distance[-ty];

        add_constraint(a_idx, b_idx, ty, dist);
        c_cluster[n_constr - 1] = c;  // join peer's existing cluster
        seen_global.insert(key);
        ++n_added;
      }
    }
  }

  // Cluster atoms that exceed the peer's ghost-shell halo would
  // silently leave the peer with an incomplete cluster matrix and
  // re-introduce the cross-rank c_lambda divergence the broadcast is
  // there to prevent.  Fail loudly so the user can widen the comm
  // cutoff (`comm_modify cutoff <value>`) or refactor the topology.
  int n_skipped_all = n_skipped;
  if (nprocs > 1)
    MPI_Allreduce(&n_skipped, &n_skipped_all, 1, MPI_INT, MPI_SUM, world);
  if (n_skipped_all > 0)
    error->all(FLERR, Error::NOLASTLINE,
               "Fix ilves: {} cluster constraints contain atoms outside the communication cutoff. "
               "Increase the cutoff or remove constraints so clusters become smaller.",
               n_skipped_all);

  return n_added;
}

/* ----------------------------------------------------------------------
   Free comm-graph allocations.  Called at the start of every
   build_comm_graph() and at destruction.
------------------------------------------------------------------------- */

void FixIlves::free_comm_graph()
{
  // The arrays are reallocated on next use via memory->grow, so just
  // reset the count metadata to indicate "no live data" without
  // actually destroying.  The destructor handles final teardown.
  n_owner_peer_entries = 0;
}

/* ----------------------------------------------------------------------
   (Re)allocate ilv_bond_* per-fix bond storage to handle nmax atoms with
   ilv_bond_per_atom columns each.  Called from refresh_ilv_bond_data.
------------------------------------------------------------------------- */

void FixIlves::grow_ilv_bond(int nmax)
{
  if (nmax <= ilv_nmax_alloc) return;
  const bigint sz = (bigint) nmax * (bigint) ilv_bond_per_atom;
  memory->grow(ilv_num_bond,  nmax, "ilves:ilv_num_bond");
  memory->grow(ilv_bond_atom, sz,   "ilves:ilv_bond_atom");
  memory->grow(ilv_bond_type, sz,   "ilves:ilv_bond_type");
  ilv_nmax_alloc = nmax;
}

/* ----------------------------------------------------------------------
   Refresh ilv_bond_* from atom->bond_* for local atoms, then forward_comm
   to populate the ghost slots on neighbor ranks.  Called from
   build_constraint_list at every reneighbor.  After this call, both local
   and ghost atoms in [0, nmax) have their bond storage available via the
   ilv_* arrays -- enabling Schwarz overlap in the constraint-list build.

   Note: under newton on bond, atom->bond_* has each bond stored only at
   the lower-tag endpoint.  Copying that as-is preserves the storage
   convention; the constraint-list dedup (`if (ti > tj) continue;`)
   then sees each unique bond exactly once globally regardless of how
   many ranks have one of its endpoints in their halo.
------------------------------------------------------------------------- */

void FixIlves::refresh_ilv_bond_data()
{
  const int nmax = atom->nmax;
  const int nlocal_now = atom->nlocal;
  const int bpa = ilv_bond_per_atom;
  if (bpa <= 0) return;            // no bonds in this atom_style

  grow_ilv_bond(nmax);

  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;

  for (int i = 0; i < nlocal_now; ++i) {
    const int nb = num_bond[i];
    ilv_num_bond[i] = nb;
    const bigint base = (bigint) i * (bigint) bpa;
    for (int b = 0; b < nb && b < bpa; ++b) {
      ilv_bond_atom[base + b] = bond_atom[i][b];
      // store the absolute type (positive); fix_ilves may have negated
      // atom->bond_type[][] via negate_constrained_topology, but for the
      // type lookup we want the underlying bond type, always positive.
      int bt = bond_type[i][b];
      ilv_bond_type[base + b] = (bt < 0) ? -bt : bt;
    }
    // pad unused slots with sentinels so ghost reads of those slots are safe
    for (int b = nb; b < bpa; ++b) {
      ilv_bond_atom[base + b] = 0;
      ilv_bond_type[base + b] = 0;
    }
  }
  // ghost slots will be filled by forward_comm below; initialize to 0
  for (int i = nlocal_now; i < nmax; ++i) {
    ilv_num_bond[i] = 0;
    const bigint base = (bigint) i * (bigint) bpa;
    for (int b = 0; b < bpa; ++b) {
      ilv_bond_atom[base + b] = 0;
      ilv_bond_type[base + b] = 0;
    }
  }

  // forward_comm to populate ghosts
  comm_mode = 2;
  comm->forward_comm(this);
  comm_mode = 0;
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
  // Do not early-return on n_constr == 0: this routine calls forward_comm
  // and solve_constraints (which contain collective MPI calls).  Empty
  // ranks must still participate or the others deadlock.

  // save current f and v, then zero them so unconstrained_update reduces
  // to xshake = x (the current positions).  Flat 3*nlocal-element
  // scratch buffers; correct_coordinates runs at most a few times per
  // run (setup only), so a local std::vector beats holding member
  // buffers of nmax extent between calls.
  std::vector<double> f_save(3 * nlocal);
  std::vector<double> v_save(3 * nlocal);
  for (int i = 0; i < nlocal; ++i) {
    f_save[3*i+0] = f[i][0]; f_save[3*i+1] = f[i][1]; f_save[3*i+2] = f[i][2];
    v_save[3*i+0] = v[i][0]; v_save[3*i+1] = v[i][1]; v_save[3*i+2] = v[i][2];
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
    f[i][0] = f_save[3*i+0]; f[i][1] = f_save[3*i+1]; f[i][2] = f_save[3*i+2];
    v[i][0] = v_save[3*i+0]; v[i][1] = v_save[3*i+1]; v[i][2] = v_save[3*i+2];
  }

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
  // NOTE: do NOT early-return on n_constr == 0.  Under the distributed-
  // topology model, some ranks may have zero local constraints while
  // others have many.  The forward_comm(v) and MPI_Allreduce inside the
  // iteration loop below are collective; skipping the loop on empty
  // ranks would deadlock the others (same as solve_constraints).

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
  // the bond-length constraint).  Ghost velocities are kept fresh between
  // outer iterations by forward_comm(v) below.
  auto velocity_rhs = [&](int k) -> double {
    int a = c_atom1[k];
    int b = c_atom2[k];
    double vxd = v[a][0] - v[b][0];
    double vyd = v[a][1] - v[b][1];
    double vzd = v[a][2] - v[b][2];
    return vxd*c_rx[k] + vyd*c_ry[k] + vzd*c_rz[k];
  };

  // Schwarz iteration over the velocity projection: each rank solves its
  // local per-cluster slice exactly, applies the update to its local
  // atoms, then forward_comms v so partner ranks see the updates before
  // re-evaluating their own rhs.  For clusters that fit on a single
  // rank, the first iteration drives the residual to zero exactly and
  // we exit.  For clusters that span ranks (e.g. all-bonds backbone),
  // multiple iterations converge across the rank boundary.
  //
  // Tolerance is the relative bond-rate residual |gk|/c_rsq[k] (units
  // 1/time).  Reuse the position-solve tolerance: the velocity-projection
  // matrix is identical in structure to the position-solve matrix and
  // converges in the same iter count.
  const double tol_g = tolerance;

  for (int iter = 0; iter < max_iter; ++iter) {
    double max_relres = 0.0;

    for (int c = 0; c < n_clusters; ++c) {
      const int beg = cluster_offset[c];
      const int end = cluster_offset[c + 1];
      const int n_c = end - beg;

      // assemble symmetric A and rhs g
      for (int i = 0; i < n_c * n_c; ++i) lu_A[i] = 0.0;

      for (int s = 0; s < n_c; ++s) {
        int k = c_perm[beg + s];
        lu_A[s*n_c + s] = diag_factor(k) * c_rsq[k];
        double gk = velocity_rhs(k);
        lu_b[s] = gk;
        double relres = fabs(gk) / c_rsq[k];
        if (relres > max_relres) max_relres = relres;
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
      // atom.  Net effect for atom_p: dv = -w_p * mu_k / m_p * r_k.  For
      // 2-atom this reduces to v[a] -= mu * (1/m_a) * r_k and v[b] +=
      // mu * (1/m_b) * r_k as in the original code.
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

    // forward_comm v so partner ranks (and PBC ghosts) see our updates
    // before the next rhs evaluation.
    comm_mode = 1;
    comm->forward_comm(this);
    comm_mode = 0;

    // synchronize convergence across ranks (same pattern as solve_constraints)
    double global_relres = max_relres;
    if (comm->nprocs > 1)
      MPI_Allreduce(&max_relres, &global_relres, 1, MPI_DOUBLE, MPI_MAX, world);
    if (global_relres < tol_g) break;
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

  refresh_constraint_geometry();

  // forward-communicate atom->v so ghost velocities reflect the values
  // just-set on the owning rank's final_integrate.  Needed at every
  // rank count -- not only for MPI ghosts but also for PBC self-images
  // at np == 1.  Some constraints have c_atom2 pointing at a PBC ghost
  // (closest_image returned a ghost index); without this, correct_velocities
  // reads a stale v[ghost] from the last forward_comm of v (typically none
  // since LAMMPS routinely forward_comms positions, not velocities) and
  // the projection injects energy step by step.
  comm_mode = 1;
  comm->forward_comm(this);
  comm_mode = 0;

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
  bytes += (double) atom->nmax * 3 * sizeof(double);   // xshake
  if (store_flag) bytes += (double) maxstore * 3 * sizeof(double);
  bytes += (double) max_constr * (5*sizeof(int) + 3*sizeof(double));
  // Per-angle-type metadata: angle_btype1, angle_btype2, angle_linear, etc.
  // already accounted for in their own allocations (small, O(nangletypes)).
  // No replicated global topology is held -- the constraint list is built
  // from this rank's local atom storage only.
  // Cholesky factor cache (only allocated for the "fast" variant).
  bytes += (double) chol_pool_alloc * sizeof(double);
  return bytes;
}

void FixIlves::grow_arrays(int nmax)
{
  memory->destroy(xshake);
  memory->create(xshake, nmax, 3, "ilves:xshake");
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
  // comm_mode == 2: pack ilv_bond_* (num_bond + bpa partner tags + bpa types)
  if (comm_mode == 2) {
    const int bpa = ilv_bond_per_atom;
    for (int i = 0; i < n; ++i) {
      int j = list[i];
      buf[m++] = (double) ilv_num_bond[j];
      const bigint base = (bigint) j * (bigint) bpa;
      for (int b = 0; b < bpa; ++b) buf[m++] = (double) ilv_bond_atom[base + b];
      for (int b = 0; b < bpa; ++b) buf[m++] = (double) ilv_bond_type[base + b];
    }
    return m;
  }
  // comm_mode == 5: pack cluster_tag (one tagint per atom, encoded as double).
  if (comm_mode == 5) {
    for (int i = 0; i < n; ++i) {
      int j = list[i];
      buf[m++] = (double) cluster_tag[j];
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
  if (comm_mode == 2) {
    const int bpa = ilv_bond_per_atom;
    for (int i = first; i < last; ++i) {
      ilv_num_bond[i] = (int) buf[m++];
      const bigint base = (bigint) i * (bigint) bpa;
      for (int b = 0; b < bpa; ++b) ilv_bond_atom[base + b] = (tagint) buf[m++];
      for (int b = 0; b < bpa; ++b) ilv_bond_type[base + b] = (int) buf[m++];
    }
    return;
  }
  if (comm_mode == 5) {
    for (int i = first; i < last; ++i) {
      cluster_tag[i] = (tagint) buf[m++];
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
   reverse_comm: pack ghost lang_fbuf entries (the per-atom stiff angle
   force contribution from this rank's middle-owners writing into ghost
   slots).  comm_mode == 3 packs lang_fbuf (3 doubles per atom);
   comm_mode == 6 packs lang_vbuf (6 doubles per atom, per-atom virial
   tensor) for the stress/atom path.  unpack adds the received
   contributions into the owner rank's local buffer entries.
------------------------------------------------------------------------- */

int FixIlves::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  if (comm_mode == 6) {
    for (int i = first; i < last; ++i) {
      buf[m++] = lang_vbuf[6*i+0];
      buf[m++] = lang_vbuf[6*i+1];
      buf[m++] = lang_vbuf[6*i+2];
      buf[m++] = lang_vbuf[6*i+3];
      buf[m++] = lang_vbuf[6*i+4];
      buf[m++] = lang_vbuf[6*i+5];
    }
    return m;
  }
  for (int i = first; i < last; ++i) {
    buf[m++] = lang_fbuf[3*i+0];
    buf[m++] = lang_fbuf[3*i+1];
    buf[m++] = lang_fbuf[3*i+2];
  }
  return m;
}

void FixIlves::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m = 0;
  if (comm_mode == 6) {
    for (int i = 0; i < n; ++i) {
      const int j = list[i];
      lang_vbuf[6*j+0] += buf[m++];
      lang_vbuf[6*j+1] += buf[m++];
      lang_vbuf[6*j+2] += buf[m++];
      lang_vbuf[6*j+3] += buf[m++];
      lang_vbuf[6*j+4] += buf[m++];
      lang_vbuf[6*j+5] += buf[m++];
    }
    return;
  }
  for (int i = 0; i < n; ++i) {
    const int j = list[i];
    lang_fbuf[3*j+0] += buf[m++];
    lang_fbuf[3*j+1] += buf[m++];
    lang_fbuf[3*j+2] += buf[m++];
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
   compute_scalar: total restraint energy.  Only nonzero during
   minimization, where min_post_force accumulates the harmonic
   restraint sum sum_k 0.5*kbond*(|r_ab|-d_k)^2 into ebond.
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

  // Newton iter stats: combine across ranks BEFORE the rank-0 print
  // (Allreduce is collective; must be called on every rank).
  bigint newton_nsum = newton_iter_sum;
  bigint newton_nmax = newton_iter_max;
  if (newton_solve_count > 0 && comm->nprocs > 1) {
    MPI_Allreduce(MPI_IN_PLACE, &newton_nsum, 1, MPI_LMP_BIGINT, MPI_SUM, world);
    MPI_Allreduce(MPI_IN_PLACE, &newton_nmax, 1, MPI_LMP_BIGINT, MPI_MAX, world);
  }

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
    // Newton iteration counts: averaged across calls in this stats window,
    // plus the max single-call iter count.  MPI_Allreduce to combine across
    // ranks: sum (then divide by solve-count which is identical on all
    // ranks), and max.
    if (newton_solve_count > 0) {
      // newton_nsum across ranks; each rank does the same number of solves,
      // so divide by (nprocs * solve_count) for the per-call average.
      const double n_avg = (double) newton_nsum /
                           (double) (comm->nprocs * (double) newton_solve_count);
      mesg += fmt::format("Newton iter: avg {:.1f}, max {} (over {} solves x {} ranks)\n",
                          n_avg, newton_nmax, newton_solve_count, comm->nprocs);
    }
    utils::logmesg(lmp, mesg);
  }

  // Reset stats accumulators for the next window.
  newton_iter_sum = 0;
  newton_iter_max = 0;
  newton_solve_count = 0;

  bigint nt = update->ntimestep;
  next_output = nt + output_every;
  if (nt % output_every != 0)
    next_output = (nt / output_every) * output_every + output_every;
}

/* ----------------------------------------------------------------------
   init_topology: under the distributed-topology model we do not gather
   a global bond/angle table.  Each rank uses its own local atom storage
   (atom->bond_*, atom->angle_*); cross-rank constraint atoms are reached
   via the LAMMPS ghost shell.  The only init-time work here is
   computing angle_distance[at] for each constrained angle type from a
   per-rank scan of local angles, with MPI_Allreduce(MAX) to agree on
   the canonical (b1, b2) flanking bond types per angle type.

   Newton on or off bond is supported.  Under newton off both endpoints
   store every bond, so any rank with a local middle atom resolves the
   flanking types directly.  Under newton on a bond is stored only at
   the lower-tag endpoint -- if that endpoint is remote on this rank,
   the middle's flanking-bond lookup returns 0 and the angle_btype1/2
   cache (filled here from whichever rank does have the storage) is
   used as the fallback in build_constraint_list.
------------------------------------------------------------------------- */

void FixIlves::init_topology()
{
  if (!has_angle) {
    for (int at = 1; at <= atom->nangletypes; ++at) {
      angle_distance[at] = 0.0;
      angle_r1[at] = angle_r2[at] = 0.0;
      angle_btype1[at] = angle_btype2[at] = 0;
    }
    record_topology_baseline();
    return;
  }

  const int natypes  = atom->nangletypes;
  const int nlocal_now = atom->nlocal;
  tagint *tag = atom->tag;
  int *mask = atom->mask;
  int *num_bond     = atom->num_bond;
  int **nb_type     = atom->bond_type;
  tagint **nb_atom  = atom->bond_atom;
  int *num_angle    = atom->num_angle;
  int **na_type     = atom->angle_type;
  tagint **na_atom1 = atom->angle_atom1;
  tagint **na_atom2 = atom->angle_atom2;
  tagint **na_atom3 = atom->angle_atom3;

  // Look up the bond type between atom i (local) and a partner with
  // tag tj.  Uses the member lookup_local_bond_type which checks both
  // endpoints' bond storage -- so it works under newton on (where the
  // bond is at the lower-tag endpoint, which may or may not be local
  // here) as long as either endpoint is local on this rank.  Returns 0
  // if neither endpoint is local or the bond does not exist.
  auto lookup_local_btype = [&](int i, tagint tj) -> int {
    if (i < 0 || i >= nlocal_now) return 0;
    return lookup_local_bond_type(tag[i], tj);
  };

  // Per-rank candidate (b1, b2) for each constrained angle type, plus
  // conflict flag.  b1 <= b2 by convention.  Zero means "this rank has
  // not seen a constrained angle of this type yet".
  std::vector<int> local_b1(natypes + 1, 0);
  std::vector<int> local_b2(natypes + 1, 0);
  std::vector<int> local_conflict(natypes + 1, 0);

  for (int i = 0; i < nlocal_now; ++i) {
    if (!(mask[i] & groupbit)) continue;
    const int nang = num_angle[i];
    for (int m = 0; m < nang; ++m) {
      int at = na_type[i][m];
      if (at == 0) continue;
      if (at < 0) at = -at;
      if (at > natypes) continue;
      if (!angle_flag[at]) continue;
      if (na_atom2[i][m] != tag[i]) continue;   // middle-atom owner only

      tagint tA = na_atom1[i][m];
      tagint tC = na_atom3[i][m];

      // Both flanking bonds must be at this local middle atom (newton off).
      int bt1 = lookup_local_btype(i, tA);
      int bt2 = lookup_local_btype(i, tC);
      if (bt1 == 0 || bt2 == 0) continue;

      // Check that the bonds and the angle's endpoints are in our group.
      // (Endpoints may be ghosts; mask is still valid for them.)
      int iA = atom->map(tA);
      int iC = atom->map(tC);
      if (iA < 0 || iC < 0) continue;
      if (!(mask[iA] & groupbit) || !(mask[iC] & groupbit)) continue;

      // Both flanking bond types must be selected (bond_flag, type_flag, mass_list).
      if (!bond_selected_for_atoms(i, iA, bt1)) continue;
      if (!bond_selected_for_atoms(i, iC, bt2)) continue;

      int bmin = MIN(bt1, bt2);
      int bmax = MAX(bt1, bt2);
      if (local_b1[at] == 0) {
        local_b1[at] = bmin;
        local_b2[at] = bmax;
      } else if (local_b1[at] != bmin || local_b2[at] != bmax) {
        local_conflict[at] = 1;
      }
    }
  }

  // Cross-rank consensus.  Each rank may have b1=b2=0 (no local constrained
  // angles of this type) or a nonzero pair.  We Allreduce MAX over the
  // pairs (treating zero as "no info").  Then a second Allreduce MAX over
  // an in-place conflict flag also fires if two ranks disagreed on a
  // nonzero pair.
  std::vector<int> global_b1(natypes + 1, 0);
  std::vector<int> global_b2(natypes + 1, 0);
  std::vector<int> global_conflict(natypes + 1, 0);
  if (comm->nprocs > 1) {
    MPI_Allreduce(local_b1.data(), global_b1.data(), natypes + 1, MPI_INT, MPI_MAX, world);
    MPI_Allreduce(local_b2.data(), global_b2.data(), natypes + 1, MPI_INT, MPI_MAX, world);

    // Now detect a cross-rank disagreement: a rank with a nonzero
    // (b1, b2) different from the agreed-on global pair is a conflict.
    for (int at = 1; at <= natypes; ++at) {
      if (local_b1[at] != 0 &&
          (local_b1[at] != global_b1[at] || local_b2[at] != global_b2[at]))
        local_conflict[at] = 1;
    }
    MPI_Allreduce(local_conflict.data(), global_conflict.data(),
                  natypes + 1, MPI_INT, MPI_MAX, world);
  } else {
    for (int at = 0; at <= natypes; ++at) {
      global_b1[at] = local_b1[at];
      global_b2[at] = local_b2[at];
      global_conflict[at] = local_conflict[at];
    }
  }

  for (int at = 1; at <= natypes; ++at) {
    if (global_conflict[at])
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix ilves: angle type {} has mixed type bonds", at);

    if (!angle_flag[at]) {
      angle_distance[at] = 0.0;
      angle_r1[at] = angle_r2[at] = 0.0;
      angle_btype1[at] = angle_btype2[at] = 0;
      continue;
    }
    if (global_b1[at] == 0) {
      // User requested constraints on this angle type via `a`, but no
      // rank found a constrained angle of that type whose flanking
      // bonds it could resolve locally.  The AC virtual constraint
      // gets silently dropped -- warn so the user notices.
      if (comm->me == 0)
        error->warning(FLERR,
                       "Fix ilves: angle type {} was selected but could not determine its "
                       "bond types; the angle's A-C virtual constraint is disabled", at);
      angle_distance[at] = 0.0;
      angle_r1[at] = angle_r2[at] = 0.0;
      angle_btype1[at] = angle_btype2[at] = 0;
      continue;
    }
    int b1 = global_b1[at];
    int b2 = global_b2[at];
    const double theta0 = force->angle->equilibrium_angle(at);
    const double r1 = bond_distance[b1];
    const double r2 = bond_distance[b2];
    angle_r1[at] = r1;
    angle_r2[at] = r2;
    angle_distance[at] = sqrt(r1*r1 + r2*r2 - 2.0*r1*r2*cos(theta0));
    angle_btype1[at] = b1;
    angle_btype2[at] = b2;
  }

  // Asymmetric angles + partial `b` selector under newton on bond.
  // Under newton on, a bond is stored only at its lower-tag endpoint;
  // when both endpoints of a flanking bond are remote on this rank,
  // lookup_local_bond_type returns 0 and build_constraint_list falls
  // back to accepting the angle if EITHER cached flanking type
  // (b1 or b2) is user-selected.  For asymmetric types (b1 != b2)
  // with only one of bond_flag[b1] / bond_flag[b2] set, that fallback
  // is slightly over-permissive -- it may constrain an angle whose
  // actual flanking type the user did not select.
  if (force->newton_bond && comm->me == 0) {
    for (int at = 1; at <= natypes; ++at) {
      if (!angle_flag[at]) continue;
      const int b1 = angle_btype1[at];
      const int b2 = angle_btype2[at];
      if (b1 == 0 || b2 == 0 || b1 == b2) continue;
      const bool s1 = bond_flag[b1] != 0;
      const bool s2 = bond_flag[b2] != 0;
      if (s1 == s2) continue;
      error->warning(FLERR,
                     "Fix ilves: angle type {} has asymmetric bond types ({} and {}) and "
                     "only one is selected for constraining via 'b'; "
                     "this may lead to over-constraining", at, b1, b2);
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
  baseline_ready = true;
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
   Look up the type of the bond between tags ta and tb by scanning the
   bond storage of either endpoint.  Under newton off both endpoints
   store the bond; under newton on only the lower-tag endpoint does.
   The per-fix ilv_bond_* arrays (refreshed at reneighbor and
   forward_commed to ghosts) carry storage for halo atoms too, so a
   newton-on bond is found whenever the lower-tag endpoint is in this
   rank's halo.  Returns 0 if no such bond is visible locally.
------------------------------------------------------------------------- */

int FixIlves::lookup_local_bond_type(tagint ta, tagint tb)
{
  // Use the per-fix ilv_bond_* arrays (refreshed at every reneighbor by
  // refresh_ilv_bond_data) which include ghost-side bond storage too.
  // Falls back to atom->bond_* if ilv_* hasn't been populated yet
  // (e.g. on the very first build_constraint_list before reneighbor).
  const int bpa = ilv_bond_per_atom;
  if (bpa > 0 && ilv_num_bond) {
    int ia = atom->map(ta);
    if (ia >= 0 && ia < ilv_nmax_alloc) {
      const int nb = ilv_num_bond[ia];
      const bigint base = (bigint) ia * (bigint) bpa;
      for (int b = 0; b < nb && b < bpa; ++b) {
        if (ilv_bond_atom[base + b] == tb) return ilv_bond_type[base + b];
      }
    }
    int ib = atom->map(tb);
    if (ib >= 0 && ib < ilv_nmax_alloc) {
      const int nb = ilv_num_bond[ib];
      const bigint base = (bigint) ib * (bigint) bpa;
      for (int b = 0; b < nb && b < bpa; ++b) {
        if (ilv_bond_atom[base + b] == ta) return ilv_bond_type[base + b];
      }
    }
    return 0;
  }

  // Pre-refresh fallback: scan atom->bond_* directly (local only).
  const int nlocal_now = atom->nlocal;
  int *num_bond     = atom->num_bond;
  int **nb_type     = atom->bond_type;
  tagint **nb_atom  = atom->bond_atom;

  int ia = atom->map(ta);
  if (ia >= 0 && ia < nlocal_now) {
    const int nb = num_bond[ia];
    for (int b = 0; b < nb; ++b) {
      if (nb_atom[ia][b] != tb) continue;
      int bt = nb_type[ia][b];
      if (bt < 0) bt = -bt;
      return bt;
    }
  }
  int ib = atom->map(tb);
  if (ib >= 0 && ib < nlocal_now) {
    const int nb = num_bond[ib];
    for (int b = 0; b < nb; ++b) {
      if (nb_atom[ib][b] != ta) continue;
      int bt = nb_type[ib][b];
      if (bt < 0) bt = -bt;
      return bt;
    }
  }
  return 0;
}

/* ----------------------------------------------------------------------
   True if the bond between tags ta and tb is selected for constraint.
   Uses local bond storage to find the type, then applies the selectors
   (bond_flag, type_flag, mass_list).  Returns false if either atom is
   not locally available (cannot apply atom-type / mass selectors), or
   if no bond exists between them locally.
------------------------------------------------------------------------- */

bool FixIlves::bond_is_constrained(tagint ta, tagint tb)
{
  int bt = lookup_local_bond_type(ta, tb);
  if (bt == 0) return false;
  int ia = atom->map(ta);
  int ib = atom->map(tb);
  if (ia < 0 || ib < 0) return false;
  return bond_selected_for_atoms(ia, ib, bt);
}

/* ----------------------------------------------------------------------
   Build the flat constraint list from this rank's local atom bond/angle
   storage (distributed topology, no global gather).

   Bond pass.  Walk each local atom's bond storage; for every bond
   (i, partner_tag) with tag[i] < partner_tag (dedup under newton off),
   add a constraint if the partner is locally addressable (in the ghost
   shell), both atoms are in the fix group, and selectors match.

   Angle pass.  Walk each local atom's angle storage; for every angle
   where the local atom is the middle (na_atom2 == tag[i]) and the angle
   is not flagged near-linear, add the A-C virtual constraint if both
   endpoints are locally addressable, all three atoms are in the fix
   group, and both flanking bonds are constrained.  This rank's middle
   atom is the unique owner of the angle's AC constraint.

   Cross-rank constraints reach the partner atom through the LAMMPS
   ghost shell.  Constraints between two ghosts -- i.e. constraints
   wholly outside this rank's local atoms -- are not added here; they
   are added by the rank that owns their lower-tag (middle) endpoint
   locally.  Cross-rank cluster pieces converge via the Newton loop's
   forward_comm(xshake) between iterations.
------------------------------------------------------------------------- */

void FixIlves::build_constraint_list()
{
  const int nlocal_now = atom->nlocal;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  int *num_angle    = atom->num_angle;
  int **na_type     = atom->angle_type;
  tagint **na_atom1 = atom->angle_atom1;
  tagint **na_atom2 = atom->angle_atom2;
  tagint **na_atom3 = atom->angle_atom3;

  // Schwarz overlap: refresh per-fix ilv_bond_* (local + ghost) so we
  // can walk the full halo's bond storage in the bond loop below.
  refresh_ilv_bond_data();

  n_constr = 0;

  // ------------------------------------------------------------------
  // Bond pass.  Walk EVERY atom in this rank's halo (local + ghost)
  // using the per-fix ilv_bond_* arrays.  Dedup by `ti < tj` (lower-tag
  // atom adds the constraint); this works because each bond is
  // canonically stored at its lower-tag endpoint regardless of newton
  // flag, plus the refresh has populated ghost-side storage so we can
  // see ghost-ghost bonds visible in the halo.  Iteration via atom
  // indices visits each unique tag once (filter by `atom->map(ti) == i`
  // to skip non-canonical PBC ghost copies of the same tag).
  // ------------------------------------------------------------------
  const int bpa = ilv_bond_per_atom;
  const int halo_end = nlocal_now + atom->nghost;
  for (int i = 0; i < halo_end; ++i) {
    const tagint ti = tag[i];
    if (ti <= 0) continue;
    if (atom->map(ti) != i) continue;       // skip non-canonical PBC copies
    if (!(mask[i] & groupbit)) continue;

    const int nb = (i < ilv_nmax_alloc) ? ilv_num_bond[i] : 0;
    const bigint base = (bigint) i * (bigint) bpa;
    for (int b = 0; b < nb && b < bpa; ++b) {
      int bt = ilv_bond_type[base + b];
      if (bt <= 0 || bt > atom->nbondtypes) continue;
      const tagint tj = ilv_bond_atom[base + b];
      if (ti > tj) continue;                // dedup: lower-tag adds
      int j = atom->map(tj);
      if (j < 0) continue;                  // partner outside halo
      if (!(mask[j] & groupbit)) continue;
      if (!bond_selected_for_atoms(i, j, bt)) continue;

      const int a_idx = i;
      int b_idx = domain->closest_image(a_idx, j);
      add_constraint(a_idx, b_idx, bt, bond_distance[bt]);
    }
  }

  // ------------------------------------------------------------------
  // Angle pass: A-C virtual constraints.  Near-linear angle types
  // (angle_linear[at] == 1) are silently skipped -- the AC entry would
  // be rank-deficient against the two flanking bonds; the standard
  // angle_style force-field term continues to act on those angles
  // because negate_constrained_topology() skips them too.
  // ------------------------------------------------------------------
  // Dedup rule for angles (analogous to the bond rule above): with
  // newton off, the angle is stored at all three atoms.  Each rank
  // processes the angle once, at the local atom with the smallest
  // tag among (t1, t2, t3) that is local on this rank.
  //   - single-rank: all three local; the smallest-tag local processes
  //   - multi-rank with t1, t2, t3 split: each rank's only local atom
  //     of the three trivially wins -- every involved rank processes,
  //     same convergent solve via forward_comm(xshake)
  //   - if only the middle is local: it still processes but applies no
  //     force (forces go to t1/t3 endpoints).  Endpoints' ranks also
  //     process and apply force to their local copy.
  if (has_angle) {
    for (int i = 0; i < nlocal_now; ++i) {
      if (!(mask[i] & groupbit)) continue;
      const tagint ti = tag[i];
      const int nang = num_angle[i];
      for (int m = 0; m < nang; ++m) {
        int at = na_type[i][m];
        if (at == 0) continue;
        if (at < 0) at = -at;
        if (at > atom->nangletypes) continue;
        if (!angle_flag[at]) continue;
        if (angle_linear[at]) continue;                    // angle term handles it

        const tagint t1 = na_atom1[i][m];
        const tagint t2 = na_atom2[i][m];
        const tagint t3 = na_atom3[i][m];
        const int i1 = atom->map(t1);
        const int i2 = atom->map(t2);
        const int i3 = atom->map(t3);

        // Smallest-tag local-on-this-rank dedup: only the smallest-tag
        // among (t1, t2, t3) that is local on this rank processes.
        tagint smallest_local = -1;
        if (i1 >= 0 && i1 < nlocal_now)
          if (smallest_local < 0 || t1 < smallest_local) smallest_local = t1;
        if (i2 >= 0 && i2 < nlocal_now)
          if (smallest_local < 0 || t2 < smallest_local) smallest_local = t2;
        if (i3 >= 0 && i3 < nlocal_now)
          if (smallest_local < 0 || t3 < smallest_local) smallest_local = t3;
        if (ti != smallest_local) continue;

        // The constraint endpoints (t1, t3) and the middle (t2) must
        // all be reachable in this rank's halo for the constraint
        // geometry and verification to work.
        if (i1 < 0 || i2 < 0 || i3 < 0) continue;
        if (!(mask[i1] & groupbit)) continue;
        if (!(mask[i2] & groupbit)) continue;
        if (!(mask[i3] & groupbit)) continue;

        // Verify both flanking bonds are constrained.  First try the
        // local-bond-storage lookup (correct under newton off bond, and
        // also under newton on whenever the bond's lower-tag endpoint is
        // local on this rank).  If that lookup returns 0 (newton on bond
        // with the lower-tag endpoint on a different rank), fall back to
        // the per-angle-type cache built at init by init_topology: bond
        // type is angle_btype1[at] or angle_btype2[at] for any angle of
        // type at.  Apply selectors using the cached type + the angle's
        // local atom data.
        auto flank_constrained = [&](tagint tmid, int imid, tagint te, int ie) -> bool {
          int bt = lookup_local_bond_type(tmid, te);
          if (bt > 0) return bond_selected_for_atoms(imid, ie, bt);
          // Fallback under newton on with both bond endpoints remote.
          // The cache may have two candidates (angle_btype1/2); accept
          // the bond if EITHER candidate's selectors fire.  Exact for
          // symmetric angles (b1 == b2); slightly over-permissive for
          // asymmetric angles where only one of b1, b2 is user-selected.
          int b1 = angle_btype1[at];
          int b2 = angle_btype2[at];
          if (b1 > 0 && bond_selected_for_atoms(imid, ie, b1)) return true;
          if (b2 > 0 && b2 != b1 && bond_selected_for_atoms(imid, ie, b2)) return true;
          return false;
        };
        if (!flank_constrained(t2, i2, t1, i1)) continue;
        if (!flank_constrained(t2, i2, t3, i3)) continue;

        // Pick a_idx/b_idx so the lower-tag endpoint is c_atom1.  Both
        // endpoints are locally addressable (i1, i3); apply PBC
        // closest-image relative to a_idx.
        int a_idx, b_idx;
        if (t1 < t3) { a_idx = i1; b_idx = i3; }
        else         { a_idx = i3; b_idx = i1; }
        b_idx = domain->closest_image(a_idx, b_idx);

        add_constraint(a_idx, b_idx, -at, angle_distance[at]);
      }
    }
  }

  // connected-component labelling and cluster grouping
  group_by_cluster();
  // determine the global owner of each local cluster (rank holding the
  // lowest-tag atom of the cluster).
  identify_cluster_owners();
  // build the per-cluster comm graph and broadcast the full cluster
  // constraint set from owner to peers so every participating rank
  // has the same cluster view.
  build_comm_graph();
  const int n_added = broadcast_full_clusters_to_peers();
  if (n_added > 0) {
    // Constraint list grew on this rank.  Re-run group_by_cluster to
    // rebuild cluster_offset / c_perm for the augmented list, then
    // re-run identify_cluster_owners so cluster_global_id and
    // cluster_owner_rank match the new local cluster indices (used
    // by the setup() owned/remote-count stats).  The owner_peers_*
    // / owned_aug_* arrays remain indexed by the OLD cluster IDs and
    // would be stale here, but they've already been consumed by the
    // broadcast and are rebuilt fresh at the next reneighbor.
    group_by_cluster();
    identify_cluster_owners();
  }
  precompute_constraint_data();
  build_apply_order();
}
