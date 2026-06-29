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
   Contributing author: Axel Kohlmeyer (Temple U)
   Based on the symd reference implementation at https://github.com/whitead/symd
------------------------------------------------------------------------- */

#include "fix_symmetry.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "respa.h"
#include "symmetry_group.h"
#include "update.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <unordered_set>

using namespace LAMMPS_NS;
using namespace FixConst;

namespace {
constexpr char cite_symd[] = "fix symmetry: https://doi.org/10.1021/acs.jctc.2c00401\n\n"
                             "@Article{WhiteCoxSymmetricMD,\n"
                             "author = {Cox, Sam and White, Andrew D.},\n"
                             "title = {Symmetric Molecular Dynamics},\n"
                             "journal = {Journal of Chemical Theory and Computation},\n"
                             "volume = {18},\n"
                             "number = {7},\n"
                             "pages = {4077-4081},\n"
                             "year = {2022},\n"
                             "doi = {10.1021/acs.jctc.2c00401},\n"
                             "note ={PMID: 35699649},\n"
                             "url = {https://doi.org/10.1021/acs.jctc.2c00401},\n"
                             "}\n\n";

/* ----------------------------------------------------------------------
   Cartesian <-> fractional vector conversions. These are the vector
   (no boxlo offset) analogues of Domain::lamda2x / x2lamda. h and h_inv
   are stored in Voigt order: xx, yy, zz, yz, xz, xy.

       v_cart[0] = h[0]*v_frac[0] + h[5]*v_frac[1] + h[4]*v_frac[2]
       v_cart[1] =                  h[1]*v_frac[1] + h[3]*v_frac[2]
       v_cart[2] =                                   h[2]*v_frac[2]
------------------------------------------------------------------------- */

inline void frac_to_cart_vec(double out[3], const double in[3], const double h[6])
{
  out[0] = h[0] * in[0] + h[5] * in[1] + h[4] * in[2];
  out[1] = h[1] * in[1] + h[3] * in[2];
  out[2] = h[2] * in[2];
}

inline void cart_to_frac_vec(double out[3], const double in[3], const double h_inv[6])
{
  out[0] = h_inv[0] * in[0] + h_inv[5] * in[1] + h_inv[4] * in[2];
  out[1] = h_inv[1] * in[1] + h_inv[3] * in[2];
  out[2] = h_inv[2] * in[2];
}

// out = M * in for a 3x3 matrix in row-major storage M[i][j].
inline void apply_3x3(double out[3], const double M[3][3], const double in[3])
{
  out[0] = M[0][0] * in[0] + M[0][1] * in[1] + M[0][2] * in[2];
  out[1] = M[1][0] * in[0] + M[1][1] * in[1] + M[1][2] * in[2];
  out[2] = M[2][0] * in[0] + M[2][1] * in[1] + M[2][2] * in[2];
}
}    // namespace

/* ---------------------------------------------------------------------- */

FixSymmetry::FixSymmetry(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), grp(nullptr), tol(1.0e-6), symfile(nullptr), fsum_local(nullptr),
    fsum_global(nullptr), vsum_local(nullptr), vsum_global(nullptr), xasym_local(nullptr),
    xasym_global(nullptr), n_orbits(0)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix symmetry", error);

  symfile = utils::strdup(arg[3]);

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "tol") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix symmetry tol", error);
      tol = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (tol <= 0.0) error->all(FLERR, iarg + 1, "Fix symmetry tol must be > 0");
      iarg += 2;
    } else {
      error->all(FLERR, iarg, "Unknown fix symmetry keyword: {}", arg[iarg]);
    }
  }

  grp = new SymmetryGroup(lmp);
  grp->read(symfile);

  if (lmp->citeme) lmp->citeme->add(cite_symd);
}

/* ---------------------------------------------------------------------- */

FixSymmetry::~FixSymmetry()
{
  delete grp;
  delete[] symfile;
  memory->destroy(fsum_local);
  memory->destroy(fsum_global);
  memory->destroy(vsum_local);
  memory->destroy(vsum_global);
  memory->destroy(xasym_local);
  memory->destroy(xasym_global);
}

/* ---------------------------------------------------------------------- */

int FixSymmetry::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSymmetry::init()
{
  validate_box();
  build_orbit_map();

  n_orbits = grp->n_orbits();
  memory->grow(fsum_local, 3 * n_orbits, "symmetry:fsum_local");
  memory->grow(fsum_global, 3 * n_orbits, "symmetry:fsum_global");
  memory->grow(vsum_local, 3 * n_orbits, "symmetry:vsum_local");
  memory->grow(vsum_global, 3 * n_orbits, "symmetry:vsum_global");
  memory->grow(xasym_local, 3 * n_orbits, "symmetry:xasym_local");
  memory->grow(xasym_global, 3 * n_orbits, "symmetry:xasym_global");
}

/* ---------------------------------------------------------------------- */

void FixSymmetry::setup(int vflag)
{
  // Per-level force symmetrization for RESPA: each level's f contribution
  // must be individually symmetric so the velocity update from the sum
  // f_outer + f_inner stays symmetric. The pattern mirrors fix_manifoldforce.
  if (utils::strmatch(update->integrate_style, "^verlet")) {
    post_force(vflag);
  } else {
    auto *respa_ptr = dynamic_cast<Respa *>(update->integrate);
    const int nlevels = respa_ptr->nlevels;
    for (int ilevel = 0; ilevel < nlevels; ++ilevel) {
      respa_ptr->copy_flevel_f(ilevel);
      post_force_respa(vflag, ilevel, 0);
      respa_ptr->copy_f_flevel(ilevel);
    }
  }
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixSymmetry::post_force_respa(int vflag, int /*ilevel*/, int /*iloop*/)
{
  // Modify the per-level forces in atom->f at every level. The RESPA
  // integrator's framework copies atom->f to/from fix_respa->f_level
  // around the post_force_respa call, so symmetrizing atom->f here
  // is equivalent to symmetrizing that level's f_level entry. Velocity
  // symmetrization happens repeatedly (idempotent) which is harmless.
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Symmetrize per-orbit forces and velocities.

   Two-phase algorithm, all per-orbit accumulation done in fractional
   coords so the symmetry operators (which act on fractional coords)
   apply directly:

     pass 1: every local atom i in orbit o, op k contributes
             R_k^{-1} * (h_inv * f_cart_i)  to fsum_local[o].
             Same for velocities -> vsum_local[o].
     Allreduce(SUM) across ranks.
     pass 2: every local atom i in orbit o, op k overwrites its force
             with  h * R_k * (fsum_global[o] / n_ops),  same for velocity.

   Atoms not declared in the symmetry file are left untouched.
------------------------------------------------------------------------- */

void FixSymmetry::post_force(int /*vflag*/)
{
  if (n_orbits == 0) return;

  const int n_ops = grp->n_ops();
  const double inv_n = 1.0 / static_cast<double>(n_ops);
  const double *h = domain->h;
  const double *h_inv = domain->h_inv;

  const int nlocal = atom->nlocal;
  double **f = atom->f;
  double **v = atom->v;
  tagint *tag = atom->tag;

  std::fill_n(fsum_local, 3 * n_orbits, 0.0);
  std::fill_n(vsum_local, 3 * n_orbits, 0.0);

  // pass 1: project each local atom's f and v back to the asymmetric unit
  // and accumulate per orbit. For Wyckoff asym atoms (those with a
  // non-empty stabilizer), apply contributions for every op that maps
  // the asym atom onto itself in addition to op 0 -- otherwise the sum
  // would under-count by a factor of (1 + |stabilizer|) and the / n_ops
  // normalization in pass 2 would give the wrong magnitude.
  for (int i = 0; i < nlocal; ++i) {
    auto it = tag_info.find(tag[i]);
    if (it == tag_info.end()) continue;
    const int orbit = it->second.first;
    const int op_id = it->second.second;

    double f_frac[3], v_frac[3];
    cart_to_frac_vec(f_frac, f[i], h_inv);
    cart_to_frac_vec(v_frac, v[i], h_inv);

    auto accumulate = [&](int k) {
      const SymmOp &op = grp->op(k);
      double r[3];
      apply_3x3(r, op.Rinv, f_frac);
      fsum_local[3 * orbit + 0] += r[0];
      fsum_local[3 * orbit + 1] += r[1];
      fsum_local[3 * orbit + 2] += r[2];
      apply_3x3(r, op.Rinv, v_frac);
      vsum_local[3 * orbit + 0] += r[0];
      vsum_local[3 * orbit + 1] += r[1];
      vsum_local[3 * orbit + 2] += r[2];
    };
    accumulate(op_id);
    if (op_id == 0) {
      const Orbit &orb = grp->orbit(orbit);
      for (int k : orb.stabilizer) accumulate(k);
    }
  }

  MPI_Allreduce(fsum_local, fsum_global, 3 * n_orbits, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(vsum_local, vsum_global, 3 * n_orbits, MPI_DOUBLE, MPI_SUM, world);

  // pass 2: redistribute the per-orbit symmetric value to every image
  // (including the asymmetric representative, where op_id == 0 / R == I).
  for (int i = 0; i < nlocal; ++i) {
    auto it = tag_info.find(tag[i]);
    if (it == tag_info.end()) continue;
    const int orbit = it->second.first;
    const int op_id = it->second.second;
    const SymmOp &op = grp->op(op_id);

    double fsym_frac[3] = {fsum_global[3 * orbit + 0] * inv_n, fsum_global[3 * orbit + 1] * inv_n,
                           fsum_global[3 * orbit + 2] * inv_n};
    double vsym_frac[3] = {vsum_global[3 * orbit + 0] * inv_n, vsum_global[3 * orbit + 1] * inv_n,
                           vsum_global[3 * orbit + 2] * inv_n};

    double image_frac[3];

    apply_3x3(image_frac, op.R, fsym_frac);
    frac_to_cart_vec(f[i], image_frac, h);

    apply_3x3(image_frac, op.R, vsym_frac);
    frac_to_cart_vec(v[i], image_frac, h);
  }
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Re-fold image positions from the asymmetric representative. This cancels
   ULP-level drift the integrator accumulates between the asym atom and its
   images even when the symmetrized forces are exact.

   Two-phase:
     pass 1: every rank publishes the fractional position of any locally
             owned asymmetric atom into xasym_local[orbit] (zero elsewhere).
             Allreduce(SUM) collapses to the unique owner -> xasym_global.
     pass 2: every rank rewrites the positions of its local image atoms
             (op_id > 0) to  s_target + round(s_current - s_target),
             i.e. the symmetry image of the asymmetric atom *in the
             periodic copy nearest the image atom's current location*.
             That preserves whatever unit-cell offset the user chose at
             setup and is robust against the integrator transiently moving
             atoms outside [0,1)^3 between comm steps.

   Asymmetric atoms (op_id == 0) are the source of truth and are not
   modified. Atoms not declared in the symmetry file are untouched.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   For each locally owned Wyckoff asym atom, project its position and
   velocity onto the constraint manifold defined by the orbit's stabilizer.
   The projection uses the precomputed pseudo-inverse Binv and the per-orbit
   constant wyckoff_c (built at init from the user's initial placement).

       residual r = B * s + c     (zero when atom is on the Wyckoff site)
       Delta_s    = Binv * r      (smallest fractional shift back to the site)
       s_new      = s - Delta_s
       v_frac_new = (I - Binv * B) * v_frac     (drop velocity components
                    perpendicular to the constraint manifold)

   No MPI: each asym atom lives on exactly one rank, and Binv/c are
   replicated on all ranks. This must run BEFORE the image-refold phase
   of end_of_step so that the published asym positions already satisfy
   the site-symmetry constraint.
------------------------------------------------------------------------- */

void FixSymmetry::wyckoff_project_locals()
{
  const double *h = domain->h;
  const double *h_inv = domain->h_inv;
  const int nlocal = atom->nlocal;
  double **x = atom->x;
  double **v = atom->v;
  tagint *tag = atom->tag;

  for (int i = 0; i < nlocal; ++i) {
    auto it = tag_info.find(tag[i]);
    if (it == tag_info.end()) continue;
    if (it->second.second != 0) continue;    // images are handled by refold
    const int orbit = it->second.first;
    const Orbit &orb = grp->orbit(orbit);
    if (orb.stabilizer.empty()) continue;    // general position: nothing to do

    // position projection: Delta_s = Binv * (B * s + c)
    double s[3];
    domain->x2lamda(x[i], s);
    double res[3] = {wyckoff_c[3 * orbit + 0], wyckoff_c[3 * orbit + 1], wyckoff_c[3 * orbit + 2]};
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) res[a] += orb.B[a][b] * s[b];
    double ds[3] = {0.0, 0.0, 0.0};
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) ds[a] += orb.Binv[a][b] * res[b];
    double s_new[3] = {s[0] - ds[0], s[1] - ds[1], s[2] - ds[2]};
    domain->lamda2x(s_new, x[i]);

    // velocity projection: drop the perpendicular component
    //   v_frac_new = v_frac - Binv * B * v_frac
    double v_frac[3];
    cart_to_frac_vec(v_frac, v[i], h_inv);
    double Bv[3] = {0.0, 0.0, 0.0};
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) Bv[a] += orb.B[a][b] * v_frac[b];
    double BinvBv[3] = {0.0, 0.0, 0.0};
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) BinvBv[a] += orb.Binv[a][b] * Bv[b];
    double v_frac_new[3] = {v_frac[0] - BinvBv[0], v_frac[1] - BinvBv[1], v_frac[2] - BinvBv[2]};
    frac_to_cart_vec(v[i], v_frac_new, h);
  }
}

/* ---------------------------------------------------------------------- */

void FixSymmetry::end_of_step()
{
  if (n_orbits == 0) return;

  // Wyckoff projection runs first so the published asym positions in
  // pass 1 below already satisfy site symmetry.
  wyckoff_project_locals();

  const int nlocal = atom->nlocal;
  double **x = atom->x;
  tagint *tag = atom->tag;

  // pass 1: publish each locally owned asymmetric atom's fractional position.
  std::fill_n(xasym_local, 3 * n_orbits, 0.0);
  for (int i = 0; i < nlocal; ++i) {
    auto it = tag_info.find(tag[i]);
    if (it == tag_info.end()) continue;
    if (it->second.second != 0) continue;    // not the asym representative
    const int orbit = it->second.first;
    double s[3];
    domain->x2lamda(x[i], s);
    xasym_local[3 * orbit + 0] = s[0];
    xasym_local[3 * orbit + 1] = s[1];
    xasym_local[3 * orbit + 2] = s[2];
  }
  MPI_Allreduce(xasym_local, xasym_global, 3 * n_orbits, MPI_DOUBLE, MPI_SUM, world);

  // pass 2: rewrite each local image atom's position.
  for (int i = 0; i < nlocal; ++i) {
    auto it = tag_info.find(tag[i]);
    if (it == tag_info.end()) continue;
    const int op_id = it->second.second;
    if (op_id == 0) continue;
    const int orbit = it->second.first;
    const SymmOp &op = grp->op(op_id);

    const double s_asym[3] = {xasym_global[3 * orbit + 0], xasym_global[3 * orbit + 1],
                              xasym_global[3 * orbit + 2]};

    // s_target = R * s_asym + t
    double s_target[3];
    apply_3x3(s_target, op.R, s_asym);
    s_target[0] += op.t[0];
    s_target[1] += op.t[1];
    s_target[2] += op.t[2];

    // s_current of this image
    double s_current[3];
    domain->x2lamda(x[i], s_current);

    // pick the periodic copy of s_target closest to s_current
    double s_new[3];
    for (int k = 0; k < 3; ++k) {
      const double L = std::round(s_current[k] - s_target[k]);
      s_new[k] = s_target[k] + L;
    }
    domain->lamda2x(s_new, x[i]);
  }
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   verify the simulation box matches the lattice family declared in the
   symmetry file. The fix does NOT modify the box -- the user is expected
   to pick consistent NPT coupling (iso/aniso/tri); this check just
   catches an inconsistent initial setup at the earliest opportunity.

   The required relationships, in LAMMPS restricted-triclinic conventions
   (xprd along x, b-vector at xy tilt, c-vector at xz/yz tilts):

     triclinic     anything
     monoclinic    triclinic, xz = yz = 0 (only xy nonzero; b unique axis)
     orthorhombic  orthogonal (or triclinic with all tilts zero)
     tetragonal    orthorhombic + xprd == yprd
     cubic         orthorhombic + xprd == yprd == zprd
     hexagonal     triclinic, xz = yz = 0, xy = -xprd/2, yprd = xprd*sqrt(3)/2
     trigonal      hexagonal-style box (rhombohedral primitive is not auto-detected)

   Tolerance is the user-supplied `tol` scaled by the longest box edge,
   so the check tracks the same fractional-coordinate slack used to
   validate atom positions in build_orbit_map().
------------------------------------------------------------------------- */

void FixSymmetry::validate_box()
{
  const double xprd = domain->xprd;
  const double yprd = domain->yprd;
  const double zprd = domain->zprd;
  const double xy = domain->xy;
  const double xz = domain->xz;
  const double yz = domain->yz;
  const int tri = domain->triclinic;

  const double box_tol = tol * std::max({xprd, yprd, zprd});
  const LatticeFamily lat = grp->lattice();

  auto orthogonal_or_zero_tilt = [&]() {
    if (tri == 0) return true;
    return std::fabs(xy) < box_tol && std::fabs(xz) < box_tol && std::fabs(yz) < box_tol;
  };

  switch (lat) {

    case LATTICE_TRICLINIC:
      // any box is acceptable
      break;

    case LATTICE_MONOCLINIC:
      if (tri == 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix symmetry: monoclinic lattice requires a triclinic box "
                   "(create_box ... or triclinic data file)");
      if (std::fabs(xz) > box_tol || std::fabs(yz) > box_tol)
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix symmetry: monoclinic lattice requires xz = yz = 0 "
                   "(only xy tilt is allowed; got xz = {:.6g}, yz = {:.6g})",
                   xz, yz);
      break;

    case LATTICE_ORTHORHOMBIC:
      if (!orthogonal_or_zero_tilt())
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix symmetry: orthorhombic lattice requires an orthogonal box "
                   "(or triclinic with all tilts = 0; got xy = {:.6g}, xz = {:.6g}, yz = {:.6g})",
                   xy, xz, yz);
      break;

    case LATTICE_TETRAGONAL:
      if (!orthogonal_or_zero_tilt())
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix symmetry: tetragonal lattice requires an orthogonal box "
                   "(or triclinic with all tilts = 0)");
      if (std::fabs(xprd - yprd) > box_tol)
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix symmetry: tetragonal lattice requires xprd == yprd "
                   "(got {:.6g} vs {:.6g})",
                   xprd, yprd);
      break;

    case LATTICE_CUBIC:
      if (!orthogonal_or_zero_tilt())
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix symmetry: cubic lattice requires an orthogonal box "
                   "(or triclinic with all tilts = 0)");
      if (std::fabs(xprd - yprd) > box_tol || std::fabs(yprd - zprd) > box_tol)
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix symmetry: cubic lattice requires xprd == yprd == zprd "
                   "(got {:.6g}, {:.6g}, {:.6g})",
                   xprd, yprd, zprd);
      break;

    case LATTICE_HEXAGONAL:
    case LATTICE_TRIGONAL: {
      // Both are validated against a hexagonal-style triclinic cell:
      // a along x, b at 120 deg from a, c along z.
      if (tri == 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix symmetry: {} lattice requires a triclinic box with "
                   "xy = -xprd/2, yprd = xprd * sqrt(3)/2",
                   lat == LATTICE_HEXAGONAL ? "hexagonal" : "trigonal");
      if (std::fabs(xz) > box_tol || std::fabs(yz) > box_tol)
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix symmetry: {} lattice requires xz = yz = 0 "
                   "(got xz = {:.6g}, yz = {:.6g})",
                   lat == LATTICE_HEXAGONAL ? "hexagonal" : "trigonal", xz, yz);
      const double xy_expect = -0.5 * xprd;
      if (std::fabs(xy - xy_expect) > box_tol)
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix symmetry: {} lattice requires xy = -xprd/2 "
                   "(expected {:.6g}, got {:.6g})",
                   lat == LATTICE_HEXAGONAL ? "hexagonal" : "trigonal", xy_expect, xy);
      const double yprd_expect = xprd * std::sqrt(3.0) / 2.0;
      if (std::fabs(yprd - yprd_expect) > box_tol)
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix symmetry: {} lattice requires yprd = xprd * sqrt(3)/2 "
                   "(expected {:.6g}, got {:.6g})",
                   lat == LATTICE_HEXAGONAL ? "hexagonal" : "trigonal", yprd_expect, yprd);
      break;
    }
  }
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   invert the SymmetryGroup orbit table into the per-tag lookup
   tag_info[tag] = (orbit_id, op_id), then verify every declared tag is
   actually owned by some rank in the current system. Detects:
     - duplicate tag assignments (same tag listed in two orbits, or twice
       in one orbit -- which the symmetry-file group_op<->tag schema can
       in principle express even though it shouldn't)
     - declared tags that don't exist in the system
------------------------------------------------------------------------- */

void FixSymmetry::build_orbit_map()
{
  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, Error::NOLASTLINE, "Fix symmetry requires an atom map, see atom_modify");

  tag_info.clear();

  const int n_orbits_grp = grp->n_orbits();
  const int n_ops_grp = grp->n_ops();

  // pass 1: build the (tag -> (orbit, op)) lookup, detecting any duplicate
  // tag assignment between or within orbits. For atoms on Wyckoff special
  // positions, the orbit's image_tag[k] for k in stabilizer == asym_tag by
  // construction; skip those entries here so they don't trip the duplicate
  // check. post_force re-injects them via Orbit::stabilizer when needed.
  for (int o = 0; o < n_orbits_grp; ++o) {
    const Orbit &orb = grp->orbit(o);
    std::unordered_set<int> stab_set(orb.stabilizer.begin(), orb.stabilizer.end());
    for (int k = 0; k < n_ops_grp; ++k) {
      if (k != 0 && stab_set.count(k)) continue;
      tagint t = orb.image_tag[k];
      auto it = tag_info.find(t);
      if (it != tag_info.end()) {
        const auto &prev = it->second;
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix symmetry: atom tag {} is listed in two orbit positions "
                   "(orbit {} op {} and orbit {} op {})",
                   t, prev.first + 1, prev.second + 1, o + 1, k + 1);
      }
      tag_info.emplace(t, std::make_pair(o, k));
    }
  }

  // pass 2: confirm every declared tag is owned by exactly one rank.
  // We do NOT use atom->map() here -- for a tag larger than the system's
  // tag_max, map() reads out of bounds and may return a plausible-looking
  // local index, silently passing a missing-tag check. Instead, iterate
  // local atoms and mark which declared tags we found; the SUM across
  // ranks tells us how many owners exist (1 = good, 0 = missing,
  // >1 = duplicate, which should never happen in well-formed LAMMPS).
  std::vector<tagint> all_tags;
  all_tags.reserve(tag_info.size());
  for (const auto &kv : tag_info) all_tags.push_back(kv.first);
  std::sort(all_tags.begin(), all_tags.end());

  std::unordered_map<tagint, size_t> tag_to_idx;
  for (size_t i = 0; i < all_tags.size(); ++i) tag_to_idx[all_tags[i]] = i;

  std::vector<int> local_seen(all_tags.size(), 0);
  std::vector<int> global_seen(all_tags.size(), 0);
  for (int i = 0; i < atom->nlocal; ++i) {
    auto it = tag_to_idx.find(atom->tag[i]);
    if (it != tag_to_idx.end()) local_seen[it->second] = 1;
  }
  if (!all_tags.empty())
    MPI_Allreduce(local_seen.data(), global_seen.data(), all_tags.size(), MPI_INT, MPI_SUM, world);

  for (size_t i = 0; i < all_tags.size(); ++i) {
    if (global_seen[i] == 0) {
      const auto &info = tag_info[all_tags[i]];
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix symmetry: atom tag {} (orbit {} op {}) "
                 "does not exist in the system",
                 all_tags[i], info.first + 1, info.second + 1);
    } else if (global_seen[i] > 1) {
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix symmetry: atom tag {} is owned by more than one rank "
                 "(internal error, please report)",
                 all_tags[i]);
    }
  }

  // pass 3: verify the user's initial configuration is actually consistent
  // with the declared symmetry. For each orbit, gather the asymmetric atom's
  // fractional position (only the owning rank writes nonzero), then on
  // every rank check that each locally owned image atom sits at
  //   R_op * s_asym + t_op + (integer lattice translation)
  // to within `tol`. Catches mismatched orbit-to-tag assignments and
  // typos in the symmetry file before the first integration step.
  std::vector<double> s_asym_local(3 * n_orbits_grp, 0.0);
  std::vector<double> s_asym_global(3 * n_orbits_grp, 0.0);
  for (int i = 0; i < atom->nlocal; ++i) {
    auto it = tag_info.find(atom->tag[i]);
    if (it == tag_info.end()) continue;
    if (it->second.second != 0) continue;
    double s[3];
    domain->x2lamda(atom->x[i], s);
    const int o = it->second.first;
    s_asym_local[3 * o + 0] = s[0];
    s_asym_local[3 * o + 1] = s[1];
    s_asym_local[3 * o + 2] = s[2];
  }
  if (n_orbits_grp > 0)
    MPI_Allreduce(s_asym_local.data(), s_asym_global.data(), 3 * n_orbits_grp, MPI_DOUBLE, MPI_SUM,
                  world);

  for (int i = 0; i < atom->nlocal; ++i) {
    auto it = tag_info.find(atom->tag[i]);
    if (it == tag_info.end()) continue;
    const int op_id = it->second.second;
    if (op_id == 0) continue;
    const int o = it->second.first;
    const SymmOp &op = grp->op(op_id);

    const double s_asym[3] = {s_asym_global[3 * o + 0], s_asym_global[3 * o + 1],
                              s_asym_global[3 * o + 2]};
    double s_expect[3] = {
        op.R[0][0] * s_asym[0] + op.R[0][1] * s_asym[1] + op.R[0][2] * s_asym[2] + op.t[0],
        op.R[1][0] * s_asym[0] + op.R[1][1] * s_asym[1] + op.R[1][2] * s_asym[2] + op.t[1],
        op.R[2][0] * s_asym[0] + op.R[2][1] * s_asym[1] + op.R[2][2] * s_asym[2] + op.t[2]};
    double s_actual[3];
    domain->x2lamda(atom->x[i], s_actual);

    double mismatch = 0.0;
    for (int k = 0; k < 3; ++k) {
      double d = s_actual[k] - s_expect[k];
      d -= std::round(d);    // discount any periodic-image offset
      mismatch = std::max(mismatch, std::fabs(d));
    }
    if (mismatch > tol) {
      const int asym_tag = grp->orbit(o).asym_tag;
      error->one(FLERR, Error::NOLASTLINE,
                 "Fix symmetry: initial position of atom {} (orbit {} op {}) "
                 "differs from R*s_asym + t by {:.3g} in fractional coords "
                 "(asym atom {}, tol {:.3g}). The orbit assignment in the "
                 "symmetry file does not match the actual atom layout.",
                 atom->tag[i], o + 1, op_id + 1, mismatch, asym_tag, tol);
    }
  }

  // pass 4: compute the per-orbit Wyckoff constant c. For each orbit with
  // a non-empty stabilizer, the constraint is
  //   (R_k - I) * s_asym + (t_k - L_k) == 0   for every k in stabilizer
  // where L_k is the integer lattice translation chosen by the user's
  // initial placement (we round to capture it). Define
  //   c_o = sum_{k in stabilizer_o} (R_k - I)^T * (t_k - L_k)
  // Then at runtime the projection residual is  A^T r = B * s + c_o
  // and the position correction is  delta_s = Binv * (B * s + c_o).
  // Also verify that the residual at init is essentially zero -- if it
  // isn't, the user's initial asym atom is not actually on the declared
  // Wyckoff site and we should refuse to proceed.
  wyckoff_c.assign(3 * n_orbits_grp, 0.0);
  for (int o = 0; o < n_orbits_grp; ++o) {
    const Orbit &orb = grp->orbit(o);
    if (orb.stabilizer.empty()) continue;
    const double *s_asym = &s_asym_global[3 * o];
    for (int k : orb.stabilizer) {
      const SymmOp &op = grp->op(k);
      double delta[3];
      for (int i = 0; i < 3; ++i)
        delta[i] = op.R[i][0] * s_asym[0] + op.R[i][1] * s_asym[1] + op.R[i][2] * s_asym[2] +
            op.t[i] - s_asym[i];
      double L[3], res[3];
      for (int i = 0; i < 3; ++i) {
        L[i] = std::round(delta[i]);
        res[i] = delta[i] - L[i];
      }
      double max_res = std::max({std::fabs(res[0]), std::fabs(res[1]), std::fabs(res[2])});
      if (max_res > tol)
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix symmetry: asym atom {} of orbit {} is declared on Wyckoff "
                   "site (op {}) but residual |R*s + t - s| (mod 1) = {:.3g} "
                   "exceeds tol = {:.3g}",
                   orb.asym_tag, o + 1, k + 1, max_res, tol);
      const double tk_shifted[3] = {op.t[0] - L[0], op.t[1] - L[1], op.t[2] - L[2]};
      // c += (R_k - I)^T * tk_shifted
      for (int i = 0; i < 3; ++i) {
        double sum = 0.0;
        for (int j = 0; j < 3; ++j) {
          const double Rmi_T = op.R[j][i] - (j == i ? 1.0 : 0.0);
          sum += Rmi_T * tk_shifted[j];
        }
        wyckoff_c[3 * o + i] += sum;
      }
    }
  }
}
