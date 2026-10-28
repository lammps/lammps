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
#include "group.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "respa.h"
#include "symmetry_group.h"
#include "update.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using MathConst::MY_PI;

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

// position conversions like Domain::x2lamda() and Domain::lamda2x(), but inlined
inline void x2frac(double s[3], const double x[3], const double h_inv[6], const double boxlo[3])
{
  const double d0 = x[0] - boxlo[0];
  const double d1 = x[1] - boxlo[1];
  const double d2 = x[2] - boxlo[2];
  s[0] = h_inv[0] * d0 + h_inv[5] * d1 + h_inv[4] * d2;
  s[1] = h_inv[1] * d1 + h_inv[3] * d2;
  s[2] = h_inv[2] * d2;
}

inline void frac2x(double x[3], const double s[3], const double h[6], const double boxlo[3])
{
  x[0] = h[0] * s[0] + h[5] * s[1] + h[4] * s[2] + boxlo[0];
  x[1] = h[1] * s[1] + h[3] * s[2] + boxlo[1];
  x[2] = h[2] * s[2] + boxlo[2];
}

// out = M * in for a 3x3 matrix stored as 9 consecutive doubles (row-major)
inline void matvec(double out[3], const double *M, const double in[3])
{
  out[0] = M[0] * in[0] + M[1] * in[1] + M[2] * in[2];
  out[1] = M[3] * in[0] + M[4] * in[1] + M[5] * in[2];
  out[2] = M[6] * in[0] + M[7] * in[1] + M[8] * in[2];
}

// out = M * in for a 3x3 matrix M[i][j]
inline void matvec(double out[3], const double M[3][3], const double in[3])
{
  out[0] = M[0][0] * in[0] + M[0][1] * in[1] + M[0][2] * in[2];
  out[1] = M[1][0] * in[0] + M[1][1] * in[1] + M[1][2] * in[2];
  out[2] = M[2][0] * in[0] + M[2][1] * in[1] + M[2][2] * in[2];
}

// image of fractional position s under operator op: R s + t
inline void apply_op(double out[3], const SymmOp &op, const double s[3])
{
  matvec(out, op.R, s);
  out[0] += op.t[0];
  out[1] += op.t[1];
  out[2] += op.t[2];
}
}    // namespace

/* ---------------------------------------------------------------------- */

FixSymmetry::FixSymmetry(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), grp(nullptr), tol(1.0e-6), symfile(nullptr), n_orbits(0), n_special(0),
    n_trans(0), atom2member(nullptr), nmax_cache(0), shift_valid(false), buf_local(nullptr),
    buf_global(nullptr)
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

  dof_flag = 1;

  grp = new SymmetryGroup(lmp);
  grp->read(symfile);

  // in 2d systems the operators must leave the z coordinate unchanged
  if (domain->dimension == 2) {
    for (int k = 1; k < grp->n_ops(); ++k) {
      const auto &op = grp->op(k);
      if ((op.R[0][2] != 0.0) || (op.R[1][2] != 0.0) || (op.R[2][0] != 0.0) ||
          (op.R[2][1] != 0.0) || (op.R[2][2] != 1.0) ||
          (std::fabs(op.t[2] - std::round(op.t[2])) > tol))
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix symmetry: symmetry operator {} changes the z coordinate, which is not "
                   "possible for a 2d system",
                   k + 1);
    }
  }

  // number of independent translations of the whole system: the dimension
  // of the space of vectors left unchanged by all operators
  double trace = 0.0;
  const int dim = domain->dimension;
  for (int k = 0; k < grp->n_ops(); ++k)
    for (int d = 0; d < dim; ++d) trace += grp->op(k).R[d][d];
  n_trans = static_cast<int>(std::lround(trace / grp->n_ops()));

  for (double &hval : h_cart) hval = 0.0;

  if (lmp->citeme) lmp->citeme->add(cite_symd);
}

/* ---------------------------------------------------------------------- */

FixSymmetry::~FixSymmetry()
{
  delete grp;
  delete[] symfile;
  memory->destroy(atom2member);
  memory->destroy(buf_local);
  memory->destroy(buf_global);
}

/* ---------------------------------------------------------------------- */

int FixSymmetry::setmask()
{
  int mask = 0;
  mask |= POST_NEIGHBOR;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_NEIGHBOR;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSymmetry::init()
{
  // check the box against the operators
  mcart.clear();
  update_cartesian_ops();

  detect_orbits();

  memory->destroy(buf_local);
  memory->destroy(buf_global);
  memory->create(buf_local, 9 * n_orbits + 1, "symmetry:buf_local");
  memory->create(buf_global, 9 * n_orbits + 1, "symmetry:buf_global");

  warn_conflicting_fixes();
  check_fix_order();
  print_summary();
}

/* ----------------------------------------------------------------------
   Forces added by fixes that are applied after this fix are not
   symmetrized and will break the symmetry.
------------------------------------------------------------------------- */

void FixSymmetry::check_fix_order()
{
  bool after = false;
  for (auto *ifix : modify->get_fix_list()) {
    if (ifix == this) {
      after = true;
    } else if (after && (modify->get_fix_mask(ifix) & (POST_FORCE | POST_FORCE_RESPA))) {
      if (comm->me == 0)
        error->warning(FLERR,
                       "Fix {} (style {}) may change forces after fix symmetry {}; fix symmetry "
                       "must be defined after all fixes that add or change forces",
                       ifix->id, ifix->style, id);
    }
  }
}

/* ----------------------------------------------------------------------
   Print a summary of the orbits when it is new or has changed
------------------------------------------------------------------------- */

void FixSymmetry::print_summary()
{
  bigint nconstrained = 0, nfree = 0;
  for (const auto &orb : orbits) {
    nconstrained += orb.nmembers;
    nfree += orb.fixed_dim;
  }
  auto mesg = fmt::format("Fix symmetry {}: group {} with {} operators, {} atoms in {} orbits "
                          "({} on special positions) with {} independent coordinates\n",
                          id, grp->name(), grp->n_ops(), nconstrained, n_orbits, n_special, nfree);
  if (mesg != summary) {
    summary = mesg;
    if (comm->me == 0) utils::logmesg(lmp, summary);
  }
}

/* ----------------------------------------------------------------------
   Warn if any atom constrained by this fix is also controlled by a fix
   that enforces its own rigid-body or fixed bond/angle geometry (rigid,
   shake, rattle, ilves). Both constraint schemes independently overwrite
   the forces, velocities, and/or positions of the shared atoms every step
   with no coordination, so they interfere. This is only a problem when the
   two act on overlapping atoms, so we check for an actual shared atom (an
   orbit atom that is also in the other fix's group) rather than the mere
   presence of such a fix.
------------------------------------------------------------------------- */

void FixSymmetry::warn_conflicting_fixes()
{
  const int *mask = atom->mask;
  const tagint *tag = atom->tag;
  const int nlocal = atom->nlocal;

  for (auto *fix : modify->get_fix_list()) {
    if (fix == this) continue;
    std::string s = fix->style;
    if (!(utils::strmatch(s, "^rigid") || utils::strmatch(s, "^shake") ||
          utils::strmatch(s, "^rattle") || utils::strmatch(s, "^ilves")))
      continue;

    // count local orbit atoms that are also in the other fix's group
    const int gbit = fix->groupbit;
    bigint nshared = 0;
    for (int i = 0; i < nlocal; ++i)
      if ((mask[i] & gbit) && (member_of_tag.find(tag[i]) != member_of_tag.end())) ++nshared;

    bigint nshared_all = 0;
    MPI_Allreduce(&nshared, &nshared_all, 1, MPI_LMP_BIGINT, MPI_SUM, world);

    if ((nshared_all > 0) && (comm->me == 0))
      error->warning(FLERR,
                     "Fix symmetry and fix {} (style {}) both constrain {} shared "
                     "atom(s); the two constraint schemes are applied independently "
                     "and will interfere. See the fix symmetry doc page.",
                     fix->id, s, nshared_all);
  }
}

/* ---------------------------------------------------------------------- */

void FixSymmetry::setup(int vflag)
{
  build_cache();
  update_cartesian_ops();

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
}

/* ----------------------------------------------------------------------
   local atoms change only when reneighboring, so the member index of each
   local atom is looked up once here and not on every step
------------------------------------------------------------------------- */

void FixSymmetry::post_neighbor()
{
  build_cache();
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

/* ----------------------------------------------------------------------
   Symmetrize forces, velocities, and positions of all orbit atoms.

   With Cartesian operators M_k = h R_k h^-1 (orthogonal matrices when the
   box is compatible with the group), member m of an orbit is the image of
   the representative under its canonical operator k(m).  The constrained
   force on the representative is the average of the forces mapped back by
   the inverse operators, restricted to the free coordinates of the
   representative (projector P, the identity on a general position):

     F_rep = P * (1/nmembers) * sum_m  M_k(m)^-1 f_m
     f_m  <- M_k(m) F_rep

   The same is done for the velocities.  This is the exact projection of
   the forces onto the constraint, so from a symmetric starting point the
   trajectory is the same as without the fix, only free of round-off.

   Positions are reset from the representatives in the same pass.  This
   removes the round-off the integrator accumulates between members of an
   orbit.  The owner of a representative on a special position first moves
   it back onto that position: the average of its images under the site
   symmetry operators, each shifted by the lattice translation that brings
   it closest.  Every other member is placed at the image R_k s + t_k of the
   representative in the periodic copy closest to its current position.
   This preserves the periodic image the atom had at setup and is robust
   against atoms that moved outside the box between reneighborings.

   All data goes into one collective sum of 9 numbers per orbit.  Atoms not
   in any orbit are left untouched.
------------------------------------------------------------------------- */

void FixSymmetry::post_force(int /*vflag*/)
{
  symmetrize(true);
}

/* ----------------------------------------------------------------------
   Minimization: only the forces are symmetrized.  Positions are owned by
   the minimizer and move along the symmetric forces.
------------------------------------------------------------------------- */

void FixSymmetry::min_setup(int /*vflag*/)
{
  build_cache();
  update_cartesian_ops();
  symmetrize(false);
}

/* ---------------------------------------------------------------------- */

void FixSymmetry::min_post_neighbor()
{
  build_cache();
}

/* ---------------------------------------------------------------------- */

void FixSymmetry::min_post_force(int /*vflag*/)
{
  symmetrize(false);
}

/* ----------------------------------------------------------------------
   Symmetrize the forces and, if full is true, the velocities and positions
------------------------------------------------------------------------- */

void FixSymmetry::symmetrize(bool full)
{
  if (n_orbits == 0) return;
  update_cartesian_ops();

  const int nlocal = atom->nlocal;
  double **x = atom->x;
  double **f = atom->f;
  double **v = atom->v;
  const double *h = domain->h;
  const double *h_inv = domain->h_inv;
  const double *boxlo = domain->boxlo;
  double s[3], image[3];

  std::fill_n(buf_local, 9 * n_orbits, 0.0);

  // pass 1: map forces and velocities back to the representative,
  // collect the positions of the representatives

  for (int i = 0; i < nlocal; ++i) {
    const int m = atom2member[i];
    if (m < 0) continue;
    const int k = member_op[m];
    const double *minv = &mcart_inv[9 * k];
    double *sum = buf_local + 9 * member_orbit[m];
    double r[3];
    matvec(r, minv, f[i]);
    sum[0] += r[0];
    sum[1] += r[1];
    sum[2] += r[2];
    if (!full) continue;
    matvec(r, minv, v[i]);
    sum[3] += r[0];
    sum[4] += r[1];
    sum[5] += r[2];

    if (k == 0) {
      const auto &orb = orbits[member_orbit[m]];
      x2frac(s, x[i], h_inv, boxlo);
      if (orb.special) {
        double avg[3] = {s[0], s[1], s[2]};
        for (int kk : orb.stabilizer) {
          apply_op(image, grp->op(kk), s);
          for (int d = 0; d < 3; ++d) avg[d] += image[d] - std::round(image[d] - s[d]);
        }
        const double scale = 1.0 / static_cast<double>(orb.stabilizer.size() + 1);
        for (int d = 0; d < 3; ++d) s[d] = avg[d] * scale;
        frac2x(x[i], s, h, boxlo);
      }
      sum[6] = s[0];
      sum[7] = s[1];
      sum[8] = s[2];
    }
  }

  MPI_Allreduce(buf_local, buf_global, 9 * n_orbits, MPI_DOUBLE, MPI_SUM, world);

  // pass 2: assign the images of the averaged force and velocity and
  // of the position of the representative.  The lattice translation that
  // selects the periodic copy of each image can only change when atoms are
  // remapped into the box during reneighboring, so it is computed only on
  // the first call after that.

  for (int i = 0; i < nlocal; ++i) {
    const int m = atom2member[i];
    if (m < 0) continue;
    const int k = member_op[m];
    const auto &orb = orbits[member_orbit[m]];
    const double *sum = buf_global + 9 * member_orbit[m];
    const double scale = 1.0 / static_cast<double>(orb.nmembers);
    double frep[3] = {sum[0] * scale, sum[1] * scale, sum[2] * scale};
    double vrep[3] = {sum[3] * scale, sum[4] * scale, sum[5] * scale};

    // remove components the special position does not allow
    if (orb.special) {
      double frac[3], proj[3];
      cart_to_frac_vec(frac, frep, h_inv);
      matvec(proj, orb.proj, frac);
      frac_to_cart_vec(frep, proj, h);
      cart_to_frac_vec(frac, vrep, h_inv);
      matvec(proj, orb.proj, frac);
      frac_to_cart_vec(vrep, proj, h);
    }

    const double *mk = &mcart[9 * k];
    matvec(f[i], mk, frep);
    if (!full) continue;
    matvec(v[i], mk, vrep);

    if (k != 0) {
      double *shift = &image_shift[3 * i];
      apply_op(image, grp->op(k), sum + 6);
      if (!shift_valid) {
        x2frac(s, x[i], h_inv, boxlo);
        for (int d = 0; d < 3; ++d) shift[d] = std::round(s[d] - image[d]);
      }
      for (int d = 0; d < 3; ++d) image[d] += shift[d];
      frac2x(x[i], image, h, boxlo);
    }
  }
  if (full) shift_valid = true;
}

/* ----------------------------------------------------------------------
   Degrees of freedom removed by the symmetry constraints from the atoms in
   group tgroup.  An orbit of n atoms has as many degrees of freedom as its
   representative has free coordinates, so dimension*n - free are removed
   if the representative is in the group.  When the group consists of
   exactly the orbit atoms, the motion of the center of mass is restricted
   to the translations allowed by the operators (none if the group
   contains an inversion), so the dimension - n_trans degrees of freedom
   of the center of mass the temperature compute removes are added back.
------------------------------------------------------------------------- */

bigint FixSymmetry::dof(int tgroup)
{
  if (n_orbits == 0) return 0;
  const int gbit = group->bitmask[tgroup];
  const int dim = domain->dimension;
  const int *mask = atom->mask;
  const tagint *tag = atom->tag;
  const int nlocal = atom->nlocal;

  bigint counts[3] = {0, 0, 0}, counts_all[3] = {0, 0, 0};
  for (int i = 0; i < nlocal; ++i) {
    auto it = member_of_tag.find(tag[i]);
    if (mask[i] & gbit) {
      if (it == member_of_tag.end()) {
        ++counts[1];    // atom in group but not constrained
      } else if (member_op[it->second] == 0) {
        const auto &orb = orbits[member_orbit[it->second]];
        counts[0] += dim * orb.nmembers - orb.fixed_dim;
      }
    } else if (it != member_of_tag.end()) {
      ++counts[2];    // constrained atom not in group
    }
  }
  MPI_Allreduce(counts, counts_all, 3, MPI_LMP_BIGINT, MPI_SUM, world);

  bigint nremoved = counts_all[0];
  if ((counts_all[1] == 0) && (counts_all[2] == 0)) nremoved -= dim - n_trans;
  return nremoved;
}

/* ----------------------------------------------------------------------
   Look up the orbit member index of all local atoms.  Also verify that
   no orbit atom was lost or deleted since the orbits were set up.
------------------------------------------------------------------------- */

void FixSymmetry::build_cache()
{
  if (atom->nmax > nmax_cache) {
    nmax_cache = atom->nmax;
    memory->destroy(atom2member);
    memory->create(atom2member, nmax_cache, "symmetry:atom2member");
    image_shift.resize(3 * nmax_cache);
  }
  shift_valid = false;

  const int nlocal = atom->nlocal;
  const tagint *tag = atom->tag;
  bigint nfound = 0;
  for (int i = 0; i < nlocal; ++i) {
    auto it = member_of_tag.find(tag[i]);
    if (it == member_of_tag.end()) {
      atom2member[i] = -1;
    } else {
      atom2member[i] = it->second;
      ++nfound;
    }
  }

  bigint nfound_all = 0;
  MPI_Allreduce(&nfound, &nfound_all, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (nfound_all != static_cast<bigint>(member_orbit.size()))
    error->all(FLERR, Error::NOLASTLINE,
               "Fix symmetry {} found {} of its {} atoms; atoms were lost or deleted", id,
               nfound_all, member_orbit.size());
}

/* ----------------------------------------------------------------------
   Compute the Cartesian operators M_k = h R_k h^-1 and their inverses.
   Only done when the box has changed.
------------------------------------------------------------------------- */

void FixSymmetry::update_cartesian_ops()
{
  const double *h = domain->h;
  if (!mcart.empty() && (memcmp(h, h_cart, sizeof(h_cart)) == 0)) return;
  memcpy(h_cart, h, sizeof(h_cart));

  const double *hi = domain->h_inv;
  const double H[3][3] = {{h[0], h[5], h[4]}, {0.0, h[1], h[3]}, {0.0, 0.0, h[2]}};
  const double Hinv[3][3] = {{hi[0], hi[5], hi[4]}, {0.0, hi[1], hi[3]}, {0.0, 0.0, hi[2]}};

  const int nops = grp->n_ops();
  mcart.resize(9 * nops);
  mcart_inv.resize(9 * nops);

  auto transform = [&](double *out, const double R[3][3]) {
    double tmp[3][3];
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        tmp[i][j] = R[i][0] * Hinv[0][j] + R[i][1] * Hinv[1][j] + R[i][2] * Hinv[2][j];
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        out[3 * i + j] = H[i][0] * tmp[0][j] + H[i][1] * tmp[1][j] + H[i][2] * tmp[2][j];
  };

  for (int k = 0; k < nops; ++k) {
    transform(&mcart[9 * k], grp->op(k).R);
    transform(&mcart_inv[9 * k], grp->op(k).Rinv);
  }

  check_metric(H);
}

/* ----------------------------------------------------------------------
   The operators are symmetries of the lattice only if they preserve the
   lengths of and angles between the box vectors (the columns of H), i.e.
   R^T G R = G with the metric G = H^T H.  Otherwise the Cartesian
   operators are not rotations or reflections and the symmetrized forces
   would not conserve energy.
------------------------------------------------------------------------- */

void FixSymmetry::check_metric(const double H[3][3])
{
  double G[3][3], gmax = 0.0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      G[i][j] = H[0][i] * H[0][j] + H[1][i] * H[1][j] + H[2][i] * H[2][j];
      gmax = std::max(gmax, std::fabs(G[i][j]));
    }

  for (int k = 1; k < grp->n_ops(); ++k) {
    const auto &R = grp->op(k).R;
    double dev = 0.0;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        double rgr = 0.0;
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b) rgr += R[a][i] * G[a][b] * R[b][j];
        dev = std::max(dev, std::fabs(rgr - G[i][j]));
      }
    if (dev > tol * gmax) {
      const double a = sqrt(G[0][0]), b = sqrt(G[1][1]), c = sqrt(G[2][2]);
      const double rad2deg = 180.0 / MY_PI;
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix symmetry: symmetry operator {} is not compatible with the simulation box: "
                 "a = {:.8g}, b = {:.8g}, c = {:.8g}, alpha = {:.8g}, beta = {:.8g}, "
                 "gamma = {:.8g}",
                 k + 1, a, b, c, acos(G[1][2] / (b * c)) * rad2deg,
                 acos(G[0][2] / (a * c)) * rad2deg, acos(G[0][1] / (a * b)) * rad2deg);
    }
  }
}

/* ----------------------------------------------------------------------
   Find the symmetry orbits of the atoms in the fix group.

   All ranks gather the fractional positions, types, and masses of all
   group atoms and do the same deterministic search, so the result is
   identical everywhere without further communication.  Atoms are visited
   in the order of their IDs.  Each atom not yet in an orbit becomes the
   representative of a new orbit, and the atom at its image under each
   operator is looked up with a grid of cells with a size of a few tol.
   The image must be an atom of the fix group of the same type and mass
   within tol in fractional coordinates.  An atom on a special position is
   its own image under the operators of its site symmetry, so its orbit has
   fewer members than the group has operators.
------------------------------------------------------------------------- */

void FixSymmetry::detect_orbits()
{
  const int nlocal = atom->nlocal;
  const int *mask = atom->mask;
  const tagint *tag = atom->tag;
  const int *type = atom->type;
  double **x = atom->x;
  const int nops = grp->n_ops();

  orbits.clear();
  n_orbits = n_special = 0;
  member_of_tag.clear();
  member_orbit.clear();
  member_op.clear();

  // a group with only the identity does not constrain anything

  if (nops == 1) return;

  // gather IDs, fractional positions, types, and masses of all group atoms

  constexpr int NDATA = 5;
  int nmine = 0;
  for (int i = 0; i < nlocal; ++i)
    if (mask[i] & groupbit) ++nmine;

  const bigint nall_big = group->count(igroup);
  if (nall_big > MAXSMALLINT / NDATA)
    error->all(FLERR, Error::NOLASTLINE, "Fix symmetry {}: too many atoms in fix group", id);
  const int nall = static_cast<int>(nall_big);

  std::vector<tagint> mytags, alltags(nall);
  std::vector<double> mydata, alldata(NDATA * nall);
  mytags.reserve(nmine);
  mydata.reserve(NDATA * nmine);
  for (int i = 0; i < nlocal; ++i) {
    if (!(mask[i] & groupbit)) continue;
    double s[3];
    domain->x2lamda(x[i], s);
    mytags.push_back(tag[i]);
    for (double sval : s) mydata.push_back(sval - std::floor(sval));
    mydata.push_back(type[i]);
    mydata.push_back(atom->rmass ? atom->rmass[i] : atom->mass[type[i]]);
  }

  const int nprocs = comm->nprocs;
  std::vector<int> counts(nprocs), displs(nprocs, 0);
  MPI_Allgather(&nmine, 1, MPI_INT, counts.data(), 1, MPI_INT, world);
  for (int p = 1; p < nprocs; ++p) displs[p] = displs[p - 1] + counts[p - 1];
  MPI_Allgatherv(mytags.data(), nmine, MPI_LMP_TAGINT, alltags.data(), counts.data(), displs.data(),
                 MPI_LMP_TAGINT, world);
  for (int p = 0; p < nprocs; ++p) {
    counts[p] *= NDATA;
    displs[p] *= NDATA;
  }
  MPI_Allgatherv(mydata.data(), NDATA * nmine, MPI_DOUBLE, alldata.data(), counts.data(),
                 displs.data(), MPI_DOUBLE, world);
  mytags.clear();
  mydata.clear();

  // order of atom IDs

  std::vector<int> order(nall);
  for (int j = 0; j < nall; ++j) order[j] = j;
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    return alltags[a] < alltags[b];
  });

  // sort atoms into cells of a grid in fractional coordinates.  A cell is
  // larger than tol, so all atoms within tol of a position are in the same
  // or a neighboring cell.

  const int ncell = static_cast<int>(std::max(1.0, std::min(1048576.0, std::floor(0.25 / tol))));
  auto cell_of = [ncell](double sval) {
    return std::min(static_cast<int>(sval * ncell), ncell - 1);
  };
  auto cell_key = [ncell](int a, int b, int c) {
    return (static_cast<int64_t>(a) * ncell + b) * ncell + c;
  };
  std::unordered_map<int64_t, std::vector<int>> cells;
  for (int j = 0; j < nall; ++j) {
    const double *sj = &alldata[NDATA * j];
    cells[cell_key(cell_of(sj[0]), cell_of(sj[1]), cell_of(sj[2]))].push_back(j);
  }

  // neighbor cell offsets, without duplicates for grids with fewer than 3 cells
  std::vector<int> offsets = {0};
  if (ncell > 1) offsets.push_back(1);
  if (ncell > 2) offsets.push_back(-1);

  // index of the atom within tol of fractional position s, -1 if none, -2 if more than one
  auto find_atom = [&](const double *s) {
    int found = -1;
    const int c[3] = {cell_of(s[0]), cell_of(s[1]), cell_of(s[2])};
    for (int da : offsets)
      for (int db : offsets)
        for (int dc : offsets) {
          auto it = cells.find(cell_key((c[0] + da + ncell) % ncell, (c[1] + db + ncell) % ncell,
                                        (c[2] + dc + ncell) % ncell));
          if (it == cells.end()) continue;
          for (int j : it->second) {
            const double *sj = &alldata[NDATA * j];
            bool close = true;
            for (int d = 0; d < 3; ++d) {
              const double delta = sj[d] - s[d];
              if (std::fabs(delta - std::round(delta)) > tol) close = false;
            }
            if (close) found = (found == -1) ? j : -2;
          }
        }
    return found;
  };

  // build the orbits

  std::vector<int> orbit_of(nall, -1);
  std::vector<tagint> rep_tag;
  for (int rep : order) {
    if (orbit_of[rep] >= 0) continue;
    const int o = n_orbits++;
    rep_tag.push_back(alltags[rep]);
    orbits.emplace_back();
    auto &orb = orbits.back();
    const double *srep = &alldata[NDATA * rep];
    const int first = static_cast<int>(member_orbit.size());
    std::vector<int> count;

    for (int k = 0; k < nops; ++k) {
      double image[3];
      apply_op(image, grp->op(k), srep);
      for (double &sval : image) sval -= std::floor(sval);
      const int j = find_atom(image);
      if (j == -1)
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix symmetry {}: no atom of group {} at the image of atom ID {} under symmetry "
                   "operator {} (fractional coordinates {:.8g} {:.8g} {:.8g}, tolerance {:.3g})",
                   id, group->names[igroup], alltags[rep], k + 1, image[0], image[1], image[2],
                   tol);
      if (j == -2)
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix symmetry {}: more than one atom at the image of atom ID {} under symmetry "
                   "operator {} within tolerance {:.3g}",
                   id, alltags[rep], k + 1, tol);

      const double *sj = &alldata[NDATA * j];
      if (sj[3] != srep[3])
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix symmetry {}: atom ID {} at the image of atom ID {} under symmetry operator "
                   "{} has a different atom type: {} instead of {}",
                   id, alltags[j], alltags[rep], k + 1, static_cast<int>(sj[3]),
                   static_cast<int>(srep[3]));
      if (std::fabs(sj[4] - srep[4]) > 1.0e-10 * std::fabs(srep[4]))
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix symmetry {}: atom ID {} at the image of atom ID {} under symmetry operator "
                   "{} has a different mass: {} instead of {}",
                   id, alltags[j], alltags[rep], k + 1, sj[4], srep[4]);

      if (orbit_of[j] < 0) {
        orbit_of[j] = o;
        member_of_tag[alltags[j]] = static_cast<int>(member_orbit.size());
        member_orbit.push_back(o);
        member_op.push_back(k);
        count.push_back(1);
      } else if (orbit_of[j] == o) {
        ++count[member_of_tag[alltags[j]] - first];
      } else {
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix symmetry {}: atom ID {} at the image of atom ID {} under symmetry operator "
                   "{} was already found as an image of atom ID {}; the positions are not "
                   "symmetric within tolerance {:.3g}",
                   id, alltags[j], alltags[rep], k + 1, rep_tag[orbit_of[j]], tol);
      }
      if ((k > 0) && (j == rep)) orb.stabilizer.push_back(k);
    }

    // every atom of an orbit is the image of the representative under the
    // same number of operators: the size of the site symmetry group
    const int nstab = static_cast<int>(orb.stabilizer.size()) + 1;
    for (int c : count)
      if (c != nstab)
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix symmetry {}: the images of atom ID {} are not symmetric within tolerance "
                   "{:.3g}",
                   id, alltags[rep], tol);

    orb.nmembers = static_cast<int>(count.size());
    orb.special = (nstab > 1);
    if (orb.special) ++n_special;

    // projector onto the coordinates left free by the site symmetry: the
    // average of the site symmetry operators.  Its trace is the number of
    // free coordinates.
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) orb.proj[i][j] = grp->op(0).R[i][j];
    for (int k : orb.stabilizer)
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) orb.proj[i][j] += grp->op(k).R[i][j];
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) orb.proj[i][j] /= static_cast<double>(nstab);
    double trace = orb.proj[0][0] + orb.proj[1][1];
    if (domain->dimension == 3) trace += orb.proj[2][2];
    orb.fixed_dim = static_cast<int>(std::lround(trace));
  }
}
