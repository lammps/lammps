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
   Local-topology variant of the ILVES bond/angle constraint solver.

   Discovers constraint clusters purely from each rank's own bond/angle
   storage.  No MPI_Allgatherv of the global topology, no natoms-sized
   tag-to-cluster table.  See fix_ilves_local.h for the restrictions
   this introduces.

   Cluster ownership rule:
     * A "star center" is a local atom with > 1 selected constrained bond
       in its local bond list.  Under newton_bond=off (required by this
       variant) all of a center's bonds are visible locally, so a star
       center on this rank fully describes the cluster.
     * The star center's rank owns the cluster.  Leaf-only atoms (count
       == 1) on a rank are silent unless the cluster is a 1+1 pair (no
       center), in which case the rank holding the lower-tag atom owns.
   This mirrors the "rank holding the central atom processes the cluster"
   pattern used by fix shake.  Star + 1+1 cover water HOH, methyl/methane
   groups, isolated single bonds -- the canonical ILVES use cases.
------------------------------------------------------------------------- */

#include "fix_ilves_local.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "update.h"

#include <cmath>
#include <vector>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixIlvesLocal::FixIlvesLocal(LAMMPS *lmp, int narg, char **arg) :
    FixIlves(lmp, narg, arg), rbuf(nullptr), maxbuf(0)
{
  // 3 doubles per atom for reverse_comm of force / velocity contributions
  comm_reverse = 3;
  grow_rbuf(atom->nmax);
}

FixIlvesLocal::~FixIlvesLocal()
{
  memory->destroy(rbuf);
}

double FixIlvesLocal::memory_usage()
{
  return FixIlves::memory_usage() + (double) maxbuf * 3 * sizeof(double);
}

void FixIlvesLocal::grow_arrays(int nmax)
{
  FixIlves::grow_arrays(nmax);
  grow_rbuf(nmax);
}

void FixIlvesLocal::grow_rbuf(int nmax)
{
  if (nmax <= maxbuf) return;
  maxbuf = nmax;
  memory->grow(rbuf, maxbuf, 3, "ilves/local:rbuf");
}

/* ----------------------------------------------------------------------
   reverse_comm pack/unpack: send rbuf[ghost] contribution from this rank
   to the owner's local atom, where it accumulates into rbuf[local].
------------------------------------------------------------------------- */

int FixIlvesLocal::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; ++i) {
    buf[m++] = rbuf[i][0];
    buf[m++] = rbuf[i][1];
    buf[m++] = rbuf[i][2];
  }
  return m;
}

void FixIlvesLocal::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m = 0;
  for (int i = 0; i < n; ++i) {
    int j = list[i];
    rbuf[j][0] += buf[m++];
    rbuf[j][1] += buf[m++];
    rbuf[j][2] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   Count selected constrained bonds in local atom i's bond list.

   Reads bond_type / bond_atom directly, applies the b/t/m selectors,
   and requires both endpoints to be in the fix group.  Bond types may
   already have been negated by negate_constrained_topology() at this
   point, so abs() is applied before lookup.
------------------------------------------------------------------------- */

int FixIlvesLocal::count_constrained_bonds_local(int i)
{
  int **nb_type    = atom->bond_type;
  tagint **nb_atom = atom->bond_atom;
  int *num_bond    = atom->num_bond;
  int *mask        = atom->mask;

  if (!(mask[i] & groupbit)) return 0;

  int count = 0;
  for (int b = 0; b < num_bond[i]; ++b) {
    int bt = nb_type[i][b];
    if (bt == 0) continue;
    if (bt < 0) bt = -bt;
    if (bt > atom->nbondtypes) continue;

    int j = atom->map(nb_atom[i][b]);
    if (j < 0) continue;
    if (!(mask[j] & groupbit)) continue;
    if (!bond_selected_for_atoms(i, j, bt)) continue;
    ++count;
  }
  return count;
}

/* ----------------------------------------------------------------------
   Apply constraint forces for cluster-owner solves.

   The single-rank cluster ownership model requires that constraint
   contributions to ghost atoms on this rank (the cluster's "leaves")
   reach the owner rank of each ghost.  The reverse_comm machinery
   propagates this rank's rbuf[ghost_idx] -> owner's rbuf[local_idx],
   which we then accumulate into atom->f.

   PBC self-image ghosts share a tag with a local atom on this rank;
   atom->map(tag) returns the local index in that case, so writes land
   directly without a (degenerate) reverse-comm round-trip.
------------------------------------------------------------------------- */

void FixIlvesLocal::apply_constraint_forces(int vflag)
{
  // Allocate fstore lazily; mirrors FixIlves::apply_constraint_forces.
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

  // Zero rbuf for every atom slot (local + ghost) on this rank.
  for (int i = 0; i < atom->nmax; ++i)
    rbuf[i][0] = rbuf[i][1] = rbuf[i][2] = 0.0;

  const double inv_dtfsq = 1.0 / dtfsq;
  double bond_v[6];

  for (int k = 0; k < n_constr; ++k) {
    double scale = c_lambda[k] * inv_dtfsq;
    double fx = scale * c_rx[k];
    double fy = scale * c_ry[k];
    double fz = scale * c_rz[k];
    int a = c_atom1[k];
    int b = c_atom2[k];

    int a_route = atom->map(atom->tag[a]);
    int b_route = atom->map(atom->tag[b]);
    if (a_route < 0) a_route = a;
    if (b_route < 0) b_route = b;

    rbuf[a_route][0] += fx; rbuf[a_route][1] += fy; rbuf[a_route][2] += fz;
    rbuf[b_route][0] -= fx; rbuf[b_route][1] -= fy; rbuf[b_route][2] -= fz;

    if (evflag) {
      int atomlist[2];
      int count = 0;
      if (a_route < nlocal) atomlist[count++] = a_route;
      if (b_route < nlocal) atomlist[count++] = b_route;
      bond_v[0] = scale * c_rx[k]*c_rx[k];
      bond_v[1] = scale * c_ry[k]*c_ry[k];
      bond_v[2] = scale * c_rz[k]*c_rz[k];
      bond_v[3] = scale * c_rx[k]*c_ry[k];
      bond_v[4] = scale * c_rx[k]*c_rz[k];
      bond_v[5] = scale * c_ry[k]*c_rz[k];
      double fpair[1] = {scale};
      double dellist[1][3] = {{c_rx[k], c_ry[k], c_rz[k]}};
      int pairlist[1][2] = {{a_route, b_route}};
      v_tally(count, atomlist, 2.0, bond_v, nlocal, 1, pairlist, fpair, dellist);
    }
  }

  // Reverse-comm: ghost contributions on this rank land at the local
  // atoms of the actual owners.  Every rank participates, even those
  // with n_constr == 0 (they may be on the receiving end).
  comm->reverse_comm(this);

  // Accumulate rbuf into atom->f for local atoms.
  for (int i = 0; i < nlocal; ++i) {
    f[i][0] += rbuf[i][0];
    f[i][1] += rbuf[i][1];
    f[i][2] += rbuf[i][2];
    if (store_flag) {
      fstore[i][0] += rbuf[i][0];
      fstore[i][1] += rbuf[i][1];
      fstore[i][2] += rbuf[i][2];
    }
  }

  (void) vflag;
}

/* ----------------------------------------------------------------------
   RATTLE-style velocity correction for cluster-owner solves.

   Same pattern as apply_constraint_forces: write velocity deltas to
   rbuf[a_route]/rbuf[b_route], reverse_comm to deliver ghost
   contributions to their owners' local indices, then add into atom->v.
------------------------------------------------------------------------- */

void FixIlvesLocal::correct_velocities()
{
  grow_lu_workspace(largest_cluster);

  // Zero rbuf upfront so we can accumulate per-cluster updates and ranks
  // with no clusters still send/receive (the collective must match).
  for (int i = 0; i < atom->nmax; ++i)
    rbuf[i][0] = rbuf[i][1] = rbuf[i][2] = 0.0;

  for (int c = 0; c < n_clusters; ++c) {
    const int beg = cluster_offset[c];
    const int end = cluster_offset[c + 1];
    const int n_c = end - beg;

    // assemble symmetric A and rhs g  (same as FixIlves::correct_velocities)
    for (int i = 0; i < n_c * n_c; ++i) lu_A[i] = 0.0;
    for (int s = 0; s < n_c; ++s) {
      int k = c_perm[beg + s];
      int a = c_atom1[k];
      int b = c_atom2[k];
      lu_A[s*n_c + s] = (c_invma[k] + c_invmb[k]) * c_rsq[k];
      double vxd = v[a][0] - v[b][0];
      double vyd = v[a][1] - v[b][1];
      double vzd = v[a][2] - v[b][2];
      lu_b[s] = vxd*c_rx[k] + vyd*c_ry[k] + vzd*c_rz[k];
    }
    for (int s = 0; s < n_c; ++s) {
      int k = c_perm[beg + s];
      double invma_k = c_invma[k];
      double invmb_k = c_invmb[k];
      tagint tag_ak = atom->tag[c_atom1[k]];
      tagint tag_bk = atom->tag[c_atom2[k]];
      for (int t = s + 1; t < n_c; ++t) {
        int l = c_perm[beg + t];
        tagint tag_al = atom->tag[c_atom1[l]];
        tagint tag_bl = atom->tag[c_atom2[l]];
        int sig_k = 0, sig_l = 0;
        double invm_p = 0.0;
        if      (tag_ak == tag_al) { sig_k = +1; sig_l = +1; invm_p = invma_k; }
        else if (tag_ak == tag_bl) { sig_k = +1; sig_l = -1; invm_p = invma_k; }
        else if (tag_bk == tag_al) { sig_k = -1; sig_l = +1; invm_p = invmb_k; }
        else if (tag_bk == tag_bl) { sig_k = -1; sig_l = -1; invm_p = invmb_k; }
        else continue;
        double rkrl = c_rx[k]*c_rx[l] + c_ry[k]*c_ry[l] + c_rz[k]*c_rz[l];
        double val = sig_k * sig_l * rkrl * invm_p;
        lu_A[s*n_c + t] = val;
        lu_A[t*n_c + s] = val;
      }
    }

    ++chol_calls;
    int info = chol_factor_solve(n_c);
    if (info) {
      ++chol_fallbacks;
      // re-assemble and fall back to LU (same as base class)
      for (int i = 0; i < n_c * n_c; ++i) lu_A[i] = 0.0;
      for (int s = 0; s < n_c; ++s) {
        int k = c_perm[beg + s];
        int a = c_atom1[k];
        int b = c_atom2[k];
        lu_A[s*n_c + s] = (c_invma[k] + c_invmb[k]) * c_rsq[k];
        double vxd = v[a][0] - v[b][0];
        double vyd = v[a][1] - v[b][1];
        double vzd = v[a][2] - v[b][2];
        lu_b[s] = vxd*c_rx[k] + vyd*c_ry[k] + vzd*c_rz[k];
      }
      for (int s = 0; s < n_c; ++s) {
        int k = c_perm[beg + s];
        double invma_k = c_invma[k];
        double invmb_k = c_invmb[k];
        tagint tag_ak = atom->tag[c_atom1[k]];
        tagint tag_bk = atom->tag[c_atom2[k]];
        for (int t = s + 1; t < n_c; ++t) {
          int l = c_perm[beg + t];
          tagint tag_al = atom->tag[c_atom1[l]];
          tagint tag_bl = atom->tag[c_atom2[l]];
          int sig_k = 0, sig_l = 0;
          double invm_p = 0.0;
          if      (tag_ak == tag_al) { sig_k = +1; sig_l = +1; invm_p = invma_k; }
          else if (tag_ak == tag_bl) { sig_k = +1; sig_l = -1; invm_p = invma_k; }
          else if (tag_bk == tag_al) { sig_k = -1; sig_l = +1; invm_p = invmb_k; }
          else if (tag_bk == tag_bl) { sig_k = -1; sig_l = -1; invm_p = invmb_k; }
          else continue;
          double rkrl = c_rx[k]*c_rx[l] + c_ry[k]*c_ry[l] + c_rz[k]*c_rz[l];
          double val = sig_k * sig_l * rkrl * invm_p;
          lu_A[s*n_c + t] = val;
          lu_A[t*n_c + s] = val;
        }
      }
      info = lu_factor_solve(n_c);
    }
    if (info)
      error->one(FLERR, "Fix ilves/local: singular velocity-correction matrix in cluster {} "
                 "(size {}).  Check for degenerate constraint topology.", c, n_c);

    // accumulate velocity deltas into rbuf at routed indices
    for (int s = 0; s < n_c; ++s) {
      int k = c_perm[beg + s];
      double mu = lu_b[s];
      int a = c_atom1[k];
      int b = c_atom2[k];
      int a_route = atom->map(atom->tag[a]);
      int b_route = atom->map(atom->tag[b]);
      if (a_route < 0) a_route = a;
      if (b_route < 0) b_route = b;
      double da = mu * c_invma[k];
      double db = mu * c_invmb[k];
      rbuf[a_route][0] -= da * c_rx[k];
      rbuf[a_route][1] -= da * c_ry[k];
      rbuf[a_route][2] -= da * c_rz[k];
      rbuf[b_route][0] += db * c_rx[k];
      rbuf[b_route][1] += db * c_ry[k];
      rbuf[b_route][2] += db * c_rz[k];
    }
  }

  // Reverse-comm velocity deltas to owners.
  comm->reverse_comm(this);

  // Apply deltas to atom->v for local atoms.
  for (int i = 0; i < nlocal; ++i) {
    v[i][0] += rbuf[i][0];
    v[i][1] += rbuf[i][1];
    v[i][2] += rbuf[i][2];
  }
}

/* ----------------------------------------------------------------------
   Project initial atom positions onto the constraint surface at setup.
   Mirrors FixIlves::correct_coordinates but uses the local variant's
   apply_constraint_forces (so ghost-atom contributions get reverse_comm'd).
------------------------------------------------------------------------- */

void FixIlvesLocal::correct_coordinates(int vflag)
{
  if (n_constr == 0) {
    // still participate in any reverse_comm a neighbor may want to do
    // by routing through apply_constraint_forces below; n_constr==0 just
    // means we don't write anything to rbuf.
  }

  double **f_save = nullptr;
  double **v_save = nullptr;
  memory->create(f_save, nlocal, 3, "ilves/local:fsave");
  memory->create(v_save, nlocal, 3, "ilves/local:vsave");
  for (int i = 0; i < nlocal; ++i) {
    f_save[i][0] = f[i][0]; f_save[i][1] = f[i][1]; f_save[i][2] = f[i][2];
    v_save[i][0] = v[i][0]; v_save[i][1] = v[i][1]; v_save[i][2] = v[i][2];
    f[i][0] = f[i][1] = f[i][2] = 0.0;
    v[i][0] = v[i][1] = v[i][2] = 0.0;
  }

  unconstrained_update();
  comm->forward_comm(this);
  solve_constraints();
  apply_constraint_forces(vflag);

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

  for (int i = 0; i < nlocal; ++i) {
    f[i][0] = f_save[i][0]; f[i][1] = f_save[i][1]; f[i][2] = f_save[i][2];
    v[i][0] = v_save[i][0]; v[i][1] = v_save[i][1]; v[i][2] = v_save[i][2];
  }
  memory->destroy(f_save);
  memory->destroy(v_save);

  // forward-comm updated x via xshake pack
  double **xtmp = xshake;
  xshake = x;
  comm->forward_comm(this);
  xshake = xtmp;

  precompute_constraint_data();
}

/* ----------------------------------------------------------------------
   Compute angle_distance[at] = r_AC for each selected angle type at, using
   the law of cosines on the flanking bond equilibrium distances and the
   equilibrium angle.  The pair of flanking bond types for each angle type
   is discovered from local angle storage; ranks may not each see every
   angle type, so the per-rank result is reduced via MPI_MAX so all ranks
   agree.  We allreduce both b1 and b2 separately (taking MAX) and detect
   cross-rank disagreement (the same angle type appearing with different
   flanking bond pairs).
------------------------------------------------------------------------- */

void FixIlvesLocal::compute_angle_distances_local()
{
  if (!has_angle) return;

  const int natypes  = atom->nangletypes;
  int **na_type      = atom->angle_type;
  tagint **na_atom1  = atom->angle_atom1;
  tagint **na_atom2  = atom->angle_atom2;
  tagint **na_atom3  = atom->angle_atom3;
  int *num_angle     = atom->num_angle;
  int **nb_type      = atom->bond_type;
  tagint **nb_atom   = atom->bond_atom;
  int *num_bond      = atom->num_bond;
  tagint *tag        = atom->tag;
  const int nlocal_now = atom->nlocal;

  // local accumulators: for each angle type, the (min, max) pair of
  // flanking bond types as observed in local angle storage.  Initialize
  // to 0 (no observation yet).
  std::vector<int> b_local(natypes + 1, 0);
  std::vector<int> B_local(natypes + 1, 0);
  std::vector<int> conflict_local(natypes + 1, 0);

  auto local_bond_type_between = [&](int ia, tagint tb) -> int {
    if (ia < 0 || ia >= nlocal_now) return 0;
    for (int b = 0; b < num_bond[ia]; ++b) {
      int bt = nb_type[ia][b];
      if (bt == 0) continue;
      if (bt < 0) bt = -bt;
      if (nb_atom[ia][b] == tb) return bt;
    }
    return 0;
  };

  for (int i = 0; i < nlocal_now; ++i) {
    for (int m = 0; m < num_angle[i]; ++m) {
      int at = na_type[i][m];
      if (at == 0) continue;
      if (at < 0) at = -at;
      if (at > natypes || !angle_flag[at]) continue;
      // Use only middle-atom-owned angles to avoid newton-off triple
      // counting messing with the flanking-bond lookup.
      if (na_atom2[i][m] != tag[i]) continue;
      tagint t1 = na_atom1[i][m];
      tagint t3 = na_atom3[i][m];
      // flanking bond types via the middle atom's local bond list
      int bt1 = local_bond_type_between(i, t1);
      int bt2 = local_bond_type_between(i, t3);
      if (bt1 == 0 || bt2 == 0) continue;
      int bmin = MIN(bt1, bt2);
      int bmax = MAX(bt1, bt2);
      if (b_local[at] == 0) { b_local[at] = bmin; B_local[at] = bmax; }
      else if (b_local[at] != bmin || B_local[at] != bmax) conflict_local[at] = 1;
    }
  }

  // share across ranks: take MAX of (bmin, bmax) per angle type.  This is
  // consistent because all ranks that observe a given angle type must see
  // the same flanking pair (a hard cross-rank constraint of the input).
  std::vector<int> b_all(natypes + 1, 0);
  std::vector<int> B_all(natypes + 1, 0);
  std::vector<int> conflict_all(natypes + 1, 0);
  MPI_Allreduce(b_local.data(),        b_all.data(),        natypes + 1, MPI_INT, MPI_MAX, world);
  MPI_Allreduce(B_local.data(),        B_all.data(),        natypes + 1, MPI_INT, MPI_MAX, world);
  MPI_Allreduce(conflict_local.data(), conflict_all.data(), natypes + 1, MPI_INT, MPI_MAX, world);

  // detect cross-rank disagreement: if any rank reported a non-zero pair
  // that differs from the consensus MAX, flag conflict.  (A rank that
  // saw 0 entries doesn't contribute.)
  std::vector<int> b_other(natypes + 1, 0);
  std::vector<int> B_other(natypes + 1, 0);
  for (int at = 1; at <= natypes; ++at) {
    // pack -min instead of min so MAX reduction returns true min; same for max.
    b_other[at] = (b_local[at] > 0) ? -b_local[at] : -INT_MAX;
    B_other[at] = (B_local[at] > 0) ? -B_local[at] : -INT_MAX;
  }
  std::vector<int> b_min_all(natypes + 1, 0);
  std::vector<int> B_min_all(natypes + 1, 0);
  MPI_Allreduce(b_other.data(), b_min_all.data(), natypes + 1, MPI_INT, MPI_MAX, world);
  MPI_Allreduce(B_other.data(), B_min_all.data(), natypes + 1, MPI_INT, MPI_MAX, world);

  for (int at = 1; at <= natypes; ++at) {
    if (!angle_flag[at]) { angle_distance[at] = 0.0; continue; }
    int bmin_max = b_all[at];
    int bmax_max = B_all[at];
    int bmin_min = -b_min_all[at];
    int bmax_min = -B_min_all[at];
    if (bmin_max == 0) {
      // no rank observed this angle type with selectable flanking bonds
      angle_distance[at] = 0.0;
      continue;
    }
    if (bmin_max != bmin_min || bmax_max != bmax_min || conflict_all[at])
      error->all(FLERR, "Fix ilves/local: angle type {} spans bonds of mixed types", at);
    const double theta0 = force->angle->equilibrium_angle(at);
    const double r1 = bond_distance[bmin_max];
    const double r2 = bond_distance[bmax_max];
    angle_r1[at] = r1;
    angle_r2[at] = r2;
    angle_distance[at] = sqrt(r1*r1 + r2*r2 - 2.0*r1*r2*cos(theta0));
  }
}

/* ----------------------------------------------------------------------
   init_topology: validate newton_bond setting, scan local bond/angle
   storage to validate every cluster fits in local+ghost reach and is
   star-shaped (or 1+1), and compute angle_distance[].
------------------------------------------------------------------------- */

void FixIlvesLocal::init_topology()
{
  if (force->newton_bond)
    error->all(FLERR,
               "Fix ilves/local requires newton off bond.  Either set "
               "'newton off' (or 'newton on off') before fix ilves/local, "
               "or switch to fix ilves/global.");

  // angle_distance[] depends only on bond_distance[] (already filled by
  // FixIlves::init()) and on local angle storage (no ghosts needed), so
  // we can compute it now.  Cluster-reachability and star-topology
  // validation, by contrast, need ghost atoms and therefore deferred
  // until the first post_neighbor / build_constraint_list call.
  compute_angle_distances_local();
}

/* ----------------------------------------------------------------------
   Decide whether the (ta, tb) bond is selected for constraint, using
   only local bond storage.  Required by negate_constrained_topology()
   in the FixIlves base class to detect which angles to negate (both
   flanking bonds must be themselves constrained).

   Both endpoints must be reachable locally for this variant; the init
   validation has already enforced that.  We walk whichever endpoint is
   local to find the bond type.
------------------------------------------------------------------------- */

bool FixIlvesLocal::bond_is_constrained(tagint ta, tagint tb)
{
  int ia = atom->map(ta);
  int ib = atom->map(tb);
  if (ia < 0 || ib < 0) return false;

  const int nlocal_now = atom->nlocal;
  int **nb_type    = atom->bond_type;
  tagint **nb_atom = atom->bond_atom;
  int *num_bond    = atom->num_bond;

  auto find_bt = [&](int local_idx, tagint partner_tag) -> int {
    if (local_idx < 0 || local_idx >= nlocal_now) return 0;
    for (int b = 0; b < num_bond[local_idx]; ++b) {
      int bt = nb_type[local_idx][b];
      if (bt == 0) continue;
      if (bt < 0) bt = -bt;
      if (nb_atom[local_idx][b] == partner_tag) return bt;
    }
    return 0;
  };

  // try whichever side is locally owned
  int bt = find_bt(ia, tb);
  if (bt == 0) bt = find_bt(ib, ta);
  if (bt == 0) return false;
  if (bt > atom->nbondtypes) return false;

  return bond_selected_for_atoms(ia, ib, bt);
}

/* ----------------------------------------------------------------------
   build_constraint_list: scan local atoms, identify cluster centers
   (atoms with > 1 constrained bond) and 1+1 pairs, and add their
   constraints to the flat per-rank list.

   Cluster ownership:
     * Star clusters (count > 1):  this rank owns iff center is local.
     * 1+1 clusters (count == 1 on both endpoints): owned by the rank
       holding the lower-tag atom; we add the constraint only when the
       lower-tag atom is local AND this iteration's i is that atom.
------------------------------------------------------------------------- */

void FixIlvesLocal::build_constraint_list()
{
  const int nlocal_now = atom->nlocal;
  int **nb_type      = atom->bond_type;
  tagint **nb_atom   = atom->bond_atom;
  int *num_bond      = atom->num_bond;
  int **na_type      = atom->angle_type;
  tagint **na_atom1  = atom->angle_atom1;
  tagint **na_atom2  = atom->angle_atom2;
  tagint **na_atom3  = atom->angle_atom3;
  int *num_angle     = atom->num_angle;
  int *mask          = atom->mask;
  tagint *tag        = atom->tag;

  // mark per-atom flags as we go (zero first)
  for (int i = 0; i < atom->nmax; ++i) ilves_flag[i] = 0;
  n_constr = 0;

  // Validate cluster reachability and star-topology for every local
  // atom.  This is cheap (linear in num_bond per atom) and must run
  // every post_neighbor because atoms migrate.  We detect:
  //   * a constraint partner that is neither local nor ghost (cluster
  //     too big for this rank's communication shell), and
  //   * adjacent branching points (multi-center cluster) that this
  //     variant cannot handle correctly.
  for (int i = 0; i < nlocal_now; ++i) {
    if (!(mask[i] & groupbit)) continue;
    int my_count = 0;
    int partner_local_with_multibond = 0;
    for (int b = 0; b < num_bond[i]; ++b) {
      int bt = nb_type[i][b];
      if (bt == 0) continue;
      if (bt < 0) bt = -bt;
      if (bt > atom->nbondtypes) continue;
      tagint pt = nb_atom[i][b];
      int j = atom->map(pt);
      if (j < 0) {
        // bond partner not reachable: only an error if this bond type is
        // explicitly selected for constraint (bond_flag[bt]).  Type/mass
        // selectors are silently skipped because evaluating them needs
        // the partner's atom data, which we don't have without a ghost.
        if (bond_flag[bt])
          error->one(FLERR,
                     "Fix ilves/local: cluster atom (tag {}) not reachable from rank {} "
                     "as either a local atom or a ghost.  The constraint cluster spans "
                     "beyond this rank's local subdomain plus the communication cutoff.  "
                     "Switch to fix ilves/global, or extend the communication cutoff via "
                     "`comm_modify cutoff <value>` so all clusters fit.", pt, comm->me);
        continue;
      }
      if (!(mask[j] & groupbit)) continue;
      if (!bond_selected_for_atoms(i, j, bt)) continue;
      ++my_count;
      if (j < nlocal_now) {
        int pc = count_constrained_bonds_local(j);
        if (pc > 1) ++partner_local_with_multibond;
      }
    }
    if (my_count > 1 && partner_local_with_multibond > 0) {
      error->one(FLERR,
                 "Fix ilves/local: atom (tag {}) is part of a multi-center "
                 "constraint cluster (non-star topology, e.g. an interior atom of a "
                 "constrained 3+-bond chain).  Use fix ilves/global for clusters "
                 "with more than one branching point.", tag[i]);
    }
  }

  for (int i = 0; i < nlocal_now; ++i) {
    if (!(mask[i] & groupbit)) continue;

    // collect i's constrained bonds (partner index + bond type)
    struct LocBond { int p; int bt; };
    std::vector<LocBond> bonds_i;
    bonds_i.reserve(num_bond[i]);
    for (int b = 0; b < num_bond[i]; ++b) {
      int bt = nb_type[i][b];
      if (bt == 0) continue;
      if (bt < 0) bt = -bt;
      if (bt > atom->nbondtypes) continue;
      int j = atom->map(nb_atom[i][b]);
      if (j < 0) continue;       // already errored above
      if (!(mask[j] & groupbit)) continue;
      if (!bond_selected_for_atoms(i, j, bt)) continue;
      bonds_i.push_back({j, bt});
    }

    const int kc = (int) bonds_i.size();
    if (kc == 0) continue;

    if (kc > 1) {
      // i is a star center -- this rank owns the cluster
      for (auto &lb : bonds_i) {
        int j = domain->closest_image(i, lb.p);
        add_constraint(i, j, lb.bt, bond_distance[lb.bt]);
        ilves_flag[i] = 1;
        ilves_flag[j] = 1;
      }
      // angle (A-C virtual) constraints centered at i
      if (has_angle) {
        for (int m = 0; m < num_angle[i]; ++m) {
          int at = na_type[i][m];
          if (at == 0) continue;
          if (at < 0) at = -at;
          if (at > atom->nangletypes) continue;
          if (!angle_flag[at]) continue;
          // require i to be the middle atom
          if (na_atom2[i][m] != tag[i]) continue;
          tagint t1 = na_atom1[i][m];
          tagint t3 = na_atom3[i][m];
          int i1 = atom->map(t1);
          int i3 = atom->map(t3);
          if (i1 < 0 || i3 < 0) continue;
          if (!(mask[i1] & groupbit) || !(mask[i3] & groupbit)) continue;
          if (!bond_is_constrained(tag[i], t1)) continue;
          if (!bond_is_constrained(tag[i], t3)) continue;
          int a_idx = domain->closest_image(i, i1);
          int b_idx = domain->closest_image(i, i3);
          // canonical ordering by tag so the pair matches whatever the
          // star bond constraints already inserted: c_atom1 is lower-tag.
          int a_out, b_out;
          if (tag[a_idx] < tag[b_idx]) { a_out = a_idx; b_out = b_idx; }
          else                         { a_out = b_idx; b_out = a_idx; }
          // place at i if tag[i] is lowest (it isn't -- a/b are leaves of
          // the cluster centered at i, but for the virtual A-C constraint
          // the ownership atom must still be a local one).  For star
          // clusters owned by i, all three atoms are reachable; we put
          // c_atom1 = whichever of (a_out, b_out) is local (preferring
          // the lower-tag one).  This keeps the base class's local-atom
          // routing for force application correct.
          int a_local = (a_out < nlocal_now) ? a_out : -1;
          int b_local = (b_out < nlocal_now) ? b_out : -1;
          if (a_local < 0 && b_local < 0) {
            // both leaves are ghost -- this rank is the cluster owner so
            // it must still process the constraint.  c_atom1 = first
            // ghost; apply_constraint_forces will route via atom->map().
            add_constraint(a_out, b_out, -at, angle_distance[at]);
          } else if (a_local >= 0) {
            add_constraint(a_local, b_out, -at, angle_distance[at]);
          } else {
            add_constraint(b_local, a_out, -at, angle_distance[at]);
          }
          ilves_flag[a_out] = 1;
          ilves_flag[b_out] = 1;
          ilves_flag[i]     = 1;
        }
      }
    } else {
      // kc == 1: leaf of star or part of 1+1 cluster
      int j = bonds_i[0].p;
      int bt = bonds_i[0].bt;
      // determine partner's degree to distinguish
      // (a) star with remote center : partner is ghost or partner is local
      //                                with count > 1 -- skip in both cases.
      // (b) 1+1 cluster : partner is local with count == 1 -- claim if
      //                   tag[i] < tag[j].
      if (j >= nlocal_now) continue;        // partner is ghost: center is remote, skip
      int pc = count_constrained_bonds_local(j);
      if (pc != 1) continue;                // partner is local star center, partner owns
      // both are local 1+1 leaves: owner = lower-tag atom
      if (tag[i] >= tag[j]) continue;
      int b_idx = domain->closest_image(i, j);
      add_constraint(i, b_idx, bt, bond_distance[bt]);
      ilves_flag[i]     = 1;
      ilves_flag[b_idx] = 1;
    }
  }

  // connected-component labelling and cluster grouping (base class)
  group_by_cluster();
  // precompute reference vectors, |r|^2, and inverse masses (base class)
  precompute_constraint_data();
}
