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
   Global-topology variant of the ILVES bond/angle constraint solver.
   See fix_ilves.cpp for the algorithm description and the FixIlves base
   class.  This file supplies the constraint-discovery half of the fix:
   it gathers the full bond/angle table onto every rank at init() and
   builds the per-rank constraint list from that replicated table.
------------------------------------------------------------------------- */

#include "fix_ilves_global.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"

#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <vector>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixIlvesGlobal::FixIlvesGlobal(LAMMPS *lmp, int narg, char **arg) :
    FixIlves(lmp, narg, arg), global_topology_ready(false)
{
}

/* ----------------------------------------------------------------------
   init_topology: gather the global bond/angle table, build the
   tag-to-cluster map, and fill angle_distance[] for selected angle types
   using the replicated angle table.  Runs once from FixIlves::init().
------------------------------------------------------------------------- */

void FixIlvesGlobal::init_topology()
{
  gather_global_topology();

  // angle distances: for each constrained angle type, scan the global
  // angle table for the bond types of its flanking bonds.  Computing this
  // from the global table avoids problems where a given rank's local
  // angle is of a different angle type or where the relevant bond is
  // stored on a remote rank with newton_bond=on.
  if (has_angle) {
    for (int at = 1; at <= atom->nangletypes; ++at) {
      if (!angle_flag[at]) { angle_distance[at] = 0.0; continue; }

      int b1 = 0, b2 = 0;
      int conflict = 0;
      const int na = (int) ga1.size();
      for (int gi = 0; gi < na; ++gi) {
        if (ga_type[gi] != at) continue;
        tagint tA = ga1[gi];
        tagint tB = ga2[gi];     // middle atom
        tagint tC = ga3[gi];
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
        int bt1 = lookup_bt(tA, tB);
        int bt2 = lookup_bt(tB, tC);
        if (bt1 == 0 || bt2 == 0) continue;
        int bmin = MIN(bt1, bt2);
        int bmax = MAX(bt1, bt2);
        if (b1 == 0) { b1 = bmin; b2 = bmax; }
        else if (bmin != b1 || bmax != b2) { conflict = 1; break; }
      }
      if (conflict)
        error->all(FLERR, "Fix ilves/global: angle type {} spans bonds of mixed types", at);
      if (b1 == 0) { angle_distance[at] = 0.0; continue; }
      const double theta0 = force->angle->equilibrium_angle(at);
      const double r1 = bond_distance[b1];
      const double r2 = bond_distance[b2];
      angle_r1[at] = r1;
      angle_r2[at] = r2;
      angle_distance[at] = sqrt(r1*r1 + r2*r2 - 2.0*r1*r2*cos(theta0));
    }
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

void FixIlvesGlobal::gather_global_topology()
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
  sendb.reserve(3 * nlocal_now);     // approximate
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

  // -----------------------------------------------------------------
  // dedup bonds: sort by (min_tag, max_tag), keep unique
  // -----------------------------------------------------------------
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
  // unique by (a, b) -- the type is consistent across ranks for a given bond
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
  // locally-stored angle.  atom2 = middle atom.
  //
  // Pre-filter by angle_flag[at]: every consumer of ga_ (the angle
  // distance computation in init_topology and the virtual A-C constraint
  // scan in build_constraint_list) immediately discards entries with a
  // non-selected type, so gathering them is wasted MPI bandwidth and
  // memory.  Skipping them at the source keeps ga_ down to only the
  // angle types we will actually constrain.
  // -----------------------------------------------------------------
  std::vector<tagint> senda;
  for (int i = 0; i < nlocal_now; ++i) {
    for (int m = 0; m < num_angle[i]; ++m) {
      int at = na_type[i][m];
      if (at == 0) continue;
      if (at < 0) at = -at;
      if (at > atom->nangletypes || !angle_flag[at]) continue;
      // dedupe: only middle atom owner packs the angle entry to avoid
      // duplicates from newton-off triple storage.
      if (na_atom2[i][m] != tag[i]) continue;
      // canonical order: outer atoms by tag ordering so dedup works.
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
  // Build a global cluster-id table.  Run union-find over all tags
  // involved in any bond or angle (linked via bonds: a-b in a bond, and
  // a-b, b-c in an angle).  After this, tag_cluster[t] gives the
  // representative tag of the cluster that tag t belongs to.
  //
  // Sparse: only tags that participate in some bond/angle are inserted
  // into the map.  Membership tests on uninvolved tags get the iterator
  // end() and short-circuit.
  // -----------------------------------------------------------------
  tag_cluster.clear();
  tag_cluster.reserve(gb_a.size() * 2 + ga1.size() * 3);

  // initialize: every involved tag is its own parent
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
      it->second = pit->second;     // path compression
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
  // flatten: tag_cluster[t] is now the cluster's representative tag
  for (auto &kv : tag_cluster) kv.second = find(kv.first);

  if (comm->me == 0)
    utils::logmesg(lmp, "Fix ilves/global: gathered global topology with {} bonds and {} selected angles\n",
                   gb_a.size(), ga1.size());

  global_topology_ready = true;
}

/* ----------------------------------------------------------------------
   Helper: given two tags, return true if the (ta, tb) bond is selected
   for constraint.  Looks up the bond type via the global topology table.
   Both atoms must be either local or available as ghosts on this rank
   for the type/mass selectors to work.  Returns false otherwise.
------------------------------------------------------------------------- */

bool FixIlvesGlobal::bond_is_constrained(tagint ta, tagint tb)
{
  tagint tmin = MIN(ta, tb);
  tagint tmax = MAX(ta, tb);
  // binary search in gb_a/gb_b (sorted by (a, b))
  int lo = 0, hi = (int) gb_a.size();
  while (lo < hi) {
    int mid = (lo + hi) / 2;
    if (gb_a[mid] < tmin || (gb_a[mid] == tmin && gb_b[mid] < tmax)) lo = mid + 1;
    else hi = mid;
  }
  if (lo == (int) gb_a.size() || gb_a[lo] != tmin || gb_b[lo] != tmax) return false;

  int ia = atom->map(ta);
  int ib = atom->map(tb);
  if (ia < 0 || ib < 0) return false;     // can't check selectors without atoms
  return bond_selected_for_atoms(ia, ib, gb_type[lo]);
}

/* ----------------------------------------------------------------------
   build the flat constraint list from the replicated bond/angle tables.
   See FixIlves base class for the shared cluster grouping and precompute
   stages that this method calls at the end.

   - bond constraints: pair (i, partner) for every bond where:
       both atoms are in fix group
       AND bond type is in bond_flag, OR
           either atom type is in type_flag, OR
           either atom mass matches mass_list
   - angle "virtual" A-C constraints: derived from angle entries where
       angle type is in angle_flag
       AND both flanking bonds (A-B, B-C) are also selected (per above)
       AND all three atoms are in the fix group
------------------------------------------------------------------------- */

void FixIlvesGlobal::build_constraint_list()
{
  const int nlocal_now = atom->nlocal;
  int *mask = atom->mask;
  tagint *tag = atom->tag;

  // mark per-atom flags as we go (zero first)
  for (int i = 0; i < atom->nmax; ++i) ilves_flag[i] = 0;
  n_constr = 0;

  // Lazy gather of global topology on first build (in case init_topology
  // was bypassed for some reason; normally already done at init()).
  if (!global_topology_ready) gather_global_topology();

  // -----------------------------------------------------------------
  // Identify clusters whose atoms intersect my local atoms.  For each
  // cluster representative tag in this set, I will include EVERY
  // constraint in that cluster, even constraints between two ghost
  // atoms.  This is essential for correct cluster-local solves when
  // a cluster spans multiple ranks (e.g. a water HOH where O is on one
  // rank and H1, H2 on others -- every participating rank needs the
  // full 3-constraint cluster to compute the same Lagrange multipliers).
  // -----------------------------------------------------------------
  // Set of cluster representative tags whose clusters intersect my local
  // atoms.  Sparse (only involved tags participate in any cluster).
  std::unordered_map<tagint, char> my_cluster;
  my_cluster.reserve(nlocal_now);
  for (int i = 0; i < nlocal_now; ++i) {
    tagint t = tag[i];
    auto it = tag_cluster.find(t);
    if (it != tag_cluster.end()) my_cluster[it->second] = 1;
  }

  // helper to test whether a tag is in any of my clusters
  auto in_my_cluster = [&](tagint t) -> bool {
    auto it = tag_cluster.find(t);
    if (it == tag_cluster.end()) return false;
    return my_cluster.count(it->second) > 0;
  };

  // helper: pick c_atom1 and c_atom2 indices given that at least one of
  // (ta, tb) is in my cluster.  Returns false if neither atom is even
  // available as a ghost on this rank (constraint cannot be processed).
  // The chosen c_atom1 is always a local atom if any local atom is in
  // the pair; otherwise it's a ghost.  c_atom2 is the partner, with
  // PBC closest-image relative to c_atom1.
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
      // neither local -- pick whichever has lower tag as c_atom1.  This
      // is a "ghost-only" constraint that this rank includes solely for
      // cluster completeness.  Forces will not be applied to either
      // atom by this rank (apply_constraint_forces tests < nlocal).
      if (tag[ja] < tag[jb]) { a_out = ja; b_out = jb; }
      else                   { a_out = jb; b_out = ja; }
    }
    b_out = domain->closest_image(a_out, b_out);
    return true;
  };

  // -----------------------------------------------------------------
  // Pre-phase: for near-linear angle types, identify the bond legs we
  // will SKIP in Phase A (the one between B and the higher-tag endpoint
  // of {A,C}).  We do this here, before the bond loop, so Phase A can
  // cheaply test each candidate bond against the drop set.
  //
  // Stored as (lo_tag, hi_tag) keys with lo < hi to match the canonical
  // ordering of gb_a / gb_b.  Only bonds whose flanking angle has
  // angle_linear[at]==1 AND whose other flanking bond is constrained
  // get a slot here -- mirroring the conditions Phase B uses below.
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

    // Drop the redundant leg of any near-linear angle that touches both
    // ends of this bond.  The B-M virtual constraint (emitted in Phase
    // B) plus the retained A-B leg and A-C virtual together determine
    // |B-C| geometrically, so adding the B-C bond again would make the
    // Jacobian rank-deficient (the very condition we're trying to fix).
    if (!dropped_bonds.empty() && dropped_bonds.count({ta, tb})) continue;

    int a_idx, b_idx;
    if (!pick_atoms(ta, tb, a_idx, b_idx)) continue;

    // group + selector checks
    if (!(mask[a_idx] & groupbit) || !(mask[b_idx] & groupbit)) continue;
    if (!bond_selected_for_atoms(a_idx, b_idx, bt)) continue;

    add_constraint(a_idx, b_idx, bt, bond_distance[bt]);
    ilves_flag[a_idx] = 1;
    ilves_flag[b_idx] = 1;
  }

  // -----------------------------------------------------------------
  // Phase B: angle (A-C "virtual") constraints from clusters I'm in.
  // For near-linear angle types we additionally emit a 3-atom B-M
  // virtual constraint (handled by the alternate add_constraint signature)
  // so that the Jacobian stays well-conditioned at theta near 180 deg.
  // -----------------------------------------------------------------
  if (has_angle) {
    const int na_global = (int) ga1.size();
    for (int gi = 0; gi < na_global; ++gi) {
      int at = ga_type[gi];
      if (at <= 0 || at > atom->nangletypes) continue;
      if (!angle_flag[at]) continue;

      tagint t1 = ga1[gi];     // outer atoms canonical (sorted)
      tagint t2 = ga2[gi];     // middle
      tagint t3 = ga3[gi];

      if (!in_my_cluster(t1) && !in_my_cluster(t2) && !in_my_cluster(t3)) continue;

      int i1 = atom->map(t1);
      int i2 = atom->map(t2);
      int i3 = atom->map(t3);
      if (i1 < 0 || i2 < 0 || i3 < 0) continue;

      // require all three atoms in the fix group
      if (!(mask[i1] & groupbit)) continue;
      if (!(mask[i2] & groupbit)) continue;
      if (!(mask[i3] & groupbit)) continue;

      // both flanking bonds (i2-i1 and i2-i3) must be themselves selected
      if (!bond_is_constrained(t2, t1)) continue;
      if (!bond_is_constrained(t2, t3)) continue;

      int a_idx, b_idx;
      if (!pick_atoms(t1, t3, a_idx, b_idx)) continue;

      add_constraint(a_idx, b_idx, -at, angle_distance[at]);
      ilves_flag[a_idx] = 1;
      ilves_flag[b_idx] = 1;
      ilves_flag[i2]    = 1;

      if (angle_linear[at]) {
        // 3-atom B-M virtual constraint.  c_atom1 = B (central),
        // c_atom2 = A (lower tag), c_atom3 = C (higher tag).  Use
        // closest_image so A and C are PBC-consistent with B's copy.
        int iB = i2;
        int iA = atom->map(t1);
        int iC = atom->map(t3);
        // Route ownership to a local copy of B if available so the
        // constraint's c_atom1 is local (consistent with the 2-atom
        // ownership rule); if B is ghost-only, fall back to its ghost.
        // This rank still owns the constraint matrix entry for cluster
        // completion in the redundant-solve global variant.
        int iB_local = atom->map(t2);
        if (iB_local >= 0) iB = iB_local;
        // Canonicalize ghost images of A and C relative to B's copy so
        // x[A], x[C] are in the same periodic image as x[B].
        iA = domain->closest_image(iB, iA);
        iC = domain->closest_image(iB, iC);
        add_constraint(iB, iA, iC, -at, angle_dBM[at]);
        ilves_flag[iB] = 1;
        ilves_flag[iA] = 1;
        ilves_flag[iC] = 1;
      }
    }
  }

  // -----------------------------------------------------------------
  // connected-component labelling and cluster grouping (base class)
  // -----------------------------------------------------------------
  group_by_cluster();

  // -----------------------------------------------------------------
  // precompute reference vectors, |r|^2, and inverse masses (base class)
  // -----------------------------------------------------------------
  precompute_constraint_data();
}
