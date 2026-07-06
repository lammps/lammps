/* ----------------------------------------------------------------------
   pair_style template/offset

   Adds a constant energy `offset` to the system for each distinct match of a
   molecule template.  Registered with
       pair_coeff * * <mol-id> <offset>
   (one line per template; templates accumulate).  Applies to all atom types.

   The offset is a per-configuration constant, so the style contributes energy
   only -- forces and per-atom virial are zero.  A single-molecule (connected)
   template is matched with the msevb superimpose kernel; the number of distinct
   (non-overlapping) matches is counted and `offset * count` is added to the
   potential energy.

   Typical use is via hybrid/overlay alongside the real force field.  Because
   fix msevb evaluates every EVB diabatic state through the normal pair pipeline,
   the offset is automatically included in each state's energy -- giving a per-
   species diabatic shift without any bookkeeping inside the fix.
------------------------------------------------------------------------- */

#include "pair_template_offset.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "molecule.h"
#include "utils.h"

#include <algorithm>
#include <cstring>
#include <set>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairTemplateOffset::PairTemplateOffset(LAMMPS *lmp) : Pair(lmp), cut_global(0.0)
{
  single_enable = 0;       // no pairwise energy/force decomposition
  restartinfo = 0;         // nothing to write to restart
  respa_enable = 0;
  gref.type = nullptr;
  gref.nspecial_flat = nullptr;
  gref.special_flat = nullptr;
  gref.maxspecial = 0;
  gref.maxtag = 0;
}

/* ---------------------------------------------------------------------- */

PairTemplateOffset::~PairTemplateOffset()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairTemplateOffset::allocate()
{
  allocated = 1;
  const int n = atom->ntypes;
  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 0;
  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
}

/* ----------------------------------------------------------------------
   pair_style template/offset [cutoff]
   The optional cutoff only sets the (unused) pairwise cutoff; matching does not
   use the neighbor list.  Defaults to 0.0.
------------------------------------------------------------------------- */

void PairTemplateOffset::settings(int narg, char **arg)
{
  if (narg > 1) error->all(FLERR, "Illegal pair_style template/offset command");
  cut_global = (narg == 1) ? utils::numeric(FLERR, arg[0], false, lmp) : 0.0;
}

/* ----------------------------------------------------------------------
   pair_coeff * * <mol-id> <offset>
------------------------------------------------------------------------- */

void PairTemplateOffset::coeff(int narg, char **arg)
{
  if (!allocated) allocate();
  if (narg != 4) error->all(FLERR, "pair_coeff template/offset requires '* * <mol-id> <offset>'");
  if (strcmp(arg[0], "*") != 0 || strcmp(arg[1], "*") != 0)
    error->all(FLERR, "pair_style template/offset requires 'pair_coeff * *'");

  int imol = atom->find_molecule(arg[2]);
  if (imol < 0)
    error->all(FLERR,
               "pair_style template/offset: molecule template '{}' not found; define it with the "
               "molecule command first",
               arg[2]);

  Entry e;
  e.mol_id = arg[2];
  e.mol_index = imol;
  e.offset = utils::numeric(FLERR, arg[3], false, lmp);
  entries.push_back(e);

  const int n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 1;
}

/* ---------------------------------------------------------------------- */

void PairTemplateOffset::init_style()
{
  if (atom->molecular == Atom::ATOMIC)
    error->all(FLERR, "pair_style template/offset requires a molecular atom style");
  // No neighbor list is needed: matching uses the bond (1-2 special) topology.
}

/* ---------------------------------------------------------------------- */

double PairTemplateOffset::init_one(int /*i*/, int /*j*/)
{
  return cut_global;
}

/* ----------------------------------------------------------------------
   Gather the whole-system 1-2 topology (type and 1-2 neighbor tags, indexed by
   global atom tag) so that a template match whose atoms are spread over several
   MPI ranks is still found.  Every rank ends up with the same snapshot.
------------------------------------------------------------------------- */

void PairTemplateOffset::build_global_topology()
{
  const int nlocal = atom->nlocal;
  tagint *tag = atom->tag;

  tagint maxtag = 0;
  for (int i = 0; i < nlocal; i++) maxtag = MAX(maxtag, tag[i]);
  MPI_Allreduce(MPI_IN_PLACE, &maxtag, 1, MPI_LMP_TAGINT, MPI_MAX, world);

  const int mspc = atom->maxspecial;
  const size_t sz = (size_t) maxtag + 1;

  g_type.assign(sz, 0);
  g_nspecial.assign(sz * 3, 0);
  g_special.assign(sz * MAX(mspc, 1), 0);

  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  int *type = atom->type;

  for (int i = 0; i < nlocal; i++) {
    const tagint t = tag[i];
    g_type[t] = type[i];
    for (int k = 0; k < 3; k++) g_nspecial[(size_t) t * 3 + k] = nspecial[i][k];
    for (int k = 0; k < mspc; k++) g_special[(size_t) t * mspc + k] = special[i][k];
  }

  if (comm->nprocs > 1) {
    MPI_Allreduce(MPI_IN_PLACE, g_type.data(), sz, MPI_INT, MPI_SUM, world);
    MPI_Allreduce(MPI_IN_PLACE, g_nspecial.data(), sz * 3, MPI_INT, MPI_SUM, world);
    MPI_Allreduce(MPI_IN_PLACE, g_special.data(), sz * MAX(mspc, 1), MPI_LMP_TAGINT, MPI_SUM, world);
  }

  gref.type = g_type.data();
  gref.nspecial_flat = g_nspecial.data();
  gref.special_flat = g_special.data();
  gref.maxspecial = mspc;
  gref.maxtag = maxtag;
}

/* ----------------------------------------------------------------------
   Count the distinct matches of a connected molecule template in the current
   (global) topology.  A match is a set of real atoms whose types and 1-2 bond
   graph reproduce the template; matches are deduplicated by their atom set, so
   symmetric templates are counted once.
------------------------------------------------------------------------- */

int PairTemplateOffset::count_matches(Molecule *mol)
{
  const int na = mol->natoms;
  const tagint maxtag = gref.maxtag;

  // A one-atom (or unbonded) template is just a count of atoms of that type.
  int seed0 = 0, seed1 = -1;
  if (na >= 2 && mol->nspecial[0][0] > 0) seed1 = (int) (mol->special[0][0] - 1);
  if (seed1 < 0) {
    const int t0 = mol->type[0];
    int cnt = 0;
    for (tagint t = 1; t <= maxtag; t++)
      if (gref.type[t] == t0) cnt++;
    return cnt;
  }

  const int type0 = mol->type[seed0];
  const int type1 = mol->type[seed1];
  std::vector<int> is_edge(na, 0);    // full match: no edge atoms
  std::vector<tagint> glove(na, 0);
  std::set<std::vector<tagint>> matches;

  for (tagint t = 1; t <= maxtag; t++) {
    if (gref.type[t] != type0) continue;
    const int nb = gref.nspecial_flat[(size_t) t * 3 + 0];
    for (int b = 0; b < nb; b++) {
      const tagint u = gref.special_flat[(size_t) t * gref.maxspecial + b];
      if (u < 1 || u > maxtag) continue;
      if (gref.type[u] != type1) continue;
      if (!msevb_superimpose(lmp, mol, is_edge.data(), seed0, seed1, t, u, glove.data(), nullptr,
                             &gref))
        continue;
      std::vector<tagint> key(glove.begin(), glove.end());
      std::sort(key.begin(), key.end());
      matches.insert(key);
    }
  }
  return (int) matches.size();
}

/* ---------------------------------------------------------------------- */

void PairTemplateOffset::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  if (!entries.empty()) {
    build_global_topology();

    double e_offset = 0.0;
    for (const auto &e : entries)
      e_offset += e.offset * count_matches(atom->molecules[e.mol_index]);

    // The match count (and therefore e_offset) is identical on every rank; add
    // it once, on rank 0, so the sum LAMMPS takes over ranks equals e_offset.
    if (eflag_global && comm->me == 0) eng_vdwl += e_offset;
  }

  // No forces: the offset is a per-configuration constant.
  if (vflag_fdotr) virial_fdotr_compute();
}
