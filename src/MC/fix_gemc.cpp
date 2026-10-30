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
   Contributing author: Andrew Hong, Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include "fix_gemc.h"

#include "angle.h"
#include "atom.h"
#include "atom_vec.h"
#include "bond.h"
#include "comm.h"
#include "compute.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "improper.h"
#include "kspace.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "neighbor.h"
#include "pair.h"
#include "random_park.h"
#include "universe.h"
#include "update.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <set>
#include <string>
#include <unordered_map>

using namespace LAMMPS_NS;
using namespace FixConst;
using MathConst::MY_PI;

static constexpr double RESTART_VERSION = 4.0;
static constexpr double MAXDLOGVOLRATIO = 1.0;
static constexpr double DEG2RAD = MY_PI / 180.0;
static constexpr double BIG = 1.0e20;
static constexpr int RVOUS = 0;    // 0 for irregular, 1 for all2all communication in rendezvous

/* ---------------------------------------------------------------------- */

FixGEMC::FixGEMC(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), idmol(nullptr), onemol(nullptr), c_pe(nullptr),
    comm_replica(MPI_COMM_NULL), random_universe(nullptr), random_world(nullptr),
    random_proc(nullptr)
{
  if (narg < 11) utils::missing_cmd_args(FLERR, "fix gemc", error);

  // must have exactly two boxes

  if (universe->nworlds != 2)
    error->universe_all(FLERR, "Must use exactly two partitions with fix gemc");

  // various fix flags

  time_integrate = 0;
  global_freq = 1;
  time_depend = 1;
  restart_global = 1;
  vector_flag = 1;
  size_vector = 11;
  extvector = 0;

  // box size changes with volume MC moves

  box_change |= BOX_CHANGE_SIZE;
  if (domain->triclinic) box_change |= BOX_CHANGE_SHAPE;

  // required user args

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  ntranslate = utils::inumeric(FLERR, arg[4], false, lmp);
  nexchange = utils::inumeric(FLERR, arg[5], false, lmp);
  nvolume = utils::inumeric(FLERR, arg[6], false, lmp);
  box_temp = utils::numeric(FLERR, arg[7], false, lmp);
  displace = utils::numeric(FLERR, arg[8], false, lmp);
  max_dlogvolratio = utils::numeric(FLERR, arg[9], false, lmp);
  seed = utils::inumeric(FLERR, arg[10], false, lmp);

  if (nevery <= 0) error->all(FLERR, 3, "Illegal fix gemc N value {}: must be > 0", nevery);
  if (ntranslate < 0) error->all(FLERR, 4, "Illegal fix gemc M value {}: must be >= 0", ntranslate);
  if (nexchange < 0) error->all(FLERR, 5, "Illegal fix gemc X value {}: must be >= 0", nexchange);
  if (nvolume < 0) error->all(FLERR, 6, "Illegal fix gemc V value {}: must be >= 0", nvolume);
  if (ntranslate + nexchange + nvolume <= 0)
    error->all(FLERR, "Fix gemc requires at least one type of Monte Carlo move");
  if (box_temp <= 0.0)
    error->all(FLERR, 7, "Illegal fix gemc temperature {}: must be > 0", box_temp);
  if (displace < 0.0)
    error->all(FLERR, 8, "Illegal fix gemc displace value {}: must be >= 0", displace);
  if (max_dlogvolratio <= 0.0)
    error->all(FLERR, 9, "Illegal fix gemc maxdlogvolratio value {}: must be > 0",
               max_dlogvolratio);
  if ((seed <= 0) || (seed > MAXSMALLINT - 3 - universe->nprocs))
    error->all(FLERR, 10, "Illegal fix gemc random number seed {}", seed);

  // optional keywords

  full_flag = 0;
  molflag = 0;
  maxangle = 30.0 * DEG2RAD;
  tune_every = 0;
  tune_trans = tune_vol = 0.4;
  int iarg = 11;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "full_energy") == 0) {
      full_flag = 1;
      iarg++;
    } else if (strcmp(arg[iarg], "mol") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix gemc mol", error);
      if (atom->find_molecule(arg[iarg + 1]) == -1)
        error->all(FLERR, iarg + 1, "Molecule template ID {} for fix gemc does not exist",
                   arg[iarg + 1]);
      delete[] idmol;
      idmol = utils::strdup(arg[iarg + 1]);
      molflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg], "maxangle") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix gemc maxangle", error);
      maxangle = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if ((maxangle <= 0.0) || (maxangle > 180.0))
        error->all(FLERR, iarg + 1, "Illegal fix gemc maxangle value {}: must be > 0 and <= 180",
                   maxangle);
      maxangle *= DEG2RAD;
      iarg += 2;
    } else if (strcmp(arg[iarg], "tune") == 0) {
      if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "fix gemc tune", error);
      tune_every = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      tune_trans = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      tune_vol = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      if (tune_every <= 0)
        error->all(FLERR, iarg + 1, "Illegal fix gemc tune N value {}: must be > 0", tune_every);
      if ((tune_trans <= 0.0) || (tune_trans >= 1.0))
        error->all(FLERR, iarg + 2, "Illegal fix gemc tune translation acceptance ratio {}",
                   tune_trans);
      if ((tune_vol <= 0.0) || (tune_vol >= 1.0))
        error->all(FLERR, iarg + 3, "Illegal fix gemc tune volume acceptance ratio {}", tune_vol);
      iarg += 4;
    } else {
      error->all(FLERR, iarg, "Unknown fix gemc keyword: {}", arg[iarg]);
    }
  }
  local_flag = 0;
  ghosts_stale = 0;
  charge_warned = 0;
  group_charged = 0;
  natoms_per_molecule = 1;

  if (molflag) {
    if (!atom->molecule_flag)
      error->all(FLERR, "Fix gemc keyword mol requires an atom style with molecule IDs");
    if (atom->molecular == Atom::TEMPLATE)
      error->all(FLERR, "Fix gemc does not support atom style template");
  }

  // set up comm_replica = communicator between the same ranks of both boxes
  // only rank 0 of each box uses it, rank 0 of box 1 is rank 0 of comm_replica

  MPI_Comm_rank(world, &me);
  myworld = universe->iworld;
  MPI_Comm_split(universe->uworld, me, myworld, &comm_replica);

  // random number generators: unique to each rank, synchronized across one box,
  // and synchronized across both boxes. all seeds are different.

  random_universe = new RanPark(lmp, seed);
  random_world = new RanPark(lmp, seed);
  double worldcoord[3] = {(double) myworld, 1.0, 0.0};
  random_world->reset(seed, worldcoord);
  random_proc = new RanPark(lmp, seed);
  random_proc->reset_rank(seed, universe->me);

  ntranslation_attempts = ntranslation_successes = 0.0;
  nrotation_attempts = nrotation_successes = 0.0;
  nvolume_attempts = nvolume_successes = 0.0;
  nexchange_attempts = nexchange_successes = 0.0;
  for (auto &n : nlast) n = 0.0;
  for (auto &n : tune_last) n = 0.0;
  ntune_calls = 0;

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;

  natom_lower = natom_total = 0;
  natom_local = 0;
  grid_valid = 0;
  nstored = 0;
  nbin[0] = nbin[1] = nbin[2] = 1;
  binlo[0] = binlo[1] = binlo[2] = 0.0;
  bininv[0] = bininv[1] = bininv[2] = 1.0;
  ghost_skin = move_limit = 0.0;
  nbase = use_map = 0;
  maxtag_box = 0;
  pending_changes = 0;
  local_warned = 0;
  rvous_count = 0;
  rvous_bad = rvous_charged = 0;
  logvolratio = voltot = 0.0;
  triclinic = 0;
  xlo = ylo = zlo = xhi = yhi = zhi = boxxy = boxxz = boxyz = 0.0;
}

/* ---------------------------------------------------------------------- */

FixGEMC::~FixGEMC()
{
  delete random_proc;
  delete random_world;
  delete random_universe;
  delete[] idmol;
  if (comm_replica != MPI_COMM_NULL) MPI_Comm_free(&comm_replica);
}

/* ---------------------------------------------------------------------- */

int FixGEMC::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGEMC::init()
{
  if (domain->dimension != 3) error->all(FLERR, "Fix gemc requires a 3d system");
  if (domain->nonperiodic) error->all(FLERR, "Fix gemc requires a fully periodic box");
  if (!atom->tag_enable) error->all(FLERR, "Fix gemc requires atom IDs");
  if (atom->rmass_flag)
    error->all(FLERR, "Fix gemc does not support atom styles with per-atom masses");
  if (!atom->mass) error->all(FLERR, "Fix gemc requires per-type masses");

  // exchanges transfer only type, charge, group membership, velocity, and
  // topology, so atom styles with any other per-atom property that moves
  // with an atom would lose it

  // clang-format off
  static const std::set<std::string> transferred = {
    "q", "molecule", "num_bond", "bond_type", "bond_atom",
    "num_angle", "angle_type", "angle_atom1", "angle_atom2", "angle_atom3",
    "num_dihedral", "dihedral_type", "dihedral_atom1", "dihedral_atom2", "dihedral_atom3",
    "dihedral_atom4", "num_improper", "improper_type", "improper_atom1", "improper_atom2",
    "improper_atom3", "improper_atom4", "nspecial", "special"};
  // clang-format on
  if (atom->avec->bonus_flag)
    error->all(FLERR, "Fix gemc does not support atom style {}", atom->atom_style);
  for (const auto &field : atom->avec->fields_exchange)
    if (!transferred.count(field))
      error->all(FLERR, "Fix gemc does not support atom style {} with per-atom property {}",
                 atom->atom_style, field);
  if (atom->nivector || atom->ndvector || atom->niarray || atom->ndarray)
    error->all(FLERR, "Fix gemc does not support custom per-atom properties");
  if (force->pair && force->pair->tail_flag && !force->pair->reinitflag)
    error->all(FLERR, "Fix gemc with pair_modify tail yes is not supported by pair style {}",
               force->pair_style);
  for (const auto &ifix : modify->get_fix_list())
    if (ifix->rigid_flag || utils::strmatch(ifix->style, "^shake") ||
        utils::strmatch(ifix->style, "^rattle"))
      error->all(FLERR, "Fix gemc does not support constraints like fix {}", ifix->style);

  triclinic = domain->triclinic;

  // molecule template

  if (molflag) {
    int imol = atom->find_molecule(idmol);
    if (imol == -1) error->all(FLERR, "Molecule template ID {} for fix gemc does not exist", idmol);
    onemol = atom->molecules[imol];
    if ((onemol->nset > 1) && (comm->me == 0))
      error->warning(FLERR,
                     "Molecule template {} for fix gemc has multiple molecules; "
                     "using only the first",
                     idmol);
    if (!onemol->xflag || !onemol->typeflag)
      error->all(FLERR, "Molecule template {} for fix gemc must define coordinates and types",
                 idmol);
    if (onemol->ntypes > atom->ntypes)
      error->all(FLERR, "Molecule template {} for fix gemc has invalid atom types", idmol);
    if (atom->molecular == Atom::MOLECULAR) {
      if ((onemol->bondflag && !atom->avec->bonds_allow) ||
          (onemol->angleflag && !atom->avec->angles_allow) ||
          (onemol->dihedralflag && !atom->avec->dihedrals_allow) ||
          (onemol->improperflag && !atom->avec->impropers_allow))
        error->all(FLERR,
                   "Molecule template {} for fix gemc has topology not allowed by the "
                   "atom style",
                   idmol);
      if (onemol->specialflag && (onemol->maxspecial > atom->maxspecial))
        error->all(FLERR,
                   "Molecule template {} for fix gemc has too many special neighbors; "
                   "use extra/special/per/atom",
                   idmol);
    }
    onemol->check_attributes();
    natoms_per_molecule = onemol->natoms;
  } else {
    natoms_per_molecule = 1;
  }
  check_molecules();

  // decide whether single-atom energies can replace full energy evaluations
  // for translation and exchange moves. requires a pair style that
  // provides single() and has no many-body terms, no tail corrections,
  // no neighbor list exclusions, and no fixes contributing to the potential energy.
  // displacements must stay within the ghost atom cutoff minus the pair cutoff
  // so that ghost atoms cover all interactions of a displaced atom.

  local_flag = 0;
  if (!full_flag) {
    local_flag = 1;
    std::string reason;
    if (!force->pair) {
      reason = "no pair style";
    } else if (molflag) {
      reason = "molecules are exchanged";
    } else if (force->kspace) {
      reason = "a long-range solver is used";
    } else if (group_charged) {
      reason = "exchanged atoms are charged";
    } else if (!force->pair->single_enable) {
      reason = fmt::format("pair style {} does not support single()", force->pair_style);
    } else if (force->pair->manybody_flag) {
      reason = fmt::format("pair style {} is a many-body potential", force->pair_style);
    } else if (force->pair->tail_flag) {
      reason = "pair_modify tail yes is used";
    } else if (neighbor->nex_type || neighbor->nex_group || neighbor->nex_mol) {
      reason = "neighbor list exclusions are used";
    } else if (displace > init_skin()) {
      reason = "the displace value is larger than the ghost atom cutoff minus the pair cutoff";
    } else {
      for (const auto &ifix : modify->get_fix_list())
        if (ifix->energy_global_flag && ifix->thermo_energy)
          reason = fmt::format("fix {} contributes to the potential energy", ifix->id);
    }

    // moves of molecules always use the total energy, as documented, so no warning

    if (!reason.empty()) {
      local_flag = 0;
      if ((comm->me == 0) && !molflag)
        error->warning(FLERR, "Fix gemc uses full energy evaluations for all moves because {}",
                       reason);
    }
  }

  // total energy is taken from the thermo_pe compute

  c_pe = modify->get_compute_by_id("thermo_pe");
  if (!c_pe) error->all(FLERR, "Fix gemc could not find thermo_pe compute");

  // move type selection probabilities

  nmoves = nvolume + nexchange + ntranslate;
  pc_exchange = static_cast<double>(nexchange) / nmoves;
  pc_volume = static_cast<double>(nexchange + nvolume) / nmoves;

  beta = 1.0 / (force->boltz * box_temp);

  // update box dimensions and list of atoms in the fix group

  xlo = domain->boxlo[0];
  xhi = domain->boxhi[0];
  ylo = domain->boxlo[1];
  yhi = domain->boxhi[1];
  zlo = domain->boxlo[2];
  zhi = domain->boxhi[2];
  boxxy = domain->xy;
  boxxz = domain->xz;
  boxyz = domain->yz;

  update_gas_atoms_list();
  progress = 0;
}

/* ----------------------------------------------------------------------
   check that atoms in the fix group can be moved and exchanged:
   without mol keyword, the group atoms must not belong to molecules;
   with mol keyword, each group molecule must match the template:
   consecutive atom IDs in the order of the template atoms, same types.
   also warn once if exchanged atoms or molecules carry a net charge.
------------------------------------------------------------------------- */

void FixGEMC::check_molecules()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  double *q = atom->q;

  int flag = 0;
  int charged = 0;
  group_charged = 0;

  if (!molflag) {
    int anyq = 0;
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      if (molecule && (molecule[i] > 0)) flag = 1;
      if (q && (q[i] != 0.0)) charged = anyq = 1;
    }
    MPI_Allreduce(&anyq, &group_charged, 1, MPI_INT, MPI_MAX, world);
    int flag_all;
    MPI_Allreduce(&flag, &flag_all, 1, MPI_INT, MPI_MAX, world);
    if (flag_all)
      error->all(FLERR, "Fix gemc group contains atoms with molecule IDs; use the mol keyword");
  } else {

    // check each molecule on a rendezvous rank (molecule ID modulo number of ranks)

    int nsend = 0;
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      if (molecule[i] <= 0)
        flag = 1;
      else
        nsend++;
    }
    int flag_all;
    MPI_Allreduce(&flag, &flag_all, 1, MPI_INT, MPI_MAX, world);
    if (flag_all)
      error->all(FLERR, "Fix gemc group atoms must all belong to molecules with the mol keyword");

    int *proclist;
    memory->create(proclist, nsend, "gemc:proclist");
    auto *inbuf = (CheckRvous *) memory->smalloc((bigint) nsend * sizeof(CheckRvous), "gemc:inbuf");
    int n = 0;
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit) || (molecule[i] <= 0)) continue;
      proclist[n] = molecule[i] % comm->nprocs;
      inbuf[n].mol = molecule[i];
      inbuf[n].tag = atom->tag[i];
      inbuf[n].type = atom->type[i];
      inbuf[n].ntopo[0] = atom->num_bond ? atom->num_bond[i] : 0;
      inbuf[n].ntopo[1] = atom->num_angle ? atom->num_angle[i] : 0;
      inbuf[n].ntopo[2] = atom->num_dihedral ? atom->num_dihedral[i] : 0;
      inbuf[n].ntopo[3] = atom->num_improper ? atom->num_improper[i] : 0;
      inbuf[n].q = q ? q[i] : 0.0;
      n++;
    }
    char *buf;
    rvous_bad = rvous_charged = 0;
    comm->rendezvous(RVOUS, nsend, (char *) inbuf, sizeof(CheckRvous), 0, proclist,
                     rendezvous_check, 0, buf, 0, (void *) this);
    memory->destroy(proclist);
    memory->sfree(inbuf);

    MPI_Allreduce(&rvous_bad, &flag_all, 1, MPI_INT, MPI_MAX, world);
    if (flag_all)
      error->all(FLERR,
                 "Molecules in fix gemc group must match molecule template {}: same number and "
                 "types of atoms, with consecutive atom IDs in template order, and the same "
                 "number of bonds, angles, dihedrals, and impropers",
                 idmol);
    charged = rvous_charged;
    if (onemol->qflag) {
      double qmol = 0.0;
      for (int i = 0; i < onemol->natoms; i++) qmol += onemol->q[i];
      if (fabs(qmol) > 1.0e-6) charged = 1;
    }
  }

  int charged_all;
  MPI_Allreduce(&charged, &charged_all, 1, MPI_INT, MPI_MAX, world);
  if (charged_all && !charge_warned && (comm->me == 0))
    error->warning(FLERR,
                   "Fix gemc exchanges {} with a net charge, so the boxes will not be "
                   "charge neutral",
                   molflag ? "molecules" : "atoms");
  if (charged_all) charge_warned = 1;
}

/* ----------------------------------------------------------------------
   callback of comm->rendezvous() for check_molecules():
   check that the atoms of each molecule match the molecule template, with
   consecutive atom IDs in template order, and whether it has a net charge
------------------------------------------------------------------------- */

int FixGEMC::rendezvous_check(int n, char *inbuf, int &flag, int *& /*proclist*/,
                              char *& /*outbuf*/, void *ptr)
{
  auto *fptr = (FixGEMC *) ptr;
  auto *in = (CheckRvous *) inbuf;
  const int natoms = fptr->natoms_per_molecule;
  Molecule *onemol = fptr->onemol;

  // numbers of bonds, angles, dihedrals, and impropers of the template, counted
  // like those of the atoms, which follow the same newton_bond setting

  int ntopo_mol[4] = {0, 0, 0, 0};
  for (int i = 0; i < natoms; i++) {
    if (onemol->bondflag) ntopo_mol[0] += onemol->num_bond[i];
    if (onemol->angleflag) ntopo_mol[1] += onemol->num_angle[i];
    if (onemol->dihedralflag) ntopo_mol[2] += onemol->num_dihedral[i];
    if (onemol->improperflag) ntopo_mol[3] += onemol->num_improper[i];
  }

  std::sort(in, in + n, [](const CheckRvous &a, const CheckRvous &b) {
    return (a.mol < b.mol) || ((a.mol == b.mol) && (a.tag < b.tag));
  });

  int k = 0;
  while (k < n) {
    int kend = k;
    while ((kend < n) && (in[kend].mol == in[k].mol)) kend++;
    if (kend - k != natoms) fptr->rvous_bad = 1;
    double qmol = 0.0;
    int ntopo[4] = {0, 0, 0, 0};
    for (int m = k; m < kend; m++) {
      for (int t = 0; t < 4; t++) ntopo[t] += in[m].ntopo[t];
      int iatom = m - k;
      if (iatom < natoms) {
        if (in[m].type != onemol->type[iatom]) fptr->rvous_bad = 1;
        if (in[m].tag != in[k].tag + iatom) fptr->rvous_bad = 1;
      }
      qmol += in[m].q;
    }
    for (int t = 0; t < 4; t++)
      if (ntopo[t] != ntopo_mol[t]) fptr->rvous_bad = 1;
    if (fabs(qmol) > 1.0e-6) fptr->rvous_charged = 1;
    k = kend;
  }

  // flag = 0: no second communication

  flag = 0;
  return 0;
}

/* ----------------------------------------------------------------------
   1 if an atom in the fix group is part of a bond, angle, dihedral, or
   improper owned by this rank.  each of them is stored with one of its atoms
   (or all of them with newton_bond off), and all of its atoms are owned or
   ghost atoms of the rank that stores it.
------------------------------------------------------------------------- */

int FixGEMC::group_topology()
{
  if (!atom->molecular) return 0;
  int *mask = atom->mask;

  auto ingroup = [&](tagint itag) {
    int j = atom->map(itag);
    return (j >= 0) && (mask[j] & groupbit);
  };

  for (int i = 0; i < atom->nlocal; i++) {
    int self = mask[i] & groupbit;
    if (atom->num_bond)
      for (int m = 0; m < atom->num_bond[i]; m++)
        if (self || ingroup(atom->bond_atom[i][m])) return 1;
    if (atom->num_angle)
      for (int m = 0; m < atom->num_angle[i]; m++)
        if (ingroup(atom->angle_atom1[i][m]) || ingroup(atom->angle_atom2[i][m]) ||
            ingroup(atom->angle_atom3[i][m]))
          return 1;
    if (atom->num_dihedral)
      for (int m = 0; m < atom->num_dihedral[i]; m++)
        if (ingroup(atom->dihedral_atom1[i][m]) || ingroup(atom->dihedral_atom2[i][m]) ||
            ingroup(atom->dihedral_atom3[i][m]) || ingroup(atom->dihedral_atom4[i][m]))
          return 1;
    if (atom->num_improper)
      for (int m = 0; m < atom->num_improper[i]; m++)
        if (ingroup(atom->improper_atom1[i][m]) || ingroup(atom->improper_atom2[i][m]) ||
            ingroup(atom->improper_atom3[i][m]) || ingroup(atom->improper_atom4[i][m]))
          return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   checks and initialization that require communication between the boxes.
   done in setup() and not in init(), since only setup() is guaranteed to be
   called at the same time in both partitions (at the start of a run).
------------------------------------------------------------------------- */

void FixGEMC::setup(int /*vflag*/)
{
  // without the mol keyword, exchanged atoms must not be bonded to other atoms.
  // checked here, since the ghost atoms that hold the bond partners exist now

  if (!molflag) {
    int flag = group_topology();
    int flag_all;
    MPI_Allreduce(&flag, &flag_all, 1, MPI_INT, MPI_MAX, world);
    if (flag_all)
      error->all(FLERR,
                 "Fix gemc group contains atoms with bonds, angles, dihedrals, or impropers; "
                 "use the mol keyword");
  }

  // both boxes must use the same settings, or they would fall out of step and hang

  // this includes the state of the synchronized random number generator
  // and the group definitions, since group bitmasks are transferred between boxes

  uint32_t hash = 2166136261u;
  for (int ig = 0; ig < Group::MAX_GROUP; ig++) {
    std::string name = group->names[ig] ? group->names[ig] : "";
    name += '\n';
    for (const auto c : name) hash = (hash ^ (uint32_t) (unsigned char) c) * 16777619u;
  }

  // step size adjustments only use moves attempted during this run

  tune_last[0] = ntranslation_attempts;
  tune_last[1] = ntranslation_successes;
  tune_last[2] = nvolume_attempts;
  tune_last[3] = nvolume_successes;
  tune_last[4] = nrotation_attempts;
  tune_last[5] = nrotation_successes;
  ntune_calls = 0;

  // the maximum displacement and rotation may differ between the boxes

  // signature of atom style and molecule template: types, topology, charges

  double molsig = atom->molecular + 3.0 * atom->q_flag;
  if (molflag) {
    molsig += 7.0 * onemol->nbonds + 11.0 * onemol->nangles + 13.0 * onemol->ndihedrals +
        17.0 * onemol->nimpropers;
    for (int i = 0; i < onemol->natoms; i++) {
      molsig += 19.0 * (i + 1) * onemol->type[i];
      if (onemol->qflag) molsig += 23.0 * (i + 1) * onemol->q[i];
    }
  }

  constexpr int NCHECK = 19;
  int mismatch = 0;
  if (me == 0) {
    double mine[NCHECK] = {(double) nevery,
                           (double) ntranslate,
                           (double) nexchange,
                           (double) nvolume,
                           box_temp,
                           (double) tune_every,
                           tune_trans,
                           tune_vol,
                           max_dlogvolratio,
                           (double) seed,
                           (double) atom->ntypes,
                           (double) update->ntimestep,
                           (double) next_reneighbor,
                           (double) random_universe->state(),
                           (double) hash,
                           (double) molflag,
                           (double) natoms_per_molecule,
                           molsig,
                           (double) igroup};
    double other[NCHECK];
    MPI_Sendrecv(mine, NCHECK, MPI_DOUBLE, 1 - myworld, 0, other, NCHECK, MPI_DOUBLE, 1 - myworld,
                 0, comm_replica, MPI_STATUS_IGNORE);
    for (int i = 0; i < NCHECK; i++)
      if (mine[i] != other[i]) mismatch = 1;
  }
  MPI_Bcast(&mismatch, 1, MPI_INT, 0, world);
  if (mismatch)
    error->universe_all(FLERR,
                        "Fix gemc settings, groups, number of atom types, atom style, molecule "
                        "template, timestep, and restart status must be the same in both "
                        "partitions");

  // atoms are exchanged between the boxes, so single-atom energies are not
  // possible in either box if the exchanged atoms in one box are charged

  if (any_box(group_charged) && local_flag) {
    local_flag = 0;
    if (comm->me == 0)
      error->warning(FLERR,
                     "Fix gemc uses full energy evaluations for all moves because "
                     "exchanged atoms are charged");
  }

  // initialize log volume ratio and total volume

  double vol_i = box_volume();
  double vol_j = 0.0;
  if (me == 0)
    MPI_Sendrecv(&vol_i, 1, MPI_DOUBLE, 1 - myworld, 0, &vol_j, 1, MPI_DOUBLE, 1 - myworld, 0,
                 comm_replica, MPI_STATUS_IGNORE);
  MPI_Bcast(&vol_j, 1, MPI_DOUBLE, 0, world);

  voltot = vol_i + vol_j;
  if (myworld == 0)
    logvolratio = log(vol_i / vol_j);
  else
    logvolratio = log(vol_j / vol_i);
}

/* ----------------------------------------------------------------------
   attempt Monte Carlo translations, exchanges, and volume changes
   done before exchange, borders, reneighbor
   so that ghost atoms and neighbor lists will be correct
------------------------------------------------------------------------- */

void FixGEMC::pre_exchange()
{
  // just return if should not be called on this timestep

  if (next_reneighbor != update->ntimestep) return;

  xlo = domain->boxlo[0];
  xhi = domain->boxhi[0];
  ylo = domain->boxlo[1];
  yhi = domain->boxhi[1];
  zlo = domain->boxlo[2];
  zhi = domain->boxhi[2];
  boxxy = domain->xy;
  boxxz = domain->xz;
  boxyz = domain->yz;

  next_reneighbor = update->ntimestep + nevery;

  // the move type sequence comes from the synchronized RNG,
  // so both boxes always attempt the same type of move

  energy_stored = energy_full();

  for (int i = 0; i < nmoves; i++) {
    double imove = random_universe->uniform();
    if (imove < pc_exchange) {
      if (molflag)
        attempt_molecule_exchange_full();
      else
        attempt_atomic_exchange_full();
    } else if (imove < pc_volume) {
      attempt_volume_change_full();
    } else if (molflag) {
      // translation or rotation, chosen independently in each box
      if (random_world->uniform() < 0.5)
        attempt_molecule_translation_full();
      else
        attempt_molecule_rotation_full();
    } else {
      attempt_atomic_translation_full();
    }
  }

  // atoms inserted or removed with single-atom energies must exist
  // in the atom arrays for the following time steps

  flush_pending();

  if (tune_every && (++ntune_calls % tune_every == 0)) tune_steps();

  print_progress();
}

/* ----------------------------------------------------------------------
   adjust maximum displacement (separately for each box) and maximum
   change of log(V1/V2) (identical in both boxes, since volume moves are
   accepted jointly) towards the target acceptance ratios, based on the
   moves attempted since the last adjustment
------------------------------------------------------------------------- */

void FixGEMC::tune_steps()
{
  static constexpr int MINATTEMPTS = 20;

  double dtrans = ntranslation_attempts - tune_last[0];
  if (dtrans >= MINATTEMPTS) {
    double ratio = (ntranslation_successes - tune_last[1]) / dtrans;
    double factor = MIN(MAX(ratio / tune_trans, 0.5), 1.5);
    double maxdisp = 0.5 * min_box_width();
    if (local_flag) maxdisp = MIN(maxdisp, 0.5 * local_skin());
    displace = MIN(MIN(displace * factor, maxdisp), max_move());
    displace = MAX(displace, 1.0e-6 * maxdisp);
    tune_last[0] = ntranslation_attempts;
    tune_last[1] = ntranslation_successes;
  }

  double dvol = nvolume_attempts - tune_last[2];
  if (dvol >= MINATTEMPTS) {
    double ratio = (nvolume_successes - tune_last[3]) / dvol;
    double factor = MIN(MAX(ratio / tune_vol, 0.5), 1.5);
    max_dlogvolratio = MIN(max_dlogvolratio * factor, MAXDLOGVOLRATIO);
    max_dlogvolratio = MAX(max_dlogvolratio, 1.0e-6);
    tune_last[2] = nvolume_attempts;
    tune_last[3] = nvolume_successes;
  }

  double drot = nrotation_attempts - tune_last[4];
  if (drot >= MINATTEMPTS) {
    double ratio = (nrotation_successes - tune_last[5]) / drot;
    double factor = MIN(MAX(ratio / tune_trans, 0.5), 1.5);
    maxangle = MIN(maxangle * factor, MY_PI);
    maxangle = MAX(maxangle, 1.0e-6);
    tune_last[4] = nrotation_attempts;
    tune_last[5] = nrotation_successes;
  }
}

/* ----------------------------------------------------------------------
   print progress info to universe screen/logfile, about every 1% of a run
------------------------------------------------------------------------- */

void FixGEMC::print_progress()
{
  // rank 0 of both boxes follows the progress, since the boxes are at the
  // same timestep.  world 0 needs the atom count of world 1 for the message

  if (me != 0) return;

  double delta = update->ntimestep - update->beginstep;
  if ((delta != 0.0) && (update->beginstep != update->endstep))
    delta /= update->endstep - update->beginstep;
  int status = static_cast<int>(delta * 100.0);
  if (status <= progress) return;
  progress = status;

  bigint n_other = 0;
  MPI_Sendrecv(&natom_total, 1, MPI_LMP_BIGINT, 1 - myworld, 0, &n_other, 1, MPI_LMP_BIGINT,
               1 - myworld, 0, comm_replica, MPI_STATUS_IGNORE);

  if (universe->me != 0) return;

  double now[6] = {ntranslation_successes, ntranslation_attempts, nvolume_successes,
                   nvolume_attempts,       nexchange_successes,   nexchange_attempts};
  double d[6];
  for (int i = 0; i < 6; i++) {
    d[i] = now[i] - nlast[i];
    nlast[i] = now[i];
  }

  double vol1 = voltot / (1.0 + exp(-logvolratio));
  double vol2 = voltot / (1.0 + exp(logvolratio));
  bigint n1 = natom_total / natoms_per_molecule;
  bigint n2 = n_other / natoms_per_molecule;

  auto msg =
      fmt::format(" GEMC run progress: {:>3d}% \n  Trans: {:g}/{:g}\n"
                  "  Vol: {:g}/{:g}\n  Ex: {:g}/{:g}\n"
                  "  Replica Volume N{} Number-density:\n"
                  "   1: {:g} {:d} {:g}\n"
                  "   2: {:g} {:d} {:g}\n",
                  progress, d[0], d[1], d[2], d[3], d[4], d[5], molflag ? "molecules" : "particles",
                  vol1, n1, n1 / vol1, vol2, n2, n2 / vol2);
  if (universe->uscreen) utils::print(universe->uscreen, msg);
  if (universe->ulogfile) utils::print(universe->ulogfile, msg);
}

/* ----------------------------------------------------------------------
   update list of local atoms in fix group and their counts
------------------------------------------------------------------------- */

void FixGEMC::update_gas_atoms_list()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  if ((int) gas_pos.size() < atom->nmax) gas_pos.resize(atom->nmax, -1);
  for (int i : gas_list) gas_pos[i] = -1;
  gas_list.clear();
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) gas_add(i);
  natom_local = gas_list.size();

  bigint nmine = natom_local;
  MPI_Allreduce(&nmine, &natom_total, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  MPI_Scan(&nmine, &natom_lower, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  natom_lower -= nmine;
}

/* ----------------------------------------------------------------------
   add local index i to the list of local group atoms
   the counts on other ranks must be updated by the caller
------------------------------------------------------------------------- */

void FixGEMC::gas_add(int i)
{
  gas_pos[i] = gas_list.size();
  gas_list.push_back(i);
  natom_local = gas_list.size();
}

/* ----------------------------------------------------------------------
   remove local index i from the list of local group atoms
   the counts on other ranks must be updated by the caller
------------------------------------------------------------------------- */

void FixGEMC::gas_remove(int i)
{
  int pos = gas_pos[i];
  int last = gas_list.back();
  gas_list[pos] = last;
  gas_pos[last] = pos;
  gas_list.pop_back();
  gas_pos[i] = -1;
  natom_local = gas_list.size();
}

/* ----------------------------------------------------------------------
   return local index of a randomly chosen group atom in my box
   or -1 if not owned by this rank. must be called by all ranks of the box.
------------------------------------------------------------------------- */

int FixGEMC::pick_random_gas_atom()
{
  int i = -1;
  auto iwhichglobal = static_cast<bigint>(natom_total * random_world->uniform());
  if ((iwhichglobal >= natom_lower) && (iwhichglobal < natom_lower + natom_local))
    i = gas_list[iwhichglobal - natom_lower];
  return i;
}

/* ----------------------------------------------------------------------
   return 1 if flag (valid on rank 0) is set in either box, on all ranks of both boxes
------------------------------------------------------------------------- */

int FixGEMC::any_box(int flag)
{
  int result = 0;
  if (me == 0) {
    int other;
    MPI_Sendrecv(&flag, 1, MPI_INT, 1 - myworld, 0, &other, 1, MPI_INT, 1 - myworld, 0,
                 comm_replica, MPI_STATUS_IGNORE);
    result = (flag || other) ? 1 : 0;
  }
  MPI_Bcast(&result, 1, MPI_INT, 0, world);
  return result;
}

/* ----------------------------------------------------------------------
   joint Metropolis decision for a move that changes both boxes
   dU = my box's contribution to the effective energy change (valid on rank 0)
   overflow = 1 if my box's trial energy is unusable (valid on rank 0)
   both boxes compute the same sum and draw the same random number,
   so they always reach the same decision
------------------------------------------------------------------------- */

int FixGEMC::accept_both(double dU, int overflow)
{
  double rnd = random_universe->uniform();
  int success = 0;
  if (me == 0) {
    double mine[2] = {dU, (double) overflow};
    double other[2];
    MPI_Sendrecv(mine, 2, MPI_DOUBLE, 1 - myworld, 0, other, 2, MPI_DOUBLE, 1 - myworld, 0,
                 comm_replica, MPI_STATUS_IGNORE);
    double all_dU = mine[0] + other[0];
    if ((mine[1] == 0.0) && (other[1] == 0.0) && std::isfinite(all_dU))
      success = (all_dU <= 0.0) || (rnd < exp(-beta * all_dU));
  }
  MPI_Bcast(&success, 1, MPI_INT, 0, world);
  return success;
}

/* ----------------------------------------------------------------------
   move atoms to their owning ranks, rebuild ghost atoms and neighbor lists
------------------------------------------------------------------------- */

void FixGEMC::reset_comm()
{
  refresh_ghosts();
  if (modify->has_fixes(FixConst::PRE_NEIGHBOR)) modify->pre_neighbor();
  neighbor->build(1);
}

/* ----------------------------------------------------------------------
   move atoms to their owning ranks and rebuild ghost atoms, no neighbor lists
   atoms are reordered, so the list of local group atoms is regenerated
------------------------------------------------------------------------- */

void FixGEMC::refresh_ghosts()
{
  flush_pending();
  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  atom->nghost = 0;
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
  ghosts_stale = 0;
  grid_valid = 0;
  update_gas_atoms_list();
}

/* ----------------------------------------------------------------------
   return 1 if a translation or exchange move in my box can use energy_local()
   the box must be larger than the pair cutoff, so that an atom does
   not interact with its own periodic images, and a displacement must not
   reach beyond the ghost atoms. must be called by all ranks of the box.
------------------------------------------------------------------------- */

int FixGEMC::use_local()
{
  if (!local_flag || (min_box_width() <= force->pair->cutforce)) {
    flush_pending();
    return 0;
  }
  if (displace > local_skin()) {
    if (!local_warned && (comm->me == 0))
      error->warning(FLERR,
                     "Fix gemc uses full energy evaluations for translations and exchanges, "
                     "since the displace value {} is larger than the ghost atom cutoff minus the "
                     "pair cutoff {}",
                     displace, local_skin());
    local_warned = 1;
    flush_pending();
    return 0;
  }
  if (ghosts_stale) refresh_ghosts();
  if (!grid_valid) build_grid();
  return 1;
}

/* ----------------------------------------------------------------------
   estimate of local_skin() in init(), before neighbor lists and
   communication are set up: comm->setup() uses the larger of the neighbor
   cutoff and the user specified communication cutoff
------------------------------------------------------------------------- */

double FixGEMC::init_skin()
{
  double skin = neighbor->skin;
  if ((comm->mode == Comm::SINGLE) && (comm->cutghostuser > force->pair->cutforce + skin))
    skin = comm->cutghostuser - force->pair->cutforce;
  return skin;
}

/* ----------------------------------------------------------------------
   ghost atoms cover all atoms within the ghost cutoff of a subdomain, so an
   atom may move by the ghost cutoff minus the pair cutoff before its
   interactions are no longer covered. the ghost cutoff is in lamda units
   for triclinic boxes. with cutoffs that depend on the atom type, only the
   neighbor skin is guaranteed.
------------------------------------------------------------------------- */

double FixGEMC::local_skin()
{
  double skin = neighbor->skin;
  if (comm->mode == Comm::SINGLE) {
    double cut = comm->cutghost[0];
    if (triclinic) {
      double w[3];
      box_widths(w);
      cut = MIN(MIN(comm->cutghost[0] * w[0], comm->cutghost[1] * w[1]), comm->cutghost[2] * w[2]);
    }
    // the ghost cutoff is at least the pair cutoff plus the neighbor skin,
    // which the difference may miss by round-off

    skin = MAX(skin, cut - force->pair->cutforce);
  }
  return skin;
}

/* ----------------------------------------------------------------------
   pair energy of atom with local index i, type itype, and atom ID itag
   placed at coord with all owned, ghost, and pending atoms, except for its
   own images. only the bins next to the bin of coord are searched, since
   bins are at least as large as the pair cutoff.
   i may be a scratch index beyond the stored atoms for an atom not yet inserted.
------------------------------------------------------------------------- */

double FixGEMC::energy_local(int i, int itype, tagint itag, double *coord)
{
  double **x = atom->x;
  int *type = atom->type;
  tagint *tag = atom->tag;
  Pair *pair = force->pair;
  double *cutsqi = pair->cutsq[itype];
  double fpair;

  int ib[3];
  bin_index(coord, ib);

  double total_energy = 0.0;
  for (int bz = MAX(ib[2] - 1, 0); bz <= MIN(ib[2] + 1, nbin[2] - 1); bz++) {
    for (int by = MAX(ib[1] - 1, 0); by <= MIN(ib[1] + 1, nbin[1] - 1); by++) {
      for (int bx = MAX(ib[0] - 1, 0); bx <= MIN(ib[0] + 1, nbin[0] - 1); bx++) {
        for (int j : bins[(bz * nbin[1] + by) * nbin[0] + bx]) {
          if (tag[j] == itag) continue;
          int jtype = type[j];
          double delx = coord[0] - x[j][0];
          double dely = coord[1] - x[j][1];
          double delz = coord[2] - x[j][2];
          double rsq = delx * delx + dely * dely + delz * delz;
          if (rsq < cutsqi[jtype])
            total_energy += pair->single(i, j, itype, jtype, rsq, 1.0, 1.0, fpair);
        }
      }
    }
  }
  return total_energy;
}

/* ----------------------------------------------------------------------
   sort owned and ghost atoms into bins at least as large as the pair cutoff
   must be called by all ranks of the box right after the ghost atoms were built
------------------------------------------------------------------------- */

void FixGEMC::build_grid()
{
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double **x = atom->x;
  tagint *tag = atom->tag;

  ghost_skin = local_skin();
  move_limit = max_move();

  // grid covers all stored atoms, bins are at least as large as the pair cutoff
  // and there are not many more bins than atoms

  double lo[3], hi[3];
  for (int d = 0; d < 3; d++) {
    lo[d] = BIG;
    hi[d] = -BIG;
  }
  for (int i = 0; i < nall; i++) {
    for (int d = 0; d < 3; d++) {
      lo[d] = MIN(lo[d], x[i][d]);
      hi[d] = MAX(hi[d], x[i][d]);
    }
  }
  double binsize = force->pair->cutforce;
  bigint maxbins = MAX(1000, 8 * (bigint) nall);
  for (int d = 0; d < 3; d++) {
    if (hi[d] < lo[d]) lo[d] = hi[d] = 0.0;
    binlo[d] = lo[d];
    nbin[d] = MAX(1, static_cast<int>(MIN((hi[d] - lo[d]) / binsize, 1.0e6)));
  }
  while ((bigint) nbin[0] * nbin[1] * nbin[2] > maxbins)
    for (int d = 0; d < 3; d++) nbin[d] = MAX(1, nbin[d] / 2);
  for (int d = 0; d < 3; d++)
    bininv[d] = (hi[d] > lo[d]) ? nbin[d] / (hi[d] - lo[d]) : 1.0 / binsize;

  bins.resize((size_t) nbin[0] * nbin[1] * nbin[2]);
  for (auto &b : bins) b.clear();

  nstored = 0;
  grow_stored(nall);
  nstored = nbase = nall;

  // images of an atom ID: from the atom map, which comm->borders() has just built,
  // plus stored atoms in taghead. without an atom map, taghead has all of them.

  use_map = (atom->map_style != Atom::MAP_NONE) ? 1 : 0;
  taghead.clear();
  if (!use_map) taghead.reserve(nall);
  for (int i = nall - 1; i >= 0; i--) {
    bin_add(i);
    if (!use_map) {
      auto it = taghead.find(tag[i]);
      tagnext[i] = (it == taghead.end()) ? -1 : it->second;
      taghead[tag[i]] = i;
    }
  }
  std::fill(dacc.begin(), dacc.begin() + 3 * nall, 0.0);
  pending.clear();
  removed.clear();

  tagint maxtag = 0;
  for (int i = 0; i < nlocal; i++) maxtag = MAX(maxtag, tag[i]);
  MPI_Allreduce(&maxtag, &maxtag_box, 1, MPI_LMP_TAGINT, MPI_MAX, world);

  grid_valid = 1;
}

/* ----------------------------------------------------------------------
   bin indices of a point, points outside of the grid are assigned to the
   closest bin. points closer than a bin size differ by at most 1 in each index.
------------------------------------------------------------------------- */

void FixGEMC::bin_index(double *coord, int *ib)
{
  for (int d = 0; d < 3; d++) {
    double r = (coord[d] - binlo[d]) * bininv[d];
    ib[d] = (r < 0.0) ? 0 : ((r >= nbin[d]) ? nbin[d] - 1 : static_cast<int>(r));
  }
}

/* ----------------------------------------------------------------------
   bin of a point
------------------------------------------------------------------------- */

int FixGEMC::coord2bin(double *coord)
{
  int ib[3];
  bin_index(coord, ib);
  return (ib[2] * nbin[1] + ib[1]) * nbin[0] + ib[0];
}

/* ---------------------------------------------------------------------- */

void FixGEMC::bin_add(int i)
{
  int b = coord2bin(atom->x[i]);
  atombin[i] = b;
  atombinpos[i] = bins[b].size();
  bins[b].push_back(i);
}

/* ---------------------------------------------------------------------- */

void FixGEMC::bin_remove(int i)
{
  int b = atombin[i];
  if (b < 0) return;
  int pos = atombinpos[i];
  int last = bins[b].back();
  bins[b][pos] = last;
  atombinpos[last] = pos;
  bins[b].pop_back();
  atombin[i] = -1;
}

/* ----------------------------------------------------------------------
   make room in the atom arrays and in the per-atom data of this fix
   for local indices up to n-1
------------------------------------------------------------------------- */

void FixGEMC::grow_stored(int n)
{
  while (n > atom->nmax) atom->avec->grow(0);
  size_t nmax = atom->nmax;
  if (atombin.size() < nmax) {
    atombin.resize(nmax, -1);
    atombinpos.resize(nmax, 0);
    tagnext.resize(nmax, -1);
    dacc.resize(3 * nmax, 0.0);
  }
  if (gas_pos.size() < nmax) gas_pos.resize(nmax, -1);
}

/* ----------------------------------------------------------------------
   store an inserted atom with velocity iv or one of its images after the
   stored atoms, returns its local index
------------------------------------------------------------------------- */

int FixGEMC::stored_atom(tagint itag, int itype, int imask, double iq, double *coord, double *iv)
{
  int k = nstored;
  grow_stored(k + 1);
  atom->x[k][0] = coord[0];
  atom->x[k][1] = coord[1];
  atom->x[k][2] = coord[2];
  atom->v[k][0] = iv[0];
  atom->v[k][1] = iv[1];
  atom->v[k][2] = iv[2];
  atom->type[k] = itype;
  atom->mask[k] = imask;
  atom->tag[k] = itag;
  if (atom->q_flag) atom->q[k] = iq;
  if (atom->molecule_flag) atom->molecule[k] = 0;
  bin_add(k);
  auto it = taghead.find(itag);
  tagnext[k] = (it == taghead.end()) ? -1 : it->second;
  taghead[itag] = k;
  dacc[3 * k] = dacc[3 * k + 1] = dacc[3 * k + 2] = 0.0;
  gas_pos[k] = -1;
  nstored++;
  return k;
}

/* ----------------------------------------------------------------------
   local indices of all images of the atom with atom ID itag in the grid
------------------------------------------------------------------------- */

void FixGEMC::find_images(tagint itag, std::vector<int> &images)
{
  images.clear();

  // atom IDs of atoms inserted since the atom map was built are not in it

  if (use_map && (itag <= atom->map_tag_max)) {
    for (int j = atom->map(itag); (j >= 0) && (j < nbase) && (atom->tag[j] == itag);
         j = atom->sametag[j])
      images.push_back(j);
  }
  auto it = taghead.find(itag);
  if (it != taghead.end())
    for (int j = it->second; j >= 0; j = tagnext[j]) images.push_back(j);
}

/* ----------------------------------------------------------------------
   remove all images of the atom with atom ID itag from the grid
------------------------------------------------------------------------- */

void FixGEMC::remove_images(tagint itag)
{
  find_images(itag, image_list);
  for (int j : image_list) bin_remove(j);
  taghead.erase(itag);
}

/* ----------------------------------------------------------------------
   displace all images of the atom with atom ID itag by delta
------------------------------------------------------------------------- */

void FixGEMC::move_images(tagint itag, double *delta)
{
  find_images(itag, image_list);
  double **x = atom->x;
  for (int j : image_list) {
    bin_remove(j);
    x[j][0] += delta[0];
    x[j][1] += delta[1];
    x[j][2] += delta[2];
    bin_add(j);
  }
}

/* ----------------------------------------------------------------------
   insert an atom at coord in the box: the owning rank (owner = 1) stores it
   as a pending atom, and every rank stores its periodic images within its
   ghost cutoff. the atom is created in the atom arrays by flush_pending()
------------------------------------------------------------------------- */

void FixGEMC::insert_pending(tagint itag, int itype, int imask, double iq, double *coord,
                             double *iv, int owner)
{
  if (owner) {
    int k = stored_atom(itag, itype, imask, iq, coord, iv);
    pending.push_back(k);
    if (imask & groupbit) gas_add(k);
  }
  add_images(itag, itype, imask, iq, coord, iv);
}

/* ----------------------------------------------------------------------
   store all periodic images of the atom with atom ID itag at coord that lie
   within the ghost cutoff of my subdomain, like comm->borders() would, and
   are not stored yet. with this, the stored atoms of each rank always
   include all atoms within its ghost cutoff.
------------------------------------------------------------------------- */

void FixGEMC::add_images(tagint itag, int itype, int imask, double iq, double *coord, double *iv)
{
  // images in box or lamda coordinates, the ghost cutoff has the same units

  double *lo, *hi, period[3], p[3];
  if (triclinic) {
    domain->x2lamda(coord, p);
    lo = domain->sublo_lamda;
    hi = domain->subhi_lamda;
    period[0] = period[1] = period[2] = 1.0;
  } else {
    p[0] = coord[0];
    p[1] = coord[1];
    p[2] = coord[2];
    lo = domain->sublo;
    hi = domain->subhi;
    period[0] = domain->prd[0];
    period[1] = domain->prd[1];
    period[2] = domain->prd[2];
  }
  double *cut = comm->cutghost;

  // periodic shifts of the images that are already stored

  find_images(itag, image_list);
  std::vector<std::array<int, 3>> &stored = image_shifts;
  stored.clear();
  for (int j : image_list) {
    double q[3];
    if (triclinic)
      domain->x2lamda(atom->x[j], q);
    else {
      q[0] = atom->x[j][0];
      q[1] = atom->x[j][1];
      q[2] = atom->x[j][2];
    }
    std::array<int, 3> k;
    for (int d = 0; d < 3; d++) k[d] = static_cast<int>(lround((q[d] - p[d]) / period[d]));
    stored.push_back(k);
  }

  int kmax[3];
  for (int d = 0; d < 3; d++) kmax[d] = static_cast<int>(cut[d] / period[d]) + 1;
  for (int kz = -kmax[2]; kz <= kmax[2]; kz++) {
    for (int ky = -kmax[1]; ky <= kmax[1]; ky++) {
      for (int kx = -kmax[0]; kx <= kmax[0]; kx++) {
        double img[3] = {p[0] + kx * period[0], p[1] + ky * period[1], p[2] + kz * period[2]};
        int inside = 1;
        for (int d = 0; d < 3; d++)
          if ((img[d] < lo[d] - cut[d]) || (img[d] > hi[d] + cut[d])) inside = 0;
        if (!inside) continue;
        std::array<int, 3> k = {kx, ky, kz};
        if (std::find(stored.begin(), stored.end(), k) != stored.end()) continue;
        double ximg[3];
        if (triclinic)
          domain->lamda2x(img, ximg);
        else {
          ximg[0] = img[0];
          ximg[1] = img[1];
          ximg[2] = img[2];
        }
        stored_atom(itag, itype, imask, iq, ximg, iv);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   largest distance of a point outside of my subdomain, measured
   perpendicular to each pair of faces (0 inside)
------------------------------------------------------------------------- */

double FixGEMC::subdomain_excess(double *coord)
{
  double excess = 0.0;
  if (triclinic) {
    double lamda[3], w[3];
    domain->x2lamda(coord, lamda);
    box_widths(w);
    for (int d = 0; d < 3; d++) {
      double e = MAX(domain->sublo_lamda[d] - lamda[d], lamda[d] - domain->subhi_lamda[d]);
      excess = MAX(excess, e * w[d]);
    }
  } else {
    for (int d = 0; d < 3; d++)
      excess = MAX(excess, MAX(domain->sublo[d] - coord[d], coord[d] - domain->subhi[d]));
  }
  return excess;
}

/* ----------------------------------------------------------------------
   the stored atoms of a rank include all atoms within its ghost cutoff, so
   they cover all interactions of a point at most the ghost skin outside of
   its subdomain. an atom must also stay within reach of the neighboring
   subdomains, which comm->exchange() can send it to. bad = 1 on the rank
   that owns the atom, if a move would violate either condition. then the
   ghost atoms and the grid are rebuilt, and 1 is returned.
   itag = atom ID of the atom on the rank that owns it and 0 on all other
   ranks, is set on all ranks, since the atom may then be owned by another rank.
   must be called by all ranks of the box.
------------------------------------------------------------------------- */

int FixGEMC::need_refresh(int bad, tagint &itag)
{
  tagint mine[2] = {bad ? 1 : 0, itag};
  tagint all[2];
  MPI_Allreduce(mine, all, 2, MPI_LMP_TAGINT, MPI_MAX, world);
  itag = all[1];
  if (!all[0]) return 0;
  refresh_ghosts();
  build_grid();
  return 1;
}

/* ----------------------------------------------------------------------
   create pending atoms and delete removed atoms in the atom arrays
   must be called by all ranks of the box before the owned atoms are used
   in any other way than through the grid, and at the end of the MC moves
------------------------------------------------------------------------- */

void FixGEMC::flush_pending()
{
  if (!pending_changes) return;

  struct NewAtom {
    tagint tag;
    int type, mask;
    double q, x[3], v[3];
  };
  std::vector<NewAtom> newatoms;
  for (int k : pending) {
    if (k < 0) continue;
    double *xk = atom->x[k];
    double *vk = atom->v[k];
    newatoms.push_back({atom->tag[k],
                        atom->type[k],
                        atom->mask[k],
                        atom->q_flag ? atom->q[k] : 0.0,
                        {xk[0], xk[1], xk[2]},
                        {vk[0], vk[1], vk[2]}});
  }

  // ghost atoms and stored images are invalid from here on

  atom->nghost = 0;
  std::sort(removed.begin(), removed.end(), std::greater<int>());
  for (int i : removed) {
    atom->avec->copy(atom->nlocal - 1, i, 1);
    atom->nlocal--;
  }
  for (auto &n : newatoms) {
    atom->avec->create_atom(n.type, n.x);
    int m = atom->nlocal - 1;
    atom->mask[m] = n.mask;
    atom->tag[m] = n.tag;
    if (atom->q_flag) atom->q[m] = n.q;
    atom->v[m][0] = n.v[0];
    atom->v[m][1] = n.v[1];
    atom->v[m][2] = n.v[2];
    modify->create_attribute(m);
  }

  pending.clear();
  removed.clear();
  pending_changes = 0;
  nstored = 0;
  changed_atoms();
  update_gas_atoms_list();
}

/* ----------------------------------------------------------------------
   compute system potential energy
------------------------------------------------------------------------- */

double FixGEMC::energy_full()
{
  reset_comm();
  int eflag = 1;
  int vflag = 0;

  // clear forces so they don't accumulate over multiple
  // calls within fix gemc timestep, e.g. for fix shake

  size_t nbytes = sizeof(double) * (atom->nlocal + atom->nghost);
  if (nbytes) memset(&atom->f[0][0], 0, 3 * nbytes);

  if (modify->has_fixes(FixConst::PRE_FORCE)) modify->pre_force(vflag);

  if (force->pair) force->pair->compute(eflag, vflag);

  if (atom->molecular != Atom::ATOMIC) {
    if (force->bond) force->bond->compute(eflag, vflag);
    if (force->angle) force->angle->compute(eflag, vflag);
    if (force->dihedral) force->dihedral->compute(eflag, vflag);
    if (force->improper) force->improper->compute(eflag, vflag);
  }

  if (force->kspace) force->kspace->compute(eflag, vflag);

  if (modify->has_fixes(FixConst::POST_FORCE_ANY)) modify->post_force(vflag);

  // NOTE: all fixes with energy_global_flag set and which
  //   operate at pre_force() or post_force()
  //   and which user has enabled via fix_modify energy yes,
  //   will contribute to total MC energy via pe->compute_scalar()

  update->eflag_global = update->ntimestep;
  return c_pe->compute_scalar();
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixGEMC::write_restart(FILE *fp)
{
  int nglobalvals = 14;
  int nlocalvals = comm->nprocs;
  int ntotalvals = nglobalvals + nlocalvals;
  std::vector<double> list(ntotalvals);
  int n = 0;
  list[n++] = -RESTART_VERSION;
  list[n++] = random_world->state();
  list[n++] = random_universe->state();
  list[n++] = ubuf(next_reneighbor).d;
  list[n++] = ntranslation_attempts;
  list[n++] = ntranslation_successes;
  list[n++] = nexchange_attempts;
  list[n++] = nexchange_successes;
  list[n++] = nvolume_attempts;
  list[n++] = nvolume_successes;
  list[n++] = nrotation_attempts;
  list[n++] = nrotation_successes;
  list[n++] = ubuf(update->ntimestep).d;
  list[n++] = nlocalvals;

  double state = random_proc->state();
  MPI_Gather(&state, 1, MPI_DOUBLE, &list[n], 1, MPI_DOUBLE, 0, world);

  if (comm->me == 0) {
    int size = ntotalvals * sizeof(double);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(list.data(), sizeof(double), ntotalvals, fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
   restart files written before the version marker was added start with the
   (positive) random number generator state and contain rotation counters
------------------------------------------------------------------------- */

void FixGEMC::restart(char *buf)
{
  int n = 0;
  auto *list = (double *) buf;

  int oldformat = (list[0] > 0.0);
  int version = oldformat ? 1 : static_cast<int>(-list[0]);
  if (!oldformat) n++;

  // earlier versions stored the state of rank 0 only, ahead of the other values

  int procstate = (version < 4) ? static_cast<int>(list[n++]) : 0;

  random_world->reset(static_cast<int>(list[n++]));
  random_universe->reset(static_cast<int>(list[n++]));

  next_reneighbor = (bigint) ubuf(list[n++]).i;

  ntranslation_attempts = list[n++];
  ntranslation_successes = list[n++];
  if (oldformat) n += 2;
  nexchange_attempts = list[n++];
  nexchange_successes = list[n++];
  nvolume_attempts = list[n++];
  nvolume_successes = list[n++];
  if (version >= 3) {
    nrotation_attempts = list[n++];
    nrotation_successes = list[n++];
  }

  bigint ntimestep_restart = (bigint) ubuf(list[n++]).i;
  if (ntimestep_restart != update->ntimestep)
    error->all(FLERR, "Must not reset timestep when restarting fix gemc");

  // read stored state of RNG unique to this process
  // if nprocs changed, derive a unique state from the state of proc 0

  if (version >= 4) {
    auto nlocalvals = (int) list[n++];
    if (nlocalvals == comm->nprocs) {
      random_proc->reset(static_cast<int>(list[n + me]));
      return;
    }
    procstate = static_cast<int>(list[n]);
  }
  random_proc->reset_rank(procstate, universe->me);
}

/* ----------------------------------------------------------------------
   return cumulative attempt and success counts of this box
------------------------------------------------------------------------- */

double FixGEMC::compute_vector(int n)
{
  if (n == 0) return ntranslation_attempts;
  if (n == 1) return ntranslation_successes;
  if (n == 2) return nexchange_attempts;
  if (n == 3) return nexchange_successes;
  if (n == 4) return nvolume_attempts;
  if (n == 5) return nvolume_successes;
  if (n == 6) return displace;
  if (n == 7) return max_dlogvolratio;
  if (n == 8) return nrotation_attempts;
  if (n == 9) return nrotation_successes;
  if (n == 10) return maxangle / DEG2RAD;
  return 0.0;
}

/* ----------------------------------------------------------------------
   volume of the box (also for triclinic boxes)
------------------------------------------------------------------------- */

double FixGEMC::box_volume()
{
  return domain->xprd * domain->yprd * domain->zprd;
}

/* ----------------------------------------------------------------------
   smallest distance between opposite faces of the box
------------------------------------------------------------------------- */

double FixGEMC::min_box_width()
{
  double w[3];
  box_widths(w);
  return MIN(MIN(w[0], w[1]), w[2]);
}

/* ----------------------------------------------------------------------
   distances between opposite faces of the box
------------------------------------------------------------------------- */

void FixGEMC::box_widths(double *w)
{
  if (!triclinic) {
    w[0] = domain->xprd;
    w[1] = domain->yprd;
    w[2] = domain->zprd;
    return;
  }

  // rows of the inverse box matrix are normal to the box faces,
  // their lengths are the inverse face distances

  double *h_inv = domain->h_inv;
  w[0] = 1.0 / sqrt(h_inv[0] * h_inv[0] + h_inv[5] * h_inv[5] + h_inv[4] * h_inv[4]);
  w[1] = 1.0 / sqrt(h_inv[1] * h_inv[1] + h_inv[3] * h_inv[3]);
  w[2] = 1.0 / h_inv[2];
}

/* ----------------------------------------------------------------------
   largest distance an atom may move in a translation or rotation.
   with more than one rank per box, comm->exchange() in energy_full() can
   only move an atom to a neighboring subdomain, so the distance must stay
   below the smallest subdomain width.  with one rank, domain->pbc() can
   only wrap an atom by one box length, so the distance stays below half
   the smallest box width
------------------------------------------------------------------------- */

double FixGEMC::max_move()
{
  if (comm->nprocs == 1) return 0.5 * min_box_width();
  double w[3];
  box_widths(w);
  double wsub = BIG;
  for (int k = 0; k < 3; k++) {
    double frac = triclinic ? domain->subhi_lamda[k] - domain->sublo_lamda[k]
                            : (domain->subhi[k] - domain->sublo[k]) / domain->prd[k];
    wsub = MIN(wsub, frac * w[k]);
  }
  double wsub_all;
  MPI_Allreduce(&wsub, &wsub_all, 1, MPI_DOUBLE, MPI_MIN, world);
  return 0.5 * wsub_all;
}

/* ----------------------------------------------------------------------
   maximum displacement actually used for translations. the box does not
   change during a translation, so the move remains symmetric
------------------------------------------------------------------------- */

double FixGEMC::max_translation()
{
  return MIN(displace, max_move());
}

/* ----------------------------------------------------------------------
   uniformly distributed random point in my box, identical on all ranks of the box
------------------------------------------------------------------------- */

void FixGEMC::random_point(double *coord)
{
  double lamda[3];
  lamda[0] = random_world->uniform();
  lamda[1] = random_world->uniform();
  lamda[2] = random_world->uniform();
  domain->lamda2x(lamda, coord);
  domain->remap(coord);
}

/* ----------------------------------------------------------------------
   1 if a point inside the box belongs to the subdomain of this rank
------------------------------------------------------------------------- */

int FixGEMC::owns(double *coord)
{
  double *lo, *hi;
  double c[3];
  if (triclinic) {
    domain->x2lamda(coord, c);
    lo = domain->sublo_lamda;
    hi = domain->subhi_lamda;

    // round-off may place a point just outside the unit cube

    for (int k = 0; k < 3; k++) {
      if (c[k] >= 1.0) c[k] -= 1.0;
      if (c[k] < 0.0) c[k] += 1.0;
      if (c[k] >= 1.0) c[k] = 0.0;
    }
  } else {
    c[0] = coord[0];
    c[1] = coord[1];
    c[2] = coord[2];
    lo = domain->sublo;
    hi = domain->subhi;
  }
  return (c[0] >= lo[0]) && (c[0] < hi[0]) && (c[1] >= lo[1]) && (c[1] < hi[1]) &&
      (c[2] >= lo[2]) && (c[2] < hi[2]);
}

/* ----------------------------------------------------------------------
   return local index of the owned atom with atom ID itag or -1 if not owned
------------------------------------------------------------------------- */

int FixGEMC::local_index(tagint itag)
{
  if (atom->map_style != Atom::MAP_NONE) {
    int i = atom->map(itag);
    return (i < atom->nlocal) ? i : -1;
  }
  for (int i = 0; i < atom->nlocal; i++)
    if (atom->tag[i] == itag) return i;
  return -1;
}

/* ----------------------------------------------------------------------
   update after atoms were added or removed or charges changed
   ghost atoms are rebuilt by the next energy_full() or use_local()
------------------------------------------------------------------------- */

void FixGEMC::changed_atoms()
{
  // new or restored atoms were appended after the owned atoms and have
  // overwritten ghost atoms, so the ghost atoms are invalid until rebuilt

  atom->nghost = 0;
  ghosts_stale = 1;
  grid_valid = 0;
  if (atom->map_style != Atom::MAP_NONE) atom->map_init();
  if (force->kspace) force->kspace->qsum_qsq();
  if (force->pair && force->pair->tail_flag) force->pair->reinit();
}

/* ----------------------------------------------------------------------
   set box to lower corner lo (unchanged) with given upper bounds and tilts
------------------------------------------------------------------------- */

void FixGEMC::set_box(double xhi_new, double yhi_new, double zhi_new, double xy_new, double xz_new,
                      double yz_new)
{
  domain->boxhi[0] = xhi_new;
  domain->boxhi[1] = yhi_new;
  domain->boxhi[2] = zhi_new;
  if (triclinic) {
    domain->xy = xy_new;
    domain->xz = xz_new;
    domain->yz = yz_new;
  }
  domain->set_global_box();
  domain->set_local_box();
  comm->setup();
  if (neighbor->style) neighbor->setup_bins();
  if (force->kspace) force->kspace->setup();
  grid_valid = 0;
}

/* ----------------------------------------------------------------------
   scale positions for a change of all box lengths by factor s about the
   lower box corner lo, before the box itself is changed.
   atoms that belong to a molecule are displaced together with the center of
   mass of their molecule, so molecules are not deformed:
     unwrapped x += (s-1) (xcm - lo)
   for the wrapped coordinate x of an atom with unwrapped coordinate xu
   this is x += (s-1) (x + xcm - xu - lo).  atoms without a molecule ID
   are scaled individually, x += (s-1) (x - lo).
   returns the number of independent units (molecules plus single atoms).
   if check is set and the scaled box would be narrower than twice the largest
   molecule, nothing is changed and -1 is returned in both boxes, since bonded
   interactions use the closest periodic image of bond partners.
   must be called at the same time in both boxes if check is set.
------------------------------------------------------------------------- */

bigint FixGEMC::scale_positions(double s, int check)
{
  int nlocal = atom->nlocal;
  double **x = atom->x;
  imageint *image = atom->image;
  tagint *molecule = atom->molecule;
  double *lo = domain->boxlo;
  double sm1 = s - 1.0;

  // centers of mass of all molecules with atoms on this rank

  std::unordered_map<tagint, std::array<double, 3>> com;
  bigint nfree = 0;
  bigint nmolecules = 0;
  if (molecule) {
    nmolecules = molecule_centers(com);
    for (int i = 0; i < nlocal; i++)
      if (molecule[i] <= 0) nfree++;
  } else {
    nfree = nlocal;
  }
  bigint nfree_all;
  MPI_Allreduce(&nfree, &nfree_all, 1, MPI_LMP_BIGINT, MPI_SUM, world);

  // largest distance of an atom from the center of mass of its molecule

  if (check) {
    double rmax = 0.0;
    if (molecule) {
      double xu[3];
      for (int i = 0; i < nlocal; i++) {
        if (molecule[i] <= 0) continue;
        const auto &c = com[molecule[i]];
        domain->unmap(x[i], image[i], xu);
        double dx = xu[0] - c[0];
        double dy = xu[1] - c[1];
        double dz = xu[2] - c[2];
        rmax = MAX(rmax, dx * dx + dy * dy + dz * dz);
      }
    }
    double rmax_all;
    MPI_Allreduce(&rmax, &rmax_all, 1, MPI_DOUBLE, MPI_MAX, world);
    if (any_box(min_box_width() * s <= 4.0 * sqrt(rmax_all))) return -1;
  }

  double xu[3];
  for (int i = 0; i < nlocal; i++) {
    if (molecule && (molecule[i] > 0)) {
      const auto &c = com[molecule[i]];
      domain->unmap(x[i], image[i], xu);
      x[i][0] += sm1 * (x[i][0] + c[0] - xu[0] - lo[0]);
      x[i][1] += sm1 * (x[i][1] + c[1] - xu[1] - lo[1]);
      x[i][2] += sm1 * (x[i][2] + c[2] - xu[2] - lo[2]);
    } else {
      x[i][0] += sm1 * (x[i][0] - lo[0]);
      x[i][1] += sm1 * (x[i][1] - lo[1]);
      x[i][2] += sm1 * (x[i][2] - lo[2]);
    }
  }

  return nfree_all + nmolecules;
}

/* ----------------------------------------------------------------------
   centers of mass (from unwrapped coordinates) of all molecules with atoms
   on this rank. partial sums of each molecule are summed on a rendezvous
   rank (molecule ID modulo number of ranks) and returned to the ranks that
   hold its atoms, so the communication is proportional to the local atoms.
   returns the number of molecules in the box.
   must be called by all ranks of the box.
------------------------------------------------------------------------- */

bigint FixGEMC::molecule_centers(std::unordered_map<tagint, std::array<double, 3>> &com)
{
  int nlocal = atom->nlocal;
  tagint *molecule = atom->molecule;

  // partial sums of mass and mass-weighted unwrapped coordinates

  std::unordered_map<tagint, std::array<double, 4>> partial;
  double xu[3];
  for (int i = 0; i < nlocal; i++) {
    if (molecule[i] <= 0) continue;
    double m = atom->mass[atom->type[i]];
    domain->unmap(atom->x[i], atom->image[i], xu);
    auto &p = partial[molecule[i]];
    p[0] += m;
    p[1] += m * xu[0];
    p[2] += m * xu[1];
    p[3] += m * xu[2];
  }

  int nsend = partial.size();
  int *proclist;
  memory->create(proclist, nsend, "gemc:proclist");
  auto *inbuf = (ComRvous *) memory->smalloc((bigint) nsend * sizeof(ComRvous), "gemc:inbuf");
  int n = 0;
  for (const auto &p : partial) {
    proclist[n] = p.first % comm->nprocs;
    inbuf[n].mol = p.first;
    inbuf[n].proc = comm->me;
    for (int k = 0; k < 4; k++) inbuf[n].sum[k] = p.second[k];
    n++;
  }

  char *buf;
  rvous_count = 0;
  int nreturn = comm->rendezvous(RVOUS, nsend, (char *) inbuf, sizeof(ComRvous), 0, proclist,
                                 rendezvous_centers, 0, buf, sizeof(ComRvous), (void *) this);
  auto *outbuf = (ComRvous *) buf;
  memory->destroy(proclist);
  memory->sfree(inbuf);

  com.clear();
  com.reserve(nreturn);
  for (int k = 0; k < nreturn; k++)
    com[outbuf[k].mol] = {outbuf[k].sum[1], outbuf[k].sum[2], outbuf[k].sum[3]};
  memory->sfree(outbuf);

  bigint nmolecules;
  MPI_Allreduce(&rvous_count, &nmolecules, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  return nmolecules;
}

/* ----------------------------------------------------------------------
   callback of comm->rendezvous() for molecule_centers():
   sum partial sums of each molecule and return its center of mass
   to every rank that sent a partial sum
------------------------------------------------------------------------- */

int FixGEMC::rendezvous_centers(int n, char *inbuf, int &flag, int *&proclist, char *&outbuf,
                                void *ptr)
{
  auto *fptr = (FixGEMC *) ptr;
  Memory *memory = fptr->memory;
  auto *in = (ComRvous *) inbuf;

  std::unordered_map<tagint, std::array<double, 4>> total;
  for (int i = 0; i < n; i++) {
    auto &t = total[in[i].mol];
    for (int k = 0; k < 4; k++) t[k] += in[i].sum[k];
  }
  fptr->rvous_count = total.size();

  memory->create(proclist, n, "gemc:proclist");
  auto *out = (ComRvous *) memory->smalloc((bigint) n * sizeof(ComRvous), "gemc:outbuf");
  for (int i = 0; i < n; i++) {
    const auto &t = total[in[i].mol];
    proclist[i] = in[i].proc;
    out[i].mol = in[i].mol;
    out[i].proc = in[i].proc;
    out[i].sum[0] = t[0];
    out[i].sum[1] = t[1] / t[0];
    out[i].sum[2] = t[2] / t[0];
    out[i].sum[3] = t[3] / t[0];
  }
  outbuf = (char *) out;

  // flag = 2: new outbuf

  flag = 2;
  return n;
}
