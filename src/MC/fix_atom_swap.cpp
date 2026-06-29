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
   Contributing authors: Paul Crozier (SNL)
                         Alexander Stukowski
------------------------------------------------------------------------- */

#include "fix_atom_swap.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "compute.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "graphics.h"
#include "group.h"
#include "improper.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "random_park.h"
#include "region.h"
#include "suffix.h"
#include "update.h"

#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAtomSwap::FixAtomSwap(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), region(nullptr), idregion(nullptr), type_list(nullptr), mu(nullptr),
    qtype(nullptr), mtype(nullptr), sqrt_mass_ratio(nullptr), local_swap_iatom_list(nullptr),
    local_swap_jatom_list(nullptr), local_swap_atom_list(nullptr), random_equal(nullptr),
    random_unequal(nullptr), c_pe(nullptr), imgobjs(nullptr), imgparms(nullptr)
{
  if (narg < 10) utils::missing_cmd_args(FLERR, "fix atom/swap", error);

  dynamic_group_allow = 1;

  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 0;
  restart_global = 1;
  time_depend = 1;

  // no visualization without an atom map
  if (atom->map_style == Atom::MAP_NONE) {
    vizsteps = 0;
  } else {
    vizsteps = 1000;
  }

  // required args

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  ncycles = utils::inumeric(FLERR, arg[4], false, lmp);
  seed = utils::inumeric(FLERR, arg[5], false, lmp);
  double temperature = utils::numeric(FLERR, arg[6], false, lmp);

  if (nevery <= 0) error->all(FLERR, 3, "Illegal fix atom/swap command nevery value");
  if (ncycles < 0) error->all(FLERR, 4, "Illegal fix atom/swap command ncycles value");
  if (seed <= 0) error->all(FLERR, 5, "Illegal fix atom/swap command random seed");
  if (temperature <= 0.0) error->all(FLERR, 6, "Illegal fix atom/swap command temperature value");

  beta = 1.0 / (force->boltz * temperature);

  memory->create(type_list, atom->ntypes, "atom/swap:type_list");
  memory->create(mu, atom->ntypes + 1, "atom/swap:mu");
  for (int i = 0; i <= atom->ntypes; i++) mu[i] = 0.0;

  // read options from end of input line

  options(narg - 7, &arg[7]);

  // random number generator, same for all procs

  random_equal = new RanPark(lmp, seed);

  // random number generator, not the same for all procs

  random_unequal = new RanPark(lmp, seed);

  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;

  // zero out counters

  mc_active = 0;

  nswap_attempts = 0.0;
  nswap_successes = 0.0;

  atom_swap_nmax = 0;
  local_swap_atom_list = nullptr;
  local_swap_iatom_list = nullptr;
  local_swap_jatom_list = nullptr;

  // set comm size needed by this Fix

  if (atom->q_flag)
    comm_forward = 2;
  else
    comm_forward = 1;
}

/* ---------------------------------------------------------------------- */

FixAtomSwap::~FixAtomSwap()
{
  memory->destroy(type_list);
  memory->destroy(mu);
  memory->destroy(qtype);
  memory->destroy(mtype);
  memory->destroy(sqrt_mass_ratio);
  memory->destroy(local_swap_iatom_list);
  memory->destroy(local_swap_jatom_list);
  delete[] idregion;
  delete random_equal;
  delete random_unequal;
  memory->destroy(imgobjs);
  memory->destroy(imgparms);
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void FixAtomSwap::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR, "Illegal fix atom/swap command");

  ke_flag = 1;
  semi_grand_flag = 0;
  nswaptypes = 0;
  nmutypes = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "region") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix atom/swap command");
      region = domain->get_region_by_id(arg[iarg + 1]);
      if (!region) error->all(FLERR, "Region {} for fix atom/swap does not exist", arg[iarg + 1]);
      idregion = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "ke") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix atom/swap command");
      ke_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "semi-grand") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix atom/swap command");
      semi_grand_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "types") == 0) {
      if (iarg + 3 > narg) error->all(FLERR, "Illegal fix atom/swap command");
      iarg++;
      while (iarg < narg) {
        if (isalpha(arg[iarg][0])) break;
        if (nswaptypes >= atom->ntypes) error->all(FLERR, "Illegal fix atom/swap command");
        type_list[nswaptypes] = utils::expand_type_int(FLERR, arg[iarg], Atom::ATOM, lmp);
        nswaptypes++;
        iarg++;
      }
    } else if (strcmp(arg[iarg], "mu") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix atom/swap command");
      iarg++;
      while (iarg < narg) {
        if (isalpha(arg[iarg][0])) break;
        nmutypes++;
        if (nmutypes > atom->ntypes) error->all(FLERR, "Illegal fix atom/swap command");
        mu[nmutypes] = utils::numeric(FLERR, arg[iarg], false, lmp);
        iarg++;
      }
    } else
      error->all(FLERR, "Illegal fix atom/swap command");
  }
}

/* ---------------------------------------------------------------------- */

int FixAtomSwap::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"vizsteps") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "fix_modify atom/swap", error);
    vizsteps = utils::inumeric(FLERR, arg[1], false, lmp);
    return 2;
  }

  return 0;
}

/* ---------------------------------------------------------------------- */

int FixAtomSwap::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAtomSwap::init()
{
  if ((atom->mass != nullptr) && (atom->rmass != nullptr) && (comm->me == 0))
    error->warning(FLERR, "Fix atom/swap will use per-atom masses for velocity rescaling");

  c_pe = modify->get_compute_by_id("thermo_pe");

  int *type = atom->type;

  if (nswaptypes < 2)
    error->all(FLERR, Error::NOLASTLINE,
               "Must specify at least 2 atom types in fix atom/swap command");

  if (semi_grand_flag) {
    if (nswaptypes != nmutypes)
      error->all(FLERR, Error::NOLASTLINE, "Need nswaptypes mu values in fix atom/swap command");
  } else {
    if (nswaptypes != 2)
      error->all(FLERR, Error::NOLASTLINE,
                 "Exactly 2 atom types must be used without semi-grand keyword in fix atom/swap");
    if (nmutypes != 0)
      error->all(FLERR, Error::NOLASTLINE,
                 "Mu not allowed when not using semi-grand in fix atom/swap command");
  }

  // must have a pair style and not use INTEL package

  if (!force->pair) error->all(FLERR, Error::NOLASTLINE, "Fix atom/swap requires a pair style");
  if (force->pair && (force->pair->suffix_flag & Suffix::INTEL))
    error->all(FLERR, Error::NOLASTLINE, "Fix {} is not compatible with /intel pair styles", style);

  // check if constraints for hybrid pair styles are fulfilled

  if (utils::strmatch(force->pair_style, "^hybrid")) {
    auto *hybrid = dynamic_cast<PairHybrid *>(force->pair);
    if (hybrid) {
      for (int i = 0; i < nswaptypes - 1; ++i) {
        int type1 = type_list[i];
        for (int j = i + 1; j < nswaptypes; ++j) {
          int type2 = type_list[j];
          if (hybrid->nmap[type1][type1] != hybrid->nmap[type2][type2])
            error->all(FLERR, Error::NOLASTLINE,
                       "Pair {} substyles for atom types {} and {} are not compatible",
                       force->pair_style, type1, type2);
          for (int k = 0; k < hybrid->nmap[type1][type1]; ++k) {
            if (hybrid->map[type1][type1][k] != hybrid->map[type2][type2][k])
              error->all(FLERR, Error::NOLASTLINE,
                         "Pair {} substyles for atom types {} and {} are not compatible",
                         force->pair_style, type1, type2);
          }
        }
      }
    }
  }

  // set index and check validity of region

  if (idregion) {
    region = domain->get_region_by_id(idregion);
    if (!region)
      error->all(FLERR, Error::NOLASTLINE, "Region {} for fix atom/swap does not exist", idregion);
  }

  for (int iswaptype = 0; iswaptype < nswaptypes; iswaptype++)
    if (type_list[iswaptype] <= 0 || type_list[iswaptype] > atom->ntypes)
      error->all(FLERR, "Invalid atom type in fix atom/swap command");

  // this is only required for non-semi-grand
  // in which case, nswaptypes = 2

  if (atom->q_flag && !semi_grand_flag) {
    double qmax, qmin;
    int firstall, first;
    memory->create(qtype, nswaptypes, "atom/swap:qtype");
    for (int iswaptype = 0; iswaptype < nswaptypes; iswaptype++) {
      first = 1;
      for (int i = 0; i < atom->nlocal; i++) {
        if (atom->mask[i] & groupbit) {
          if (type[i] == type_list[iswaptype]) {
            if (first) {
              qtype[iswaptype] = atom->q[i];
              first = 0;
            } else if (qtype[iswaptype] != atom->q[i])
              error->one(FLERR, "All atoms of a swapped type must have the same charge.");
          }
        }
      }
      MPI_Allreduce(&first, &firstall, 1, MPI_INT, MPI_MIN, world);
      if (firstall)
        error->all(FLERR,
                   "At least one atom of each swapped type must be present to define charges.");
      if (first) qtype[iswaptype] = -DBL_MAX;
      MPI_Allreduce(&qtype[iswaptype], &qmax, 1, MPI_DOUBLE, MPI_MAX, world);
      if (first) qtype[iswaptype] = DBL_MAX;
      MPI_Allreduce(&qtype[iswaptype], &qmin, 1, MPI_DOUBLE, MPI_MIN, world);
      if (qmax != qmin) error->all(FLERR, "All atoms of a swapped type must have same charge.");
      qtype[iswaptype] = qmax;
    }
  }

  // if we have per-atom masses, check that rmass is consistent with type,
  // and set per-type mass to that value
  if ((atom->rmass !=  nullptr) && !semi_grand_flag) {
    double mmax, mmin;
    int firstall, first;
    memory->create(mtype, nswaptypes, "atom/swap:mtype");
    for (int iswaptype = 0; iswaptype < nswaptypes; iswaptype++) {
      first = 1;
      for (int i = 0; i < atom->nlocal; i++) {
        if (atom->mask[i] & groupbit) {
          if (type[i] == type_list[iswaptype]) {
            if (first > 0) {
              mtype[iswaptype] = atom->rmass[i];
              first = 0;
            } else if (mtype[iswaptype] != atom->rmass[i])
              first = -1;
          }
        }
      }
      MPI_Allreduce(&first, &firstall, 1, MPI_INT, MPI_MIN, world);
      if (firstall < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "All atoms of a swapped type must have the same per-atom mass");
      if (firstall > 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "At least one atom of each swapped type must be present to define masses");
      if (first) mtype[iswaptype] = -DBL_MAX;
      MPI_Allreduce(&mtype[iswaptype], &mmax, 1, MPI_DOUBLE, MPI_MAX, world);
      if (first) mtype[iswaptype] = DBL_MAX;
      MPI_Allreduce(&mtype[iswaptype], &mmin, 1, MPI_DOUBLE, MPI_MIN, world);
      if (mmax != mmin)
        error->all(FLERR, Error::NOLASTLINE, "All atoms of a swapped type must have same mass.");
      mtype[iswaptype] = mmax;
    }
  }

  memory->create(sqrt_mass_ratio, atom->ntypes + 1, atom->ntypes + 1, "atom/swap:sqrt_mass_ratio");
  if (atom->rmass != nullptr) {
    for (int itype = 1; itype <= atom->ntypes; itype++)
      for (int jtype = 1; jtype <= atom->ntypes; jtype++) sqrt_mass_ratio[itype][jtype] = 1.0;
    for (int iswaptype = 0; iswaptype < nswaptypes; iswaptype++) {
      int itype = type_list[iswaptype];
      for (int jswaptype = 0; jswaptype < nswaptypes; jswaptype++) {
        int jtype = type_list[jswaptype];
        sqrt_mass_ratio[itype][jtype] = sqrt(mtype[iswaptype] / mtype[jswaptype]);
      }
    }
  } else {
    for (int itype = 1; itype <= atom->ntypes; itype++)
      for (int jtype = 1; jtype <= atom->ntypes; jtype++)
        sqrt_mass_ratio[itype][jtype] = sqrt(atom->mass[itype] / atom->mass[jtype]);
  }

  // check to see if itype and jtype cutoffs are the same
  // if not, reneighboring will be needed between swaps

  double **cutsq = force->pair->cutsq;
  unequal_cutoffs = false;
  for (int iswaptype = 0; iswaptype < nswaptypes; iswaptype++)
    for (int jswaptype = 0; jswaptype < nswaptypes; jswaptype++)
      for (int ktype = 1; ktype <= atom->ntypes; ktype++)
        if (cutsq[type_list[iswaptype]][ktype] != cutsq[type_list[jswaptype]][ktype])
          unequal_cutoffs = true;

  // check that no swappable atoms are in atom->firstgroup
  // swapping such an atom might not leave firstgroup atoms first

  if (atom->firstgroup >= 0) {
    int *mask = atom->mask;
    int firstgroupbit = group->bitmask[atom->firstgroup];

    int flag = 0;
    for (int i = 0; i < atom->nlocal; i++)
      if ((mask[i] == groupbit) && (mask[i] && firstgroupbit)) flag = 1;

    int flagall;
    MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);

    if (flagall) error->all(FLERR, "Cannot do atom/swap on atoms in atom_modify first group");
  }
}

/* ----------------------------------------------------------------------
   attempt Monte Carlo swaps
------------------------------------------------------------------------- */

void FixAtomSwap::pre_exchange()
{
  // just return if should not be called on this timestep

  if (next_reneighbor != update->ntimestep) return;

  mc_active = 1;

  // ensure current system is ready to compute energy

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
  if (modify->n_pre_neighbor) modify->pre_neighbor();
  neighbor->build(1);

  // energy_stored = energy of current state
  // will be updated after accepted swaps

  energy_stored = energy_full();

  // attempt Ncycle atom swaps

  int nsuccess = 0;
  if (semi_grand_flag) {
    update_semi_grand_atoms_list();
    for (int i = 0; i < ncycles; i++) nsuccess += attempt_semi_grand();
  } else {
    update_swap_atoms_list();
    for (int i = 0; i < ncycles; i++) nsuccess += attempt_swap();
  }

  // udpate MC stats

  nswap_attempts += ncycles;
  nswap_successes += nsuccess;

  next_reneighbor = update->ntimestep + nevery;

  mc_active = 0;

  // if visualization support is enabled, age vizatoms and remove expired ones
  if (vizsteps > 0) {
    std::vector<tagint> eraseme;
    for (const auto &[key, data] : vizatoms) {
      int idx = atom->map(key);
      if ((idx < 0) || (data.first < 0)) {
        eraseme.push_back(key);
        continue;
      }
      vizatoms[key] = std::make_pair(data.first - nevery, data.second);
    }
    for (const auto &key : eraseme) vizatoms.erase(key);
  }
}

/* ----------------------------------------------------------------------
   attempt a semd-grand swap of a single atom
   compare before/after energy and accept/reject the swap
   NOTE: atom charges are assumed equal and so are not updated
------------------------------------------------------------------------- */

int FixAtomSwap::attempt_semi_grand()
{
  if (nswap == 0) return 0;

  // pre-swap energy

  double energy_before = energy_stored;

  // pick a random atom and perform swap

  int itype, jtype, jswaptype;
  int i = pick_semi_grand_atom();
  if (i >= 0) {
    jswaptype = static_cast<int>(nswaptypes * random_unequal->uniform());
    jtype = type_list[jswaptype];
    itype = atom->type[i];
    while (itype == jtype) {
      jswaptype = static_cast<int>(nswaptypes * random_unequal->uniform());
      jtype = type_list[jswaptype];
    }
    atom->type[i] = jtype;
  }

  // if unequal_cutoffs, call comm->borders() and rebuild neighbor list
  // else communicate ghost atoms
  // call to comm->exchange() is a no-op but clears ghost atoms

  if (unequal_cutoffs) {
    if (domain->triclinic) domain->x2lamda(atom->nlocal);
    comm->exchange();
    comm->borders();
    if (domain->triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
    if (modify->n_pre_neighbor) modify->pre_neighbor();
    neighbor->build(1);
  } else {
    comm->forward_comm(this);
  }

  // post-swap energy

  if (force->kspace) force->kspace->qsum_qsq();
  double energy_after = energy_full();

  int success = 0;
  if (i >= 0)
    if (random_unequal->uniform() <
        exp(beta * (energy_before - energy_after + mu[jtype] - mu[itype])))
      success = 1;

  int success_all = 0;
  MPI_Allreduce(&success, &success_all, 1, MPI_INT, MPI_MAX, world);

  // swap accepted, return 1

  if (success_all) {
    update_semi_grand_atoms_list();
    energy_stored = energy_after;
    if (ke_flag) {
      if (i >= 0) {
        atom->v[i][0] *= sqrt_mass_ratio[itype][jtype];
        atom->v[i][1] *= sqrt_mass_ratio[itype][jtype];
        atom->v[i][2] *= sqrt_mass_ratio[itype][jtype];
        // record atom for which the type was swapped and store the old type
        if (vizsteps > 0) {
          vizatoms[atom->tag[i]] = std::make_pair(vizsteps,itype);
        }
      }
    }
    return 1;
  }

  // swap not accepted, return 0
  // restore the swapped atom
  // do not need to re-call comm->borders() and rebuild neighbor list
  //   since will be done on next cycle or in Verlet when this fix finishes

  if (i >= 0) atom->type[i] = itype;
  if (force->kspace) force->kspace->qsum_qsq();

  return 0;
}

/* ----------------------------------------------------------------------
   attempt a swap of a pair of atoms
   compare before/after energy and accept/reject the swap
------------------------------------------------------------------------- */

int FixAtomSwap::attempt_swap()
{
  if ((niswap == 0) || (njswap == 0)) return 0;

  // pre-swap energy

  double energy_before = energy_stored;

  // pick a random pair of atoms
  // swap their properties

  int i = pick_i_swap_atom();
  int j = pick_j_swap_atom();
  int itype = type_list[0];
  int jtype = type_list[1];

  if (i >= 0) {
    atom->type[i] = jtype;
    if (atom->q_flag) atom->q[i] = qtype[1];
    if (atom->rmass != nullptr) atom->rmass[i] = mtype[1];
  }
  if (j >= 0) {
    atom->type[j] = itype;
    if (atom->q_flag) atom->q[j] = qtype[0];
    if (atom->rmass != nullptr) atom->rmass[j] = mtype[0];
  }

  // if unequal_cutoffs, call comm->borders() and rebuild neighbor list
  // else communicate ghost atoms
  // call to comm->exchange() is a no-op but clears ghost atoms

  if (unequal_cutoffs) {
    if (domain->triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    comm->exchange();
    comm->borders();
    if (domain->triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
    if (modify->n_pre_neighbor) modify->pre_neighbor();
    neighbor->build(1);
  } else {
    comm->forward_comm(this);
  }

  // post-swap energy

  double energy_after = energy_full();

  // swap accepted, return 1
  // if ke_flag, rescale atom velocities

  if (random_equal->uniform() < exp(beta * (energy_before - energy_after))) {
    update_swap_atoms_list();
    if (ke_flag) {
      if (i >= 0) {
        atom->v[i][0] *= sqrt_mass_ratio[itype][jtype];
        atom->v[i][1] *= sqrt_mass_ratio[itype][jtype];
        atom->v[i][2] *= sqrt_mass_ratio[itype][jtype];
      }
      if (j >= 0) {
        atom->v[j][0] *= sqrt_mass_ratio[jtype][itype];
        atom->v[j][1] *= sqrt_mass_ratio[jtype][itype];
        atom->v[j][2] *= sqrt_mass_ratio[jtype][itype];
      }
      // record atoms for which the type was swapped and store the old types
      if (vizsteps > 0) {
        vizatoms[atom->tag[i]] = std::make_pair(vizsteps, jtype);
        vizatoms[atom->tag[j]] = std::make_pair(vizsteps, itype);
      }
    }
    energy_stored = energy_after;
    return 1;
  }

  // swap not accepted, return 0
  // restore the swapped itype & jtype atoms
  // do not need to re-call comm->borders() and rebuild neighbor list
  //   since will be done on next cycle or in Verlet when this fix finishes

  if (i >= 0) {
    atom->type[i] = type_list[0];
    if (atom->q_flag) atom->q[i] = qtype[0];
    if (atom->rmass != nullptr) atom->rmass[i] = mtype[0];
  }
  if (j >= 0) {
    atom->type[j] = type_list[1];
    if (atom->q_flag) atom->q[j] = qtype[1];
    if (atom->rmass != nullptr) atom->rmass[j] = mtype[1];
  }

  return 0;
}

/* ----------------------------------------------------------------------
   compute system potential energy
------------------------------------------------------------------------- */

double FixAtomSwap::energy_full()
{
  // flag that we only need to compute the global energy
  int eflag = ENERGY_GLOBAL | ENERGY_ONLY;
  int vflag = VIRIAL_NONE;

  if (modify->n_pre_force) modify->pre_force(vflag);

  if (force->pair) force->pair->compute(eflag, vflag);

  if (atom->molecular != Atom::ATOMIC) {
    if (force->bond) force->bond->compute(eflag, vflag);
    if (force->angle) force->angle->compute(eflag, vflag);
    if (force->dihedral) force->dihedral->compute(eflag, vflag);
    if (force->improper) force->improper->compute(eflag, vflag);
  }

  if (force->kspace) force->kspace->compute(eflag, vflag);

  if (modify->n_post_force_any) modify->post_force(vflag);

  update->eflag_global = update->ntimestep;
  return c_pe->compute_scalar();
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int FixAtomSwap::pick_semi_grand_atom()
{
  int i = -1;
  int iwhichglobal = static_cast<int>(nswap * random_equal->uniform());
  if ((iwhichglobal >= nswap_before) && (iwhichglobal < nswap_before + nswap_local)) {
    int iwhichlocal = iwhichglobal - nswap_before;
    i = local_swap_atom_list[iwhichlocal];
  }

  return i;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int FixAtomSwap::pick_i_swap_atom()
{
  int i = -1;
  int iwhichglobal = static_cast<int>(niswap * random_equal->uniform());
  if ((iwhichglobal >= niswap_before) && (iwhichglobal < niswap_before + niswap_local)) {
    int iwhichlocal = iwhichglobal - niswap_before;
    i = local_swap_iatom_list[iwhichlocal];
  }

  return i;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int FixAtomSwap::pick_j_swap_atom()
{
  int j = -1;
  int jwhichglobal = static_cast<int>(njswap * random_equal->uniform());
  if ((jwhichglobal >= njswap_before) && (jwhichglobal < njswap_before + njswap_local)) {
    int jwhichlocal = jwhichglobal - njswap_before;
    j = local_swap_jatom_list[jwhichlocal];
  }

  return j;
}

/* ----------------------------------------------------------------------
   update the list of gas atoms
------------------------------------------------------------------------- */

void FixAtomSwap::update_semi_grand_atoms_list()
{
  int nlocal = atom->nlocal;
  double **x = atom->x;

  if (atom->nmax > atom_swap_nmax) {
    memory->destroy(local_swap_atom_list);
    atom_swap_nmax = atom->nmax;
    memory->create(local_swap_atom_list, atom_swap_nmax, "MCSWAP:local_swap_atom_list");
  }

  nswap_local = 0;

  if (region) {
    for (int i = 0; i < nlocal; i++) {
      if (region->match(x[i][0], x[i][1], x[i][2]) == 1) {
        if (atom->mask[i] & groupbit) {
          int itype = atom->type[i];
          int iswaptype;
          for (iswaptype = 0; iswaptype < nswaptypes; iswaptype++)
            if (itype == type_list[iswaptype]) break;
          if (iswaptype == nswaptypes) continue;
          local_swap_atom_list[nswap_local] = i;
          nswap_local++;
        }
      }
    }

  } else {
    for (int i = 0; i < nlocal; i++) {
      if (atom->mask[i] & groupbit) {
        int itype = atom->type[i];
        int iswaptype;
        for (iswaptype = 0; iswaptype < nswaptypes; iswaptype++)
          if (itype == type_list[iswaptype]) break;
        if (iswaptype == nswaptypes) continue;
        local_swap_atom_list[nswap_local] = i;
        nswap_local++;
      }
    }
  }

  MPI_Allreduce(&nswap_local, &nswap, 1, MPI_INT, MPI_SUM, world);
  MPI_Scan(&nswap_local, &nswap_before, 1, MPI_INT, MPI_SUM, world);
  nswap_before -= nswap_local;
}

/* ----------------------------------------------------------------------
   update the list of gas atoms
------------------------------------------------------------------------- */

void FixAtomSwap::update_swap_atoms_list()
{
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double **x = atom->x;

  if (atom->nmax > atom_swap_nmax) {
    memory->destroy(local_swap_iatom_list);
    memory->destroy(local_swap_jatom_list);
    atom_swap_nmax = atom->nmax;
    memory->create(local_swap_iatom_list, atom_swap_nmax, "MCSWAP:local_swap_iatom_list");
    memory->create(local_swap_jatom_list, atom_swap_nmax, "MCSWAP:local_swap_jatom_list");
  }

  niswap_local = 0;
  njswap_local = 0;

  if (region) {

    for (int i = 0; i < nlocal; i++) {
      if (region->match(x[i][0], x[i][1], x[i][2]) == 1) {
        if (atom->mask[i] & groupbit) {
          if (type[i] == type_list[0]) {
            local_swap_iatom_list[niswap_local] = i;
            niswap_local++;
          } else if (type[i] == type_list[1]) {
            local_swap_jatom_list[njswap_local] = i;
            njswap_local++;
          }
        }
      }
    }

  } else {
    for (int i = 0; i < nlocal; i++) {
      if (atom->mask[i] & groupbit) {
        if (type[i] == type_list[0]) {
          local_swap_iatom_list[niswap_local] = i;
          niswap_local++;
        } else if (type[i] == type_list[1]) {
          local_swap_jatom_list[njswap_local] = i;
          njswap_local++;
        }
      }
    }
  }

  MPI_Allreduce(&niswap_local, &niswap, 1, MPI_INT, MPI_SUM, world);
  MPI_Scan(&niswap_local, &niswap_before, 1, MPI_INT, MPI_SUM, world);
  niswap_before -= niswap_local;

  MPI_Allreduce(&njswap_local, &njswap, 1, MPI_INT, MPI_SUM, world);
  MPI_Scan(&njswap_local, &njswap_before, 1, MPI_INT, MPI_SUM, world);
  njswap_before -= njswap_local;
}

/* ---------------------------------------------------------------------- */

int FixAtomSwap::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;

  int *type = atom->type;
  double *q = atom->q;

  m = 0;

  if (atom->q_flag) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = type[j];
      buf[m++] = q[j];
    }
  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = type[j];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixAtomSwap::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  int *type = atom->type;
  double *q = atom->q;

  m = 0;
  last = first + n;

  if (atom->q_flag) {
    for (i = first; i < last; i++) {
      type[i] = static_cast<int>(buf[m++]);
      q[i] = buf[m++];
    }
  } else {
    for (i = first; i < last; i++) type[i] = static_cast<int>(buf[m++]);
  }
}

/* ----------------------------------------------------------------------
  return acceptance ratio
------------------------------------------------------------------------- */

double FixAtomSwap::compute_vector(int n)
{
  if (n == 0) return nswap_attempts;
  if (n == 1) return nswap_successes;
  return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixAtomSwap::memory_usage()
{
  double bytes = (double) atom_swap_nmax * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixAtomSwap::write_restart(FILE *fp)
{
  int n = 0;
  double list[6];
  list[n++] = random_equal->state();
  list[n++] = random_unequal->state();
  list[n++] = ubuf(next_reneighbor).d;
  list[n++] = nswap_attempts;
  list[n++] = nswap_successes;
  list[n++] = ubuf(update->ntimestep).d;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(list, sizeof(double), n, fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixAtomSwap::restart(char *buf)
{
  int n = 0;
  auto *list = (double *) buf;

  seed = static_cast<int>(list[n++]);
  random_equal->reset(seed);

  seed = static_cast<int>(list[n++]);
  random_unequal->reset(seed);

  next_reneighbor = (bigint) ubuf(list[n++]).i;

  nswap_attempts = static_cast<int>(list[n++]);
  nswap_successes = static_cast<int>(list[n++]);

  bigint ntimestep_restart = (bigint) ubuf(list[n++]).i;
  if (ntimestep_restart != update->ntimestep)
    error->all(FLERR, "Must not reset timestep when restarting fix atom/swap");
}

/* ----------------------------------------------------------------------
   extract variable which stores whether MC is active or not
     active = MC moves are taking place
     not active = normal MD is taking place
------------------------------------------------------------------------- */

void *FixAtomSwap::extract(const char *name, int &dim)
{
  if (strcmp(name,"mc_active") == 0) {
    dim = 0;
    return (void *) &mc_active;
  }
  return nullptr;
}

/* ----------------------------------------------------------------------
   provide graphics information to dump image to render spheres
   at the location of atoms that were involved in a reaction
------------------------------------------------------------------------- */

int FixAtomSwap::image(int *&objs, double **&parms)
{
  // no visualization without an atom map
  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, Error::NOLASTLINE,
               "Cannot use fix atom/swap in dump image without an atom map");

  memory->destroy(imgobjs);
  memory->destroy(imgparms);

  int numobjs = vizatoms.size();
  int n = 0;
  if (numobjs > 0) {
    memory->create(imgobjs, numobjs, "atom/swap:imgobjs");
    memory->create(imgparms, numobjs, 5, "atom/swap:imgparms");

    int idx;
    const auto *const *const x = atom->x;
    for (const auto &[key, data] : vizatoms) {
      idx = atom->map(key);
      if (idx < 0) continue;
      imgobjs[n] = Graphics::SPHERE;
      imgparms[n][0] = data.second; // use stored pre-swap atom type
      imgparms[n][1] = x[idx][0];
      imgparms[n][2] = x[idx][1];
      imgparms[n][3] = x[idx][2];
      imgparms[n][4] = 0.0;     // radius is set with fflag2 in dump image
      ++n;
    }
  }
  objs = imgobjs;
  parms = imgparms;
  return n;
}
