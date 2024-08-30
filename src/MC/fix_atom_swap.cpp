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
#include "group.h"
#include "improper.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "random_park.h"
#include "region.h"
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
    qtype(nullptr), sqrt_mass_ratio(nullptr), local_swap_iatom_list(nullptr),
    local_swap_jatom_list(nullptr), local_swap_atom_list(nullptr), random_equal(nullptr),
    random_unequal(nullptr), c_pe(nullptr)
{
  if (narg < 10) error->all(FLERR, "Illegal fix atom/swap command");

  dynamic_group_allow = 1;

  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 0;
  restart_global = 1;
  time_depend = 1;

  // required args

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  ncycles = utils::inumeric(FLERR, arg[4], false, lmp);
  seed = utils::inumeric(FLERR, arg[5], false, lmp);
  double temperature = utils::numeric(FLERR, arg[6], false, lmp);

  if (nevery <= 0) error->all(FLERR, "Illegal fix atom/swap command");
  if (ncycles < 0) error->all(FLERR, "Illegal fix atom/swap command");
  if (seed <= 0) error->all(FLERR, "Illegal fix atom/swap command");
  if (temperature <= 0.0) error->all(FLERR, "Illegal fix atom/swap command");

  beta = 1.0 / (force->boltz * temperature);

  memory->create(type_list, atom->ntypes, "atom/swap:type_list");
  memory->create(mu, atom->ntypes + 1, "atom/swap:mu");
  for (int i = 1; i <= atom->ntypes; i++) mu[i] = 0.0;

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
  memory->destroy(sqrt_mass_ratio);
  memory->destroy(local_swap_iatom_list);
  memory->destroy(local_swap_jatom_list);
  delete[] idregion;
  delete random_equal;
  delete random_unequal;
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

int FixAtomSwap::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAtomSwap::init()
{
  if (!atom->mass) error->all(FLERR, "Fix atom/swap requires per atom type masses");
  if (atom->rmass_flag && (comm->me == 0))
    error->warning(FLERR, "Fix atom/swap will use per-type masses for velocity rescaling");

  c_pe = modify->get_compute_by_id("thermo_pe");

  int *type = atom->type;

  if (nswaptypes < 2) error->all(FLERR, "Must specify at least 2 types in fix atom/swap command");

  if (semi_grand_flag) {
    if (nswaptypes != nmutypes)
      error->all(FLERR, "Need nswaptypes mu values in fix atom/swap command");
  } else {
    if (nswaptypes != 2)
      error->all(FLERR, "Only 2 types allowed when not using semi-grand in fix atom/swap command");
    if (nmutypes != 0)
      error->all(FLERR, "Mu not allowed when not using semi-grand in fix atom/swap command");
  }

  // set index and check validity of region

  if (idregion) {
    region = domain->get_region_by_id(idregion);
    if (!region) error->all(FLERR, "Region {} for fix setforce does not exist", idregion);
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
    }
  }

  memory->create(sqrt_mass_ratio, atom->ntypes + 1, atom->ntypes + 1, "atom/swap:sqrt_mass_ratio");
  for (int itype = 1; itype <= atom->ntypes; itype++)
    for (int jtype = 1; jtype <= atom->ntypes; jtype++)
      sqrt_mass_ratio[itype][jtype] = sqrt(atom->mass[itype] / atom->mass[jtype]);

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
  }
  if (j >= 0) {
    atom->type[j] = itype;
    if (atom->q_flag) atom->q[j] = qtype[0];
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
  }
  if (j >= 0) {
    atom->type[j] = type_list[1];
    if (atom->q_flag) atom->q[j] = qtype[1];
  }

  return 0;
}

/* ----------------------------------------------------------------------
   compute system potential energy
------------------------------------------------------------------------- */

double FixAtomSwap::energy_full()
{
  int eflag = 1;
  int vflag = 0;

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
  double total_energy = c_pe->compute_scalar();

  return total_energy;
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
    memory->sfree(local_swap_atom_list);
    atom_swap_nmax = atom->nmax;
    local_swap_atom_list =
        (int *) memory->smalloc(atom_swap_nmax * sizeof(int), "MCSWAP:local_swap_atom_list");
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
    memory->sfree(local_swap_iatom_list);
    memory->sfree(local_swap_jatom_list);
    atom_swap_nmax = atom->nmax;
    local_swap_iatom_list =
        (int *) memory->smalloc(atom_swap_nmax * sizeof(int), "MCSWAP:local_swap_iatom_list");
    local_swap_jatom_list =
        (int *) memory->smalloc(atom_swap_nmax * sizeof(int), "MCSWAP:local_swap_jatom_list");
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
  auto list = (double *) buf;

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
