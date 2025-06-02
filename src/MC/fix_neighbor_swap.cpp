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
   Contributing authors: Jacob Tavenner
------------------------------------------------------------------------- */

#include "fix_neighbor_swap.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "compute.h"
#include "compute_voronoi_atom.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "improper.h"
#include "kspace.h"
#include "math_special.h"
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

static const char cite_fix_neighbor_swap_c[] =
    "fix neighbor/swap command: doi:10.1016/j.commatsci.2022.111929\n\n"
    "@Article{Tavenner2023111929,\n"
    " author = {Jacob P. Tavenner and Mikhail I. Mendelev and John W. Lawson},\n"
    " title = {Molecular dynamics based kinetic Monte Carlo simulation for accelerated "
    "diffusion},\n"
    " journal = {Computational Materials Science},\n"
    " year = {2023},\n"
    " volume = {218},\n"
    " pages = {111929}\n"
    " url = {https://www.sciencedirect.com/science/article/pii/S0927025622006401}\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

FixNeighborSwap::FixNeighborSwap(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), region(nullptr), idregion(nullptr), type_list(nullptr), qtype(nullptr),
    c_voro(nullptr), voro_neighbor_list(nullptr), sqrt_mass_ratio(nullptr),
    local_swap_iatom_list(nullptr), random_equal(nullptr), c_pe(nullptr)
{
  if (narg < 10) utils::missing_cmd_args(FLERR, "fix neighbor/swap", error);

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
  r_0 = utils::inumeric(FLERR, arg[7], false, lmp);

  // Voro compute check

  int icompute = modify->find_compute(utils::strdup(arg[8]));
  if (icompute < 0) error->all(FLERR, "Could not find neighbor compute ID");
  c_voro = modify->compute[icompute];
  if (c_voro->local_flag == 0)
    error->all(FLERR, "Neighbor compute does not compute local info");
  if (c_voro->size_local_cols != 3)
    error->all(FLERR, "Neighbor compute does not give i, j, size as expected");

  if (nevery <= 0) error->all(FLERR, "Illegal fix neighbor/swap command nevery value");
  if (ncycles < 0) error->all(FLERR, "Illegal fix neighbor/swap command ncycles value");
  if (seed <= 0) error->all(FLERR, "Illegal fix neighbor/swap command seed value");
  if (temperature <= 0.0) error->all(FLERR, "Illegal fix neighbor/swap command temperature value");

  beta = 1.0 / (force->boltz * temperature);

  memory->create(type_list, atom->ntypes, "neighbor/swap:type_list");
  memory->create(rate_list, atom->ntypes, "neighbor/swap:rate_list");

  // read options from end of input line

  options(narg - 8, &arg[8]);

  // random number generator, same for all procs

  random_equal = new RanPark(lmp, seed);

  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;

  // zero out counters

  nswap_attempts = 0.0;
  nswap_successes = 0.0;

  atom_swap_nmax = 0;
  voro_neighbor_list = nullptr;
  local_swap_iatom_list = nullptr;
  local_swap_neighbor_list = nullptr;
  local_swap_probability = nullptr;
  local_swap_type_list = nullptr;

  // set comm size needed by this Fix

  if (atom->q_flag)
    comm_forward = 2;
  else
    comm_forward = 1;
}

/* ---------------------------------------------------------------------- */

FixNeighborSwap::~FixNeighborSwap()
{
  memory->destroy(type_list);
  memory->destroy(rate_list);
  memory->destroy(qtype);
  memory->destroy(sqrt_mass_ratio);
  memory->destroy(local_swap_iatom_list);
  memory->destroy(local_swap_neighbor_list);
  memory->destroy(local_swap_probability);
  memory->destroy(local_swap_type_list);
  delete[] idregion;
  delete random_equal;
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void FixNeighborSwap::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR, "Illegal fix neighbor/swap command\n");

  ke_flag = 1;
  diff_flag = 0;
  rates_flag = 0;
  nswaptypes = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "region") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix neighbor/swap command");
      region = domain->get_region_by_id(arg[iarg + 1]);
      if (!region) error->all(FLERR, "Region ID for fix neighbor/swap does not exist");
      idregion = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "ke") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix neighbor/swap command");
      ke_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "types") == 0) {
      if (iarg + 3 > narg) error->all(FLERR, "Illegal fix neighbor/swap command");
      if (diff_flag != 0) error->all(FLERR, "Illegal fix neighbor/swap command");
      iarg++;
      nswaptypes = 0;
      while (iarg < narg) {
        if (isalpha(arg[iarg][0])) break;
        if (nswaptypes >= atom->ntypes) error->all(FLERR, "Illegal fix neighbor/swap command");
        type_list[nswaptypes] = utils::numeric(FLERR, arg[iarg], false, lmp);
        nswaptypes++;
        iarg++;
      }
    } else if (strcmp(arg[iarg], "diff") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix neighbor/swap command");
      if (nswaptypes != 0) error->all(FLERR, "Illegal fix neighbor/swap command");
      type_list[nswaptypes] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      diff_flag = 1;
      nswaptypes++;
      iarg += 2;
    } else if (strcmp(arg[iarg], "rates") == 0) {
      if (iarg + atom->ntypes >= narg) error->all(FLERR, "Illegal fix neighbor/swap command");
      iarg++;
      int i = 0;
      while (iarg < narg) {
        if (isalpha(arg[iarg][0])) break;
        if (i >= atom->ntypes) error->all(FLERR, "Illegal fix neighbor/swap command");
        rate_list[i] = utils::numeric(FLERR, arg[iarg], false, lmp);
        i++;
        iarg++;
      }
      rates_flag = 1;
      if (i != atom->ntypes) error->all(FLERR, "Illegal fix neighbor/swap command");
    } else
      error->all(FLERR, "Illegal fix neighbor/swap command");
  }
}

/* ---------------------------------------------------------------------- */

int FixNeighborSwap::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNeighborSwap::init()
{
  c_pe = modify->get_compute_by_id("thermo_pe");

  int *type = atom->type;

  if (nswaptypes < 2 && !diff_flag)
    error->all(FLERR, "Must specify at least 2 types in fix neighbor/swap command");

  // set index and check validity of region

  if (idregion) {
    region = domain->get_region_by_id(idregion);
    if (!region) error->all(FLERR, "Region {} for fix setforce does not exist", idregion);
  }

  for (int iswaptype = 0; iswaptype < nswaptypes; iswaptype++)
    if (type_list[iswaptype] <= 0 || type_list[iswaptype] > atom->ntypes)
      error->all(FLERR, "Invalid atom type in fix neighbor/swap command");

  if (atom->q_flag) {
    double qmax, qmin;
    int firstall, first;
    memory->create(qtype, nswaptypes, "neighbor/swap:qtype");
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

  memory->create(sqrt_mass_ratio, atom->ntypes + 1, atom->ntypes + 1,
                 "neighbor/swap:sqrt_mass_ratio");
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

    if (flagall) error->all(FLERR, "Cannot do neighbor/swap on atoms in atom_modify first group");
  }
}

/* ----------------------------------------------------------------------
   attempt Monte Carlo swaps
------------------------------------------------------------------------- */

void FixNeighborSwap::pre_exchange()
{
  // just return if should not be called on this timestep

  if (next_reneighbor != update->ntimestep) return;

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
  update_iswap_atoms_list();
  for (int i = 0; i < ncycles; i++) nsuccess += attempt_swap();

  // udpate MC stats

  nswap_attempts += ncycles;
  nswap_successes += nsuccess;

  next_reneighbor = update->ntimestep + nevery;
}

/* ----------------------------------------------------------------------
   attempt a swap of a pair of atoms
   compare before/after energy and accept/reject the swap
------------------------------------------------------------------------- */

int FixNeighborSwap::attempt_swap()
{
  // int nlocal = atom->nlocal;
  tagint *id = atom->tag;

  if (niswap == 0) return 0;

  // pre-swap energy

  double energy_before = energy_stored;

  // pick a random atom i

  int i = pick_i_swap_atom();

  // get global id and position of atom i
  // get_global_i(i);

  // build nearest-neighbor list based on atom i

  build_i_neighbor_list(i);
  if (njswap <= 0) return 0;

  // pick a neighbor atom j based on i neighbor list
  jtype_selected = -1;
  int j = pick_j_swap_neighbor(i);

  int itype = type_list[0];
  int jtype = jtype_selected;

  // Accept swap if types are equal, no change to system
  if (itype == jtype) { return 1; }

  // swap their properties
  if (i >= 0) {
    atom->type[i] = jtype;
    if (atom->q_flag) atom->q[i] = qtype[jtype_selected];
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

  // if swap accepted, return 1
  // if ke_flag, rescale atom velocities

  if (random_equal->uniform() < exp(beta * (energy_before - energy_after))) {
    update_iswap_atoms_list();
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
    atom->type[i] = itype;
    if (atom->q_flag) atom->q[i] = qtype[0];
  }
  if (j >= 0) {
    atom->type[j] = jtype;
    if (atom->q_flag) atom->q[j] = qtype[jtype_selected];
  }

  return 0;
}

/* ----------------------------------------------------------------------
   compute system potential energy
------------------------------------------------------------------------- */

double FixNeighborSwap::energy_full()
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

int FixNeighborSwap::pick_i_swap_atom()
{
  tagint *id = atom->tag;
  int id_center_local = -1;
  int i = -1;

  int iwhichglobal = static_cast<int>(niswap * random_equal->uniform());
  if ((iwhichglobal >= niswap_before) && (iwhichglobal < niswap_before + niswap_local)) {
    int iwhichlocal = iwhichglobal - niswap_before;
    i = local_swap_iatom_list[iwhichlocal];
    id_center_local = id[i];
    MPI_Allreduce(&id[i], &id_center, 1, MPI_INT, MPI_MAX, world);
  } else {
    MPI_Allreduce(&id[i], &id_center, 1, MPI_INT, MPI_MAX, world);
  }

  return i;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int FixNeighborSwap::pick_j_swap_neighbor(int i)
{
  int j = -1;
  int jtype_selected_local = -1;

  // Generate random double from 0 to maximum global probability
  double selected_prob = static_cast<double>(global_probability * random_equal->uniform());

  // Find which local swap atom corresponds to probability
  if ((selected_prob >= prev_probability) &&
      (selected_prob < prev_probability + local_probability)) {
    double search_prob = selected_prob - prev_probability;
    for (int n = 0; n < njswap_local; n++) {
      if (search_prob > local_swap_probability[n]) {
        search_prob -= local_swap_probability[n];
      } else {
        j = local_swap_neighbor_list[n];
        jtype_selected_local = local_swap_type_list[n];
        MPI_Allreduce(&jtype_selected_local, &jtype_selected, 1, MPI_INT, MPI_MAX, world);
        return j;
      }
    }
    error->all(FLERR, "Did not select local neighbor swap atom");
  }

  MPI_Allreduce(&jtype_selected_local, &jtype_selected, 1, MPI_INT, MPI_MAX, world);
  return j;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

double FixNeighborSwap::get_distance(double *i, double *j)
{
  double r = sqrt(MathSpecial::square((i[0] - j[0])) + MathSpecial::square((i[1] - j[1])) +
                  MathSpecial::square((i[2] - j[2])));
  return r;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixNeighborSwap::build_i_neighbor_list(int i_center)
{
  int nghost = atom->nghost;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double **x = atom->x;
  tagint *id = atom->tag;

  // Allocate local_swap_neighbor_list size

  memory->sfree(local_swap_neighbor_list);
  atom_swap_nmax = atom->nmax;
  local_swap_neighbor_list =
      (int *) memory->smalloc(atom_swap_nmax * sizeof(int), "MCSWAP:local_swap_neighbor_list");

  memory->sfree(local_swap_probability);
  local_swap_probability = (double *) memory->smalloc(atom_swap_nmax * sizeof(double),
                                                      "MCSWAP:local_swap_probability_list");

  memory->sfree(local_swap_type_list);
  local_swap_type_list =
      (int *) memory->smalloc(atom_swap_nmax * sizeof(int), "MCSWAP:local_swap_type_list");

  // Compute voronoi and access neighbor list

  c_voro->compute_local();

  voro_neighbor_list = c_voro->array_local;
  njswap_local = 0;
  local_probability = 0.0;

  for (int n = 0; n < c_voro->size_local_rows; n++) {

    int temp_j_id = -1;
    int temp_j = -1;

    // Find local voronoi entry with selected central atom
    if ((int) voro_neighbor_list[n][0] == id_center) {
      temp_j_id = voro_neighbor_list[n][1];
      temp_j = -1;
    } else if (((int) voro_neighbor_list[n][1] == id_center) && (i_center < 0)) {
      temp_j_id = voro_neighbor_list[n][0];
      temp_j = -1;
    } else {
      continue;
    }

    // Find which local atom corresponds to neighbor
    for (int j = 0; j < nlocal; j++) {
      if (temp_j_id == id[j]) {
        temp_j = j;
        break;
      }
    }

    // If temp_j not on this processor, skip
    if (temp_j < 0) continue;

    if (region) {
      if (region->match(x[temp_j][0], x[temp_j][1], x[temp_j][2]) == 1) {
        if (atom->mask[temp_j] & groupbit) {
          if (diff_flag) {
            // Calculate distance from i to each j, adjust probability of selection

            // Get distance if own center atom
            double r = INFINITY;
            if (i_center >= 0) { double r = get_distance(x[temp_j], x[i_center]); }

            // Get local id of ghost center atom when ghost
            for (int i = nlocal; i < nlocal + nghost; i++) {
              if ((id[i] == id_center) && (get_distance(x[temp_j], x[i]) < r)) {
                r = get_distance(x[temp_j], x[i]);
              }
            }

            if (rates_flag) {
              local_swap_probability[njswap_local] =
                  rate_list[type[temp_j] - 1] * exp(-MathSpecial::square(r / r_0));
            } else {
              local_swap_probability[njswap_local] = exp(-MathSpecial::square(r / r_0));
            }
            local_probability += local_swap_probability[njswap_local];
            local_swap_type_list[njswap_local] = type[temp_j];
            local_swap_neighbor_list[njswap_local] = temp_j;
            njswap_local++;
          } else {
            for (int jswaptype = 1; jswaptype < nswaptypes; jswaptype++) {
              if (type[temp_j] == type_list[jswaptype]) {
                // Calculate distance from i to each j, adjust probability of selection
                // Get distance if own center atom
                double r = INFINITY;
                if (i_center >= 0) { double r = get_distance(x[temp_j], x[i_center]); }

                // Get local id of ghost center atom when ghost
                for (int i = nlocal; i < nlocal + nghost; i++) {
                  if ((id[i] == id_center) && (get_distance(x[temp_j], x[i]) < r)) {
                    r = get_distance(x[temp_j], x[i]);
                  }
                }

                if (rates_flag) {
                  local_swap_probability[njswap_local] =
                      rate_list[type[temp_j] - 1] * exp(-MathSpecial::square(r / r_0));
                } else {
                  local_swap_probability[njswap_local] = exp(-MathSpecial::square(r / r_0));
                }
                local_probability += local_swap_probability[njswap_local];

                local_swap_type_list[njswap_local] = jswaptype;
                local_swap_neighbor_list[njswap_local] = temp_j;
                njswap_local++;
              }
            }
          }
        }
      }
    } else {
      if (atom->mask[temp_j] & groupbit) {
        if (diff_flag) {
          // Calculate distance from i to each j, adjust probability of selection
          // Get distance if own center atom
          double r = INFINITY;
          if (i_center >= 0) { r = get_distance(x[temp_j], x[i_center]); }

          // Get local id of ghost center atoms
          for (int i = nlocal; i < nlocal + nghost; i++) {
            if ((id[i] == id_center) && (get_distance(x[temp_j], x[i]) < r))
              r = get_distance(x[temp_j], x[i]);
          }

          if (rates_flag) {
            local_swap_probability[njswap_local] =
                rate_list[type[temp_j] - 1] * exp(-MathSpecial::square(r / r_0));
          } else {
            local_swap_probability[njswap_local] = exp(-MathSpecial::square(r / r_0));
          }
          local_probability += local_swap_probability[njswap_local];

          local_swap_type_list[njswap_local] = type[temp_j];
          local_swap_neighbor_list[njswap_local] = temp_j;
          njswap_local++;
        } else {
          for (int jswaptype = 1; jswaptype < nswaptypes; jswaptype++) {
            if (type[temp_j] == type_list[jswaptype]) {
              // Calculate distance from i to each j, adjust probability of selection
              // Get distance if own center atom
              double r = INFINITY;
              if (i_center >= 0) { double r = get_distance(x[temp_j], x[i_center]); }

              // Get local id of ghost center atom when ghost
              for (int i = nlocal; i < nlocal + nghost; i++) {
                if ((id[i] == id_center) && (get_distance(x[temp_j], x[i]) < r)) {
                  r = get_distance(x[temp_j], x[i]);
                }
              }

              if (rates_flag) {
                local_swap_probability[njswap_local] =
                    rate_list[type[temp_j] - 1] * exp(-MathSpecial::square(r / r_0));
              } else {
                local_swap_probability[njswap_local] = exp(-MathSpecial::square(r / r_0));
              }
              local_probability += local_swap_probability[njswap_local];

              local_swap_type_list[njswap_local] = jswaptype;
              local_swap_neighbor_list[njswap_local] = temp_j;
              njswap_local++;
            }
          }
        }
      }
    }
  }

  MPI_Allreduce(&njswap_local, &njswap, 1, MPI_INT, MPI_SUM, world);
  MPI_Scan(&njswap_local, &njswap_before, 1, MPI_INT, MPI_SUM, world);
  njswap_before -= njswap_local;

  MPI_Allreduce(&local_probability, &global_probability, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Scan(&local_probability, &prev_probability, 1, MPI_DOUBLE, MPI_SUM, world);
  prev_probability -= local_probability;
}

/* ----------------------------------------------------------------------
   update the list of swap atoms
------------------------------------------------------------------------- */

void FixNeighborSwap::update_iswap_atoms_list()
{
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double **x = atom->x;

  if (atom->nmax > atom_swap_nmax) {
    memory->sfree(local_swap_iatom_list);
    atom_swap_nmax = atom->nmax;
    local_swap_iatom_list =
        (int *) memory->smalloc(atom_swap_nmax * sizeof(int), "MCSWAP:local_swap_iatom_list");
  }

  niswap_local = 0;

  if (region) {

    for (int i = 0; i < nlocal; i++) {
      if (region->match(x[i][0], x[i][1], x[i][2]) == 1) {
        if (atom->mask[i] & groupbit) {
          if (type[i] == type_list[0]) {
            local_swap_iatom_list[niswap_local] = i;
            niswap_local++;
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
        }
      }
    }
  }

  MPI_Allreduce(&niswap_local, &niswap, 1, MPI_INT, MPI_SUM, world);
  MPI_Scan(&niswap_local, &niswap_before, 1, MPI_INT, MPI_SUM, world);
  niswap_before -= niswap_local;
}

/* ---------------------------------------------------------------------- */

int FixNeighborSwap::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                       int * /*pbc*/)
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

void FixNeighborSwap::unpack_forward_comm(int n, int first, double *buf)
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

double FixNeighborSwap::compute_vector(int n)
{
  if (n == 0) return nswap_attempts;
  if (n == 1) return nswap_successes;
  return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixNeighborSwap::memory_usage()
{
  double bytes = (double) atom_swap_nmax * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixNeighborSwap::write_restart(FILE *fp)
{
  int n = 0;
  double list[6];
  list[n++] = random_equal->state();
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

void FixNeighborSwap::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  seed = static_cast<int>(list[n++]);
  random_equal->reset(seed);

  next_reneighbor = (bigint) ubuf(list[n++]).i;

  nswap_attempts = static_cast<int>(list[n++]);
  nswap_successes = static_cast<int>(list[n++]);

  bigint ntimestep_restart = (bigint) ubuf(list[n++]).i;
  if (ntimestep_restart != update->ntimestep)
    error->all(FLERR, "Must not reset timestep when restarting fix neighbor/swap");
}
