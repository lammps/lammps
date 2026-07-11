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
   Contributing author: David Immel (d.immel@fz-juelich.de, FZJ, Germany)
------------------------------------------------------------------------- */
#include "fix_lambda_thermostat_apip.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "my_page.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <algorithm>
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

#define PGDELTA 1

/* ---------------------------------------------------------------------- */

FixLambdaThermostatAPIP::FixLambdaThermostatAPIP(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), list(nullptr), energy_change_atom(nullptr), peratom_stats(nullptr),
    local_numneigh(nullptr), local_firstneigh(nullptr), ipage(nullptr), jlist_copy(nullptr)
{
  // set global options for fix class

  vector_flag = 1;
  size_vector = 6;
  extvector = 0;

  // set default vaues
  dtf = 0;
  update_stats = false;
  int seed = 42;
  rescaling_N_neighbours = 200;
  // output
  peratom_flag = 0;
  peratom_freq = -1;    // default from fix.cpp
  size_peratom_cols = 4;

  // parse arguments
  for (int iarg = 3; iarg < narg; iarg++) {

    if (strcmp(arg[iarg], "seed") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "fix lambda_thermostat/apip: seed requires one argument");
      seed = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (seed <= 0) { error->all(FLERR, "fix lambda_thermostat/apip seed <= 0"); }
      iarg++;

    } else if (strcmp(arg[iarg], "store_atomic_forces") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "fix lambda_thermostat/apip: store_atomic_forces requires one argument");
      peratom_flag = 1;
      peratom_freq = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (peratom_freq < 1)
        error->all(FLERR, "fix lambda_thermostat/apip: frequency of store_atomic_forces < 1");
      iarg++;

    } else if (strcmp(arg[iarg], "N_rescaling") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "fix lambda_thermostat/apip: mode number requires one argument");
      rescaling_N_neighbours = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 1;

    } else
      error->all(FLERR, "fix lambda_thermostat/apip: unknown argument {}", arg[iarg]);
  }

  // error checks

  if (!atom->apip_e_fast_flag) {
    error->all(FLERR, "fix lambda_thermostat/apip requires atomic style with e_simple.");
  }
  if (!atom->apip_e_precise_flag) {
    error->all(FLERR, "fix lambda_thermostat/apip requires atomic style with e_complex.");
  }
  if (!atom->apip_lambda_const_flag) {
    error->all(FLERR, "fix lambda_thermostat/apip requires atomic style with lambda_const.");
  }
  if (!atom->apip_lambda_flag) {
    error->all(FLERR, "fix lambda_thermostat/apip requires atomic style with lambda.");
  }
  if (!atom->apip_f_const_lambda_flag) {
    error->all(FLERR, "fix lambda_thermostat/apip requires atomic style with f_const_lambda.");
  }
  if (!atom->apip_f_dyn_lambda_flag) {
    error->all(FLERR, "fix lambda_thermostat/apip requires atomic style with f_dyn_lambda.");
  }
  if (rescaling_N_neighbours <= 1)
    error->all(FLERR, "fix lambda_thermostat/apip: rescaling_N_neighbours <= 1");

  // rng for shuffle
  random_mt = std::mt19937(seed);

  // init output values
  energy_change_kin = energy_change_pot = 0;
  nmax_energy = 0;

  reduceflag = 0;
  for (int i = 0; i < size_vector; i++) { outvec[i] = 0; }
  sum_energy_change = sum_energy_violation = 0;
  n_energy_violation = n_energy_differences = 0;

  if (peratom_flag) {
    // zero the array since dump may access it on timestep 0
    // zero the array since a variable may access it before first run

    nmax_stats = atom->nmax;
    memory->create(peratom_stats, nmax_stats, size_peratom_cols,
                   "lambda_thermostat/apip:peratom_stats");
    array_atom = peratom_stats;

    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      for (int j = 0; j < size_peratom_cols; j++) peratom_stats[i][j] = 0;
  } else {
    nmax_stats = 0;
  }

  nmax_list = 0;
  pgsize = oneatom = 0;
}

/* ---------------------------------------------------------------------- */

FixLambdaThermostatAPIP::~FixLambdaThermostatAPIP()
{
  memory->destroy(energy_change_atom);
  memory->destroy(peratom_stats);

  memory->destroy(local_numneigh);
  memory->sfree(local_firstneigh);
  memory->destroy(jlist_copy);
  delete[] ipage;
}

/* ---------------------------------------------------------------------- */

int FixLambdaThermostatAPIP::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLambdaThermostatAPIP::init()
{
  dtf = 0.5 * update->dt * force->ftm2v;

  // full neighbour list for thermostating
  neighbor->add_request(this, NeighConst::REQ_FULL);

  int counter = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style, "lambda_thermostat/apip") == 0) counter++;
  if (counter > 1)
    error->all(FLERR, "fix lambda_thermostat/apip: more than one fix lambda_thermostat/apip");

  // local neighbor list
  // create pages if first time or if neighbor pgsize/oneatom has changed

  int create = 0;
  if (ipage == nullptr) create = 1;
  if (pgsize != neighbor->pgsize) create = 1;
  if (oneatom != neighbor->oneatom) create = 1;

  if (oneatom != neighbor->oneatom || ipage == nullptr) {
    // allocate memory for copy of one ngh list
    memory->destroy(jlist_copy);
    memory->create(jlist_copy, neighbor->oneatom, "lambda_thermostat/apip:jlist_copy");
  }

  if (create) {
    delete[] ipage;
    pgsize = neighbor->pgsize;
    oneatom = neighbor->oneatom;

    int nmypage = comm->nthreads;
    ipage = new MyPage<int>[nmypage];
    for (int i = 0; i < nmypage; i++) ipage[i].init(oneatom, pgsize, PGDELTA);
  }
}

/* ---------------------------------------------------------------------- */

void FixLambdaThermostatAPIP::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
   Create neighbor list from main neighbor list with local atoms only.
   The velocity updates are dependent on the velocity.
   Thus, one would need to communicate a velocity change of a ghost atom directly to all neighbouring processors.
   This communication would kill the performance.
   Thus, update only local particles.
------------------------------------------------------------------------- */

void FixLambdaThermostatAPIP::local_neighbour_list()
{
  int i, j, ii, jj, n, inum, jnum, nlocal;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int *neighptr;
  int *mask = atom->mask;

  if (atom->nmax > nmax_list) {
    nmax_list = atom->nmax;
    memory->destroy(local_numneigh);
    memory->sfree(local_firstneigh);
    memory->create(local_numneigh, nmax_list, "lambda_thermostat/apip:numneigh");
    local_firstneigh =
        (int **) memory->smalloc(nmax_list * sizeof(int *), "lambda_thermostat/apip:firstneigh");
  }

  nlocal = atom->nlocal;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // store all local neighbours of local atoms
  // scan full neighbor list of I

  ipage->reset();

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    n = 0;
    neighptr = ipage->vget();

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      if (j < nlocal && mask[j] & groupbit) neighptr[n++] = j;
    }

    local_firstneigh[i] = neighptr;
    local_numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status()) error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
  }
}

/* ---------------------------------------------------------------------- */

double FixLambdaThermostatAPIP::calculate_kinetic_energy(int i)
{
  double *v = atom->v[i];
  double m = (atom->rmass ? atom->rmass[i] : atom->mass[atom->type[i]]);
  return 0.5 * force->mvv2e * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

/* ---------------------------------------------------------------------- */

void FixLambdaThermostatAPIP::post_force(int /*vflag*/)
{
  init_peratom_stats();
  calculate_energy_change();
}

/* ---------------------------------------------------------------------- */

void FixLambdaThermostatAPIP::end_of_step()
{
  apply_thermostat();
}

/**
  * Init per-atom array with zeros.
  */

void FixLambdaThermostatAPIP::init_peratom_stats()
{
  if ((!peratom_flag) || (update->ntimestep % peratom_freq != 0)) {
    update_stats = false;
    return;
  }

  update_stats = true;
  int nlocal;

  nlocal = atom->nlocal;

  // grow stats array if required
  if (atom->nmax > nmax_stats) {
    memory->destroy(peratom_stats);
    nmax_stats = atom->nmax;
    memory->create(peratom_stats, nmax_stats, size_peratom_cols,
                   "lambda_thermostat/apip:peratom_stats");
    array_atom = peratom_stats;
  }

  for (int i = 0; i < nlocal; i++)
    peratom_stats[i][0] = peratom_stats[i][1] = peratom_stats[i][2] = peratom_stats[i][3] = 0;
}

/**
  * Calculate the energy difference of atoms that needs to be corrected
  * by the local thermostat.
  */

void FixLambdaThermostatAPIP::calculate_energy_change()
{
  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double *e_simple = atom->apip_e_fast;
  double *e_complex = atom->apip_e_precise;
  double *lambda_const = atom->apip_lambda_const;
  double *lambda = atom->apip_lambda;
  double **f_const_lambda = atom->apip_f_const_lambda;
  double **f_dyn_lambda = atom->apip_f_dyn_lambda;

  double masstmp, dtfm, changetmp;

  // allocate memory for energy change
  if (atom->nmax > nmax_energy) {
    memory->destroy(energy_change_atom);
    nmax_energy = atom->nmax;
    memory->create(energy_change_atom, nmax_energy, "lambda_thermostat/apip:energy_change_atom");
  }

  // reset calculated changes
  energy_change_pot = 0;
  energy_change_kin = 0;
  for (int i = 0; i < nlocal; i++) { energy_change_atom[i] = 0; }

  // calculate potential energy differences
  for (int i = 0; i < nlocal; i++) {
    if ((!(mask[i] & groupbit)) || lambda_const[i] == lambda[i]) { continue; }

    if (e_simple[i] == 0 || e_complex[i] == 0)
      error->one(FLERR, "lambda = {} != {} = lambda_const and e_simple = {} and e_complex = {}",
                 lambda[i], lambda_const[i], e_simple[i], e_complex[i]);
    changetmp =
        (lambda_const[i] - lambda[i]) * e_simple[i] + (lambda[i] - lambda_const[i]) * e_complex[i];

    energy_change_atom[i] += changetmp;
    energy_change_pot += changetmp;
  }

  // calculate kinetic energy difference
  // consider all local atoms
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) { continue; }

    masstmp = (rmass ? rmass[i] : mass[type[i]]);
    dtfm = dtf / masstmp;
    // dtfm = Delta t / 2 (in corresponding units)
    changetmp = (v[i][0] * 2 * dtfm * (f_const_lambda[i][0] - f_dyn_lambda[i][0]) +
                 dtfm * dtfm *
                     (f_const_lambda[i][0] * f_const_lambda[i][0] -
                      f_dyn_lambda[i][0] * f_dyn_lambda[i][0]) +
                 v[i][1] * 2 * dtfm * (f_const_lambda[i][1] - f_dyn_lambda[i][1]) +
                 dtfm * dtfm *
                     (f_const_lambda[i][1] * f_const_lambda[i][1] -
                      f_dyn_lambda[i][1] * f_dyn_lambda[i][1]) +
                 v[i][2] * 2 * dtfm * (f_const_lambda[i][2] - f_dyn_lambda[i][2]) +
                 dtfm * dtfm *
                     (f_const_lambda[i][2] * f_const_lambda[i][2] -
                      f_dyn_lambda[i][2] * f_dyn_lambda[i][2])) *
        masstmp * force->mvv2e * 0.5;

    energy_change_atom[i] += changetmp;
    energy_change_kin += changetmp;
  }

  if (update_stats) {
    for (int i = 0; i < nlocal; i++) peratom_stats[i][4] = energy_change_atom[i];
  }
}

/* ---------------------------------------------------------------------- */

void FixLambdaThermostatAPIP::apply_thermostat()
{
  double xtmp, ytmp, ztmp, masstmp, massj;
  double m_cm, v_cm[3], beta_rescaling, k_rel, v_rel[3], radicand;
  int i, ii, inum, *ilist;
  int j, jj, jnum, *jlist;

  double **x = atom->x;
  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;

  double *lambda_const = atom->apip_lambda_const;
  double *lambda = atom->apip_lambda;

  // calculate local neighbour list without ghost atoms
  local_neighbour_list();

  inum = list->inum;
  ilist = list->ilist;

  // reset stats
  // outvec has to be calculated again
  reduceflag = 1;
  sum_energy_change = 0;
  n_energy_differences = 0;

  // perform bath collisions
  // consider all local atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    if (!(mask[i] & groupbit)) { continue; }
    if (energy_change_atom[i] == 0) { continue; }

    // update stats
    n_energy_differences++;

    // store constant information of target atom
    masstmp = (rmass ? rmass[i] : mass[type[i]]);
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // get neighbour list
    jlist = local_firstneigh[i];
    jnum = local_numneigh[i];

    if (jnum == 0)
      error->one(FLERR,
                 "fix lambda_thermostat/apip: thermostating required for particle with no local "
                 "particles in neighbour list particle: {} {} {} groupbit {}\n",
                 xtmp, ytmp, ztmp, mask[i] & groupbit);

    double e_remain = energy_change_atom[i];

    // copy ngh list
    for (jj = 0; jj < jnum; jj++) jlist_copy[jj] = jlist[jj];
    // shuffle neighbour list for random rescaling set
    std::shuffle(jlist_copy, jlist_copy + jnum, random_mt);

    // rescale velocities relative to centre of mass velocity
    const int n_ngh = MIN(rescaling_N_neighbours, jnum);
    if (n_ngh < 2)
      error->one(FLERR,
                 "fix lambda_thermostat/apip: rescaling not possible for local ngh list size {}",
                 jnum);
    // 1. calculate centre of mass velocity

    // start with own particle ...
    m_cm = masstmp;
    v_cm[0] = masstmp * v[i][0];
    v_cm[1] = masstmp * v[i][1];
    v_cm[2] = masstmp * v[i][2];
    // ... and include neighbours
    for (jj = 0; jj < n_ngh; jj++) {
      j = jlist_copy[jj];
      j &= NEIGHMASK;

      massj = (rmass ? rmass[j] : mass[type[j]]);

      m_cm += massj;
      v_cm[0] += massj * v[j][0];
      v_cm[1] += massj * v[j][1];
      v_cm[2] += massj * v[j][2];
    }
    // normalisation
    v_cm[0] /= m_cm;
    v_cm[1] /= m_cm;
    v_cm[2] /= m_cm;

    // 2. calculate beta_rescaling
    // calculate kinetic energy of relative velocity for own particle ...
    v_rel[0] = v[i][0] - v_cm[0];
    v_rel[1] = v[i][1] - v_cm[1];
    v_rel[2] = v[i][2] - v_cm[2];
    k_rel = masstmp * (v_rel[0] * v_rel[0] + v_rel[1] * v_rel[1] + v_rel[2] * v_rel[2]);
    // ... and include neighbours
    for (jj = 0; jj < n_ngh; jj++) {
      j = jlist_copy[jj];
      j &= NEIGHMASK;

      massj = (rmass ? rmass[j] : mass[type[j]]);
      v_rel[0] = v[j][0] - v_cm[0];
      v_rel[1] = v[j][1] - v_cm[1];
      v_rel[2] = v[j][2] - v_cm[2];
      k_rel += massj * (v_rel[0] * v_rel[0] + v_rel[1] * v_rel[1] + v_rel[2] * v_rel[2]);
    }
    // normalisation
    k_rel *= force->mvv2e / 2.0;
    radicand = e_remain / k_rel + 1;
    if (radicand < 0) {
      // cooling is not possible
      // e_remain is the requested energy change
      // radicand = 0 <=> e_remain = -k_rel
      // -> save the energy error
      sum_energy_violation += (-e_remain - k_rel);    // > 0
      n_energy_violation++;
      e_remain += k_rel;    // save corrected part
      // use smallest possible radicand
      radicand = 0;
    }
    beta_rescaling = sqrt(radicand) - 1;

    // 3. apply velocity changes
    // start with own particle ...
    v_rel[0] = beta_rescaling * (v[i][0] - v_cm[0]);
    v_rel[1] = beta_rescaling * (v[i][1] - v_cm[1]);
    v_rel[2] = beta_rescaling * (v[i][2] - v_cm[2]);

    v[i][0] += v_rel[0];
    v[i][1] += v_rel[1];
    v[i][2] += v_rel[2];

    if (update_stats) {
      // forces
      peratom_stats[i][0] += v_rel[0] * masstmp / dtf;
      peratom_stats[i][1] += v_rel[1] * masstmp / dtf;
      peratom_stats[i][2] += v_rel[2] * masstmp / dtf;
    }

    // ... continue with neighbours
    for (jj = 0; jj < n_ngh; jj++) {
      j = jlist_copy[jj];
      j &= NEIGHMASK;

      v_rel[0] = beta_rescaling * (v[j][0] - v_cm[0]);
      v_rel[1] = beta_rescaling * (v[j][1] - v_cm[1]);
      v_rel[2] = beta_rescaling * (v[j][2] - v_cm[2]);
      v[j][0] += v_rel[0];
      v[j][1] += v_rel[1];
      v[j][2] += v_rel[2];

      if (update_stats) {
        // forces
        massj = (rmass ? rmass[j] : mass[type[j]]);
        peratom_stats[j][0] += v_rel[0] * massj / dtf;
        peratom_stats[j][1] += v_rel[1] * massj / dtf;
        peratom_stats[j][2] += v_rel[2] * massj / dtf;
      }
    }
    sum_energy_change += fabs(e_remain);
  }

  // lambda_const is used -> reset to lambda
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) lambda_const[i] = lambda[i];
  }
}

/* ---------------------------------------------------------------------- */

double FixLambdaThermostatAPIP::compute_vector(int i)
{
  // 0 # atoms with energy differences compared to lambda_const
  // 1 # bath collisions
  // 2 total change of potential energy compared to lambda_const (sum over all atoms)
  // 3 total change of kinetic   energy compared to lambda_const (sum over all atoms)
  // 4 total energy change due to thermostat
  if (reduceflag) {
    // perform reduction only once per step (if at all)
    outvec[0] = n_energy_differences;
    outvec[1] = energy_change_pot;
    outvec[2] = energy_change_kin;
    outvec[3] = sum_energy_change;
    outvec[4] = sum_energy_violation;
    outvec[5] = n_energy_violation;
    MPI_Allreduce(MPI_IN_PLACE, &outvec, size_vector, MPI_DOUBLE, MPI_SUM, world);
    reduceflag = 0;
  }

  if (i < size_vector) return outvec[i];

  return 0;
}

/* ---------------------------------------------------------------------- */

void FixLambdaThermostatAPIP::reset_dt()
{
  dtf = 0.5 * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixLambdaThermostatAPIP::memory_usage()
{
  double bytes = 0;
  bytes += (double) nmax_energy * sizeof(double);

  bytes += (double) nmax_stats * size_peratom_cols * sizeof(double);

  bytes += (double) nmax_list * sizeof(int);
  bytes += (double) nmax_list * sizeof(int *);
  for (int i = 0; i < comm->nthreads; i++) bytes += ipage[i].size();

  return bytes;
}
