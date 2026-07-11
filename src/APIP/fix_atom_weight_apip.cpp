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

#include "fix_atom_weight_apip.h"

#include "atom.h"
#include "atom_vec_apip.h"
#include "comm.h"
#include "error.h"
#include "fix_store_atom.h"
#include "force.h"
#include "group.h"
#include "modify.h"
#include "pair.h"
#include "timer.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

FixAtomWeightAPIP::FixAtomWeightAPIP(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), time_simple_extract_name(nullptr), time_complex_extract_name(nullptr),
    time_group_extract_name(nullptr), time_lambda_extract_name(nullptr), time_group_name(nullptr),
    fixstore(nullptr), ap_timer(nullptr), fix_lambda(nullptr)
{
  if (narg < 9) error->all(FLERR, "Illegal fix balance command");

  ap_timer = new APIPtimer(lmp);

  // set defaults
  time_simple_atom = time_complex_atom = time_group_atom = time_lambda_atom = -1;
  n_simple = n_complex = n_group = n_lambda = 0;
  nevery = -1;
  rescale_work = true;

  peratom_flag = 1;
  size_peratom_cols = 0;    // vector

  time_group_i = -1;
  time_group_bit = 0;

  vector_flag = 1;
  size_vector = 4;
  extvector = 0;

  for (int i = 0; i < size_vector; i++) avg_time_atom[i] = 0;

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);

  if (strcmp(arg[4], "eam") == 0) {
    time_simple_extract_name = utils::strdup("eam/apip:time_per_atom");
  } else if (strcmp(arg[4], "ace") == 0) {
    time_simple_extract_name = utils::strdup("pace/fast/apip:time_per_atom");
  } else {
    time_simple_atom = utils::numeric(FLERR, arg[4], false, lmp);
    avg_time_atom[0] = time_simple_atom;
  }

  if (strcmp(arg[5], "ace") == 0) {
    time_complex_extract_name = utils::strdup("pace/apip:time_per_atom");
  } else {
    time_complex_atom = utils::numeric(FLERR, arg[5], false, lmp);
    avg_time_atom[1] = time_complex_atom;
  }

  if (strcmp(arg[6], "lambda/input") == 0) {
    time_group_extract_name = utils::strdup("lambda/input/apip:time_per_atom");
  } else {
    time_group_atom = utils::numeric(FLERR, arg[6], false, lmp);
    avg_time_atom[2] = time_group_atom;
  }

  if (strcmp(arg[7], "lambda/zone") == 0) {
    time_lambda_extract_name = utils::strdup("lambda/zone/apip:time_per_atom");
  } else {
    time_lambda_atom = utils::numeric(FLERR, arg[7], false, lmp);
    avg_time_atom[3] = time_lambda_atom;
  }

  // read name of group
  time_group_name = utils::strdup(arg[8]);
  time_group_i = group->find(time_group_name);
  if (time_group_i == -1)
    error->all(FLERR, "atom_weight/apip: group {} does not exist", time_group_name);
  time_group_bit = group->bitmask[time_group_i];

  // parse remaining arguments
  for (int iarg = 9; iarg < narg; iarg++) {
    if (strcmp(arg[iarg], "no_rescale") == 0) {
      rescale_work = false;
    } else {
      error->all(FLERR, "atom_weight/apip: unknown argument {}", arg[iarg]);
    }
  }

  // check arguments
  if (nevery < 1) error->all(FLERR, "atom_weight/apip: nevery > 0 required");
  if (!atom->apip_lambda_required_flag)
    error->all(FLERR, "atom_weight/apip: atomic style with lambda_required required");

  if (time_simple_extract_name || time_complex_extract_name || time_group_extract_name ||
      time_lambda_extract_name) {
    if (force->pair == nullptr)
      error->all(FLERR, "atom_weight/apip: extract requires a defined pair style");
  }

  int useless_dim = -1;
  if (time_simple_extract_name) {
    if (force->pair->extract(time_simple_extract_name, useless_dim) == nullptr)
      error->all(FLERR, "atom_weight/apip: simple time cannot be extracted with {} from {}",
                 time_simple_extract_name, force->pair_style);
  } else {
    if (time_simple_atom < 0)
      error->all(FLERR, "atom_weight/apip: time_simple_atom needs to be non-negative instead of {}",
                 time_simple_atom);
  }

  if (time_complex_extract_name) {
    if (force->pair->extract(time_complex_extract_name, useless_dim) == nullptr)
      error->all(FLERR, "atom_weight/apip: complex time cannot be extracted with {} from {}",
                 time_complex_extract_name, force->pair_style);
  } else {
    if (time_complex_atom < 0)
      error->all(FLERR,
                 "atom_weight/apip: time_complex_atom needs to be non-negative instead of {}",
                 time_complex_atom);
  }

  if (time_group_extract_name) {
    if (force->pair->extract(time_group_extract_name, useless_dim) == nullptr)
      error->all(FLERR, "atom_weight/apip: group time cannot be extracted with {} from {}",
                 time_group_extract_name, force->pair_style);
  } else {
    if (time_group_atom < 0)
      error->all(FLERR, "atom_weight/apip: time_group_atom needs to be non-negative instead of {}",
                 time_group_atom);
  }

  if (time_lambda_extract_name) {
    if (force->pair->extract(time_lambda_extract_name, useless_dim) == nullptr)
      error->all(FLERR, "atom_weight/apip: lambda time cannot be extracted with {} from {}",
                 time_lambda_extract_name, force->pair_style);
  } else {
    if (time_lambda_atom < 0)
      error->all(FLERR, "atom_weight/apip: time_lambda_atom needs to be non-negative instead of {}",
                 time_lambda_atom);
  }

  if (comm->me == 0) {
    utils::logmesg(lmp, "atomic load lambda:\n");
    if (time_simple_extract_name)
      utils::logmesg(lmp, "\tfast potential: extract {}\n", time_simple_extract_name);
    else
      utils::logmesg(lmp, "\tfast potential: const {}\n", time_simple_atom);
    if (time_complex_extract_name)
      utils::logmesg(lmp, "\tprecise potential: extract {}\n", time_complex_extract_name);
    else
      utils::logmesg(lmp, "\tprecise potential: const {}\n", time_complex_atom);
    if (time_group_extract_name)
      utils::logmesg(lmp, "\tlambda_input: extract {}\n", time_group_extract_name);
    else
      utils::logmesg(lmp, "\tlambda_input: const {}\n", time_group_atom);
    if (time_lambda_extract_name)
      utils::logmesg(lmp, "\tlambda: extract {}\n", time_lambda_extract_name);
    else
      utils::logmesg(lmp, "\tlambda: const {}\n", time_lambda_atom);
  }

  global_freq = nevery;
  peratom_freq = nevery;
}

/**
 *  Deconstructor. Delete allocated memory.
 */

FixAtomWeightAPIP::~FixAtomWeightAPIP()
{
  delete ap_timer;
  delete[] time_simple_extract_name;
  delete[] time_complex_extract_name;
  delete[] time_lambda_extract_name;
  delete[] time_group_extract_name;
  delete[] time_group_name;

  // check nfix in case all fixes have already been deleted
  if (fixstore && modify->nfix) modify->delete_fix(fixstore->id);
  fixstore = nullptr;
}

/**
  * allocate per-particle weight storage via FixStoreAtom
  * fix could already be allocated if fix balance is re-specified
  */

void FixAtomWeightAPIP::post_constructor()
{
  std::string cmd;
  cmd = id;
  cmd += "LAMBDA_WEIGHT";
  fixstore = dynamic_cast<FixStoreAtom *>(modify->get_fix_by_id(cmd));
  if (!fixstore)
    fixstore = dynamic_cast<FixStoreAtom *>(modify->add_fix(cmd + " all STORE/ATOM 1 0 0 1"));

  // do not carry weights with atoms during normal atom migration
  fixstore->disable = 1;
  vector_atom = fixstore->vstore;
}

/**
 *  For lammps.
 *  @return mask
 */

int FixAtomWeightAPIP::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= PRE_FORCE;    // for setup_pre_force only
  mask |= END_OF_STEP;
  return mask;
}

/**
 *  Initialise calculated variables and setup timer objects.
 */

void FixAtomWeightAPIP::init()
{
  int counter = 0;
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style, "atom_weight/apip") == 0) counter++;
  }
  if (counter > 1) error->all(FLERR, "More than one atom_weight/apip fix");

  // get ptr to fix lambda
  counter = 0;
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style, "lambda/apip") == 0) {
      fix_lambda = modify->fix[i];
      counter++;
    }
  }
  if (counter > 1) error->all(FLERR, "More than one fix lambda");
  if (counter == 0 && (time_lambda_extract_name || time_lambda_atom > 0))
    error->all(FLERR, "fix lambda required to approximate weight of pair style lambda/zone");

  // This fix is evaluated in pre_exchange, but needs to be evaluated before load-balancing fixes.
  for (auto *ifix : modify->get_fix_list()) {
    if (strcmp(id, ifix->id) == 0) {
      // The remaining fixes are called after fix atom_load_lambda and ,thus, are not of interest.
      break;
    } else if (ifix->box_change == BOX_CHANGE_DOMAIN) {
      error->all(FLERR, "atom_weight/apip: fix {} should come after fix {}", ifix->id, id);
    }
  }

  // check that group for time_group has not been deleted
  if (time_group_name) {
    time_group_i = group->find(time_group_name);
    if (time_group_i == -1)
      error->all(FLERR, "atom_weight/apip: group {} does not exist", time_group_name);
    time_group_bit = group->bitmask[time_group_i];
  }

  ap_timer->init();

  last_calc = -1;
}

/**
  * Set 1 as initial weight as no forces have been calculated yet.
  */

void FixAtomWeightAPIP::setup_pre_exchange()
{
  // fix balance rebalances in setup_pre_exchange.
  // setup_pre_exchange is called prior to force-calculations.
  // Thus, there are no measured times yet.
  // Fix balance with weight time 1.0 uses 1.0 as weight for each atom.
  int i, nlocal;
  double *weight;

  nlocal = atom->nlocal;

  weight = fixstore->vstore;
  for (i = 0; i < nlocal; i++) weight[i] = 1;

  // store pointer to weight for lammps
  vector_atom = fixstore->vstore;
  // fix balance can move particles and the weight information is required
  // at the end of the step for the dump output.
  // Thus, weights need to migrate with atoms.
  fixstore->disable = 0;
}

/**
  * Initial atom migration is done.
  * The atoms do not need to migrate with atoms any more.
  */

void FixAtomWeightAPIP::setup_pre_force(int /*vflag*/)
{
  fixstore->disable = 1;
  // Atoms are with their weights now.
  // -> update vector_atom
  vector_atom = fixstore->vstore;
}

/**
  *  Compute weight of particles for load balancing.
  */

void FixAtomWeightAPIP::pre_exchange()
{
  if (update->ntimestep % peratom_freq != 0) { return; }

  calc_work_per_particle();
}

/**
  *  Update output pointer for output.
  */

void FixAtomWeightAPIP::end_of_step()
{
  if (update->ntimestep % peratom_freq != 0) { return; }

  // The work is not calculated twice.
  // Call a second time as pre_exchange is not called when there is no exchange.
  calc_work_per_particle();

  // weights should not migrate with atoms
  fixstore->disable = 1;
  vector_atom = fixstore->vstore;
}

/**
 *  Calculate the work for a simple/complex atom.
 *  Times and number of atoms are extracted from pair styles.
 *  @note updates particle number and time variables of this class.
 */

void FixAtomWeightAPIP::calc_work_per_particle()
{
  // calculating twice would destroy time measurements
  if (update->ntimestep == last_calc) { return; }
  last_calc = update->ntimestep;

  char *extract_name[4] = {time_simple_extract_name, time_complex_extract_name,
                           time_group_extract_name, time_lambda_extract_name};
  double *time_pa_ptr[4] = {&time_simple_atom, &time_complex_atom, &time_group_atom,
                            &time_lambda_atom};

  double buffer[8];
  int useless_dim = -1;

  // extract times per atom from pair styles if required

  int counter = 0;
  for (int i = 0; i < 4; i++) {
    if (extract_name[i]) {
      // get time per atom
      *(time_pa_ptr[i]) = *((double *) force->pair->extract(extract_name[i], useless_dim));
      // save time per atom to buffer
      if (*(time_pa_ptr[i]) < 0) {
        // no calculations
        buffer[counter] = buffer[counter + 1] = 0;
      } else {
        // save time per calculation
        buffer[counter] = *(time_pa_ptr[i]);
        buffer[counter + 1] = 1;
      }
      counter += 2;
    }
  }

  if (counter) {
    // do not use averaged values for all processors since they depend on the number of neighbours
    // (which can vary, e.g. at a surface
    // -> only use global value if there is no local one

    MPI_Allreduce(MPI_IN_PLACE, buffer, counter, MPI_DOUBLE, MPI_SUM, world);

    for (int i = 3; i >= 0; i--) {
      if (extract_name[i]) {
        // calculate average over all processors
        avg_time_atom[i] = buffer[counter - 1] > 0 ? buffer[counter - 2] / buffer[counter - 1] : 0;
        // use average if there is no local value
        if (*(time_pa_ptr[i]) < 0) *(time_pa_ptr[i]) = avg_time_atom[i];
        counter -= 2;
      }
    }
  }

  //set weight for each particle

  double work_atom, *weight, **lambda_input_history;
  int dim, *mask, *lambda_required;

  weight = fixstore->vstore;
  mask = atom->mask;
  lambda_required = atom->apip_lambda_required;

  int nlocal = atom->nlocal;
  // assume a homogeneous time per simple and complex particle
  n_simple = n_complex = 0;
  for (int i = 0; i < nlocal; i++) {
    work_atom = 0;
    if (lambda_required[i] & ApipLambdaRequired::SIMPLE) {
      work_atom += time_simple_atom;
      n_simple++;
    }
    if (lambda_required[i] & ApipLambdaRequired::COMPLEX) {
      work_atom += time_complex_atom;
      n_complex++;
    }
    weight[i] = work_atom;
  }

  n_group = 0;
  if (time_group_atom > 0) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & time_group_bit) {
        weight[i] += time_group_atom;
        n_group++;
      }
    }
  }

  n_lambda = 0;
  if (time_lambda_atom > 0) {
    // get lambda(time averaged lambda input)
    lambda_input_history = (double **) fix_lambda->extract("fix_lambda:lambda_input_history", dim);
    const int histlen_lambda_input =
        *((int *) fix_lambda->extract("fix_lambda:lambda_input_history_len", dim));
    if (lambda_input_history == nullptr || histlen_lambda_input < 1)
      error->all(FLERR, "atom_weight/apip: extracting fix_lambda:lambda_input_history failed");

    for (int i = 0; i < nlocal; i++) {
      if (lambda_input_history[i][histlen_lambda_input + 1] != 1) {
        weight[i] += time_lambda_atom;
        n_lambda++;
      }
    }
  }

  if (rescale_work) {
    // calculate sum of vector work
    double work_atoms = 0;
    for (int i = 0; i < nlocal; i++) work_atoms += weight[i];
    // calculate rescale factor
    double rescale_factor = ap_timer->get_work() / work_atoms;
    // apply rescale factor
    for (int i = 0; i < nlocal; i++) weight[i] *= rescale_factor;
  }

  // store pointer to weight for lammps
  vector_atom = fixstore->vstore;
  // weights need to migrate with atoms
  fixstore->disable = 0;
}

/**
  * Provide average compute times for the output.
  * 1st: time per simple atom
  * 2nd: time per complex atom
  * 3rd: time per lambda_input calculation
  * 4th: time per lambda calculation
  */

double FixAtomWeightAPIP::compute_vector(int i)
{
  if (i < size_vector) return avg_time_atom[i];
  return 0;
}

/**
 *  Set everything to zero.
 *  @param[in] lmp lammps to get the Pointers class
 */

APIPtimer::APIPtimer(LAMMPS *lmp) : Pointers(lmp)
{
  for (int i = 0; i < 4; i++) {
    time[i] = 0;
    time_interval[i] = 0;
  }
}

/**
 *  Reset times to zero.
 */

void APIPtimer::init()
{
  // do not use set_time();
  // lammps resets the timers to zero, but the timers are not zero at timer initialisation time
  for (int i = 0; i < 4; i++) {
    time[i] = 0;
    time_interval[i] = 0;
  }
}

/**
 *  Save values of lammps timers.
 *  @note sets time
 */

void APIPtimer::set_time()
{
  time[0] = timer->get_wall(Timer::PAIR);
  time[1] = timer->get_wall(Timer::NEIGH);
  time[2] = timer->get_wall(Timer::BOND);
  time[3] = timer->get_wall(Timer::KSPACE);
}

/**
 *  Get time work since last call.
 *  @note sets time and time_interval
 *  @return work
 */

double APIPtimer::get_work()
{
  double work = 0;
  // store old times
  for (int i = 0; i < 4; i++) time_interval[i] = -time[i];
  set_time();
  // calculate differences to new times
  for (int i = 0; i < 4; i++) {
    time_interval[i] += time[i];
    work += time_interval[i];
  }
  // add constant time to prevent balancing with numerically zero
  // e.g. (for few atoms per processor and balancing every step)
  return work + 0.1;
}
