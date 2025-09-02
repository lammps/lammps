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

#include "fix_lambda_apip.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "fix_store_atom.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "pair_lambda_input_apip.h"
#include "pair_lambda_zone_apip.h"
#include "update.h"

#include <algorithm>
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

static const char cite_fix_lambda_c[] =
    "fix lambda command: doi.org/10.1063/5.0245877\n\n"
    "@Article{Immel25,\n"
    " author = {Immel, David and Drautz, Ralf and Sutmann, Godehard},\n"
    " title = {Adaptive-precision potentials for large-scale atomistic simulations},\n"
    " journal = {The Journal of Chemical Physics},\n"
    " volume = {162},\n"
    " number = {11},\n"
    " pages = {114119},\n"
    " year = {2025}\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

FixLambdaAPIP::FixLambdaAPIP(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), pair_lambda_input(nullptr), pair_lambda_zone(nullptr), fixstore(nullptr),
    fixstore2(nullptr), group_name_simple(nullptr), group_name_complex(nullptr),
    group_name_ignore_lambda_input(nullptr), peratom_stats(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_lambda_c);

  // set defaults
  threshold_lo = threshold_hi = threshold_width = -1;
  cut_lo = 4.0;
  cut_hi = 12.0;
  lambda_non_group = 1;    // simple
  history_last = history2_last = -1;
  history_length = history2_length = 100;
  history_used = history2_used = 0;
  min_delta_lambda = 0;

  group_bit_simple = group_bit_complex = group_bit_ignore_lambda_input = 0;

  // output
  peratom_flag = 0;
  peratom_freq = -1;    // default from fix.cpp
  dump_history_flag = false;
  size_peratom_cols = 5;
  invoked_history_update = invoked_history2_update = -1;

  if (narg < 4) error->all(FLERR, "fix lambda requires two arguments");
  threshold_lo = utils::numeric(FLERR, arg[3], false, lmp);
  threshold_hi = utils::numeric(FLERR, arg[4], false, lmp);
  threshold_width = threshold_hi - threshold_lo;

  // parse remaining arguments
  for (int iarg = 5; iarg < narg; iarg++) {
    if (strcmp(arg[iarg], "time_averaged_zone") == 0) {
      if (iarg + 4 >= narg)
        error->all(FLERR, "fix lambda: time_averaged_zone requires four arguments");
      cut_lo = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      cut_hi = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      history_length = utils::inumeric(FLERR, arg[iarg + 3], false, lmp);
      history2_length = utils::inumeric(FLERR, arg[iarg + 4], false, lmp);
      iarg += 4;
    } else if (strcmp(arg[iarg], "min_delta_lambda") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "fix lambda: min_delta_lambda requires one argument");
      min_delta_lambda = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 1;
    } else if (strcmp(arg[iarg], "lambda_non_group") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "fix lambda: lambda_non_group requires an argument");
      if (strcmp(arg[iarg + 1], "precise") == 0) {
        lambda_non_group = 0;
      } else if (strcmp(arg[iarg + 1], "fast") == 0) {
        lambda_non_group = 1;
      } else {
        lambda_non_group = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
        if (lambda_non_group < 0 || lambda_non_group > 1)
          error->all(FLERR, "fix lambda: Illegal value of lambda_non_group");
      }
      iarg++;
    } else if (strcmp(arg[iarg], "store_atomic_stats") == 0) {
      peratom_flag = 1;
    } else if (strcmp(arg[iarg], "dump_atomic_history") == 0) {
      peratom_flag = 1;
      dump_history_flag = true;
    } else if (strcmp(arg[iarg], "group_fast") == 0) {
      // read name of group
      group_name_simple = utils::strdup(arg[iarg + 1]);
      int tmp = group->find(group_name_simple);
      if (tmp == -1) error->all(FLERR, "fix lambda: group {} does not exist", group_name_simple);
      group_bit_simple = group->bitmask[tmp];
      iarg++;
    } else if (strcmp(arg[iarg], "group_precise") == 0) {
      // read name of group
      group_name_complex = utils::strdup(arg[iarg + 1]);
      int tmp = group->find(group_name_complex);
      if (tmp == -1) error->all(FLERR, "fix lambda: group {} does not exist", group_name_complex);
      group_bit_complex = group->bitmask[tmp];
      iarg++;
    } else if (strcmp(arg[iarg], "group_ignore_lambda_input") == 0) {
      // read name of group
      group_name_ignore_lambda_input = utils::strdup(arg[iarg + 1]);
      int tmp = group->find(group_name_ignore_lambda_input);
      if (tmp == -1)
        error->all(FLERR, "fix lambda: group {} does not exist", group_name_ignore_lambda_input);
      group_bit_ignore_lambda_input = group->bitmask[tmp];
      iarg++;
    } else
      error->all(FLERR, "fix lambda: unknown argument {}", arg[iarg]);
  }
  cut_hi_sq = cut_hi * cut_hi;
  cut_width = cut_hi - cut_lo;

  // verify arguments
  if (threshold_lo > threshold_hi || threshold_lo < 0)
    error->all(FLERR, "fix lambda: Illegal or missing threshold values");
  if (min_delta_lambda < 0)
    error->all(FLERR, "fix lambda: min_delta_lambda >= 0 required instead of {}", min_delta_lambda);

  if (cut_lo < 0 || cut_hi < cut_lo) error->all(FLERR, "fix lambda: Illegal cutoff values");
  if (history_length < 2 || history2_length < 2)
    error->all(FLERR, "fix lambda: history_length > 1 required");

  if (comm->me == 0 && ((group_bit_simple != 0) || (group_bit_complex != 0)) &&
      group_bit_ignore_lambda_input == 0)
    error->warning(FLERR,
                   "group_ignore_lambda_input should be used to prevent the calculation of "
                   "lambda_input for atoms that are in the groups group_fast and group_precise.");

  if (!atom->apip_lambda_const_flag) {
    error->all(FLERR, "fix lambda requires atomic style with lambda_const.");
  }
  if (!atom->apip_lambda_flag) {
    error->all(FLERR, "fix lambda requires atomic style with lambda.");
  }
  if (!atom->apip_lambda_input_flag) {
    error->all(FLERR, "fix lambda requires atomic style with lambda_input.");
  }

  comm_forward = 2;    // up to two doubles per atom
  comm_forward_flag = FORWARD_TA;

  restart_global = 1;

  if (peratom_flag) {

    if (dump_history_flag) size_peratom_cols += history_length + history2_length + 2;

    peratom_freq = 1;
    // zero the array since dump may access it on timestep 0
    // zero the array since a variable may access it before first run

    nmax_stats = atom->nmax;
    memory->create(peratom_stats, nmax_stats, size_peratom_cols, "lambda:peratom_stats");
    array_atom = peratom_stats;

    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      for (int j = 0; j < size_peratom_cols; j++) peratom_stats[i][j] = 0;
  }
}

/* ----------------------------------------------------------------------
   modify cutoff setings
------------------------------------------------------------------------- */

int FixLambdaAPIP::modify_param(int narg, char **arg)
{
  if (narg < 2) utils::missing_cmd_args(FLERR, "fix_modify lambda/apip", error);

  cut_lo = utils::numeric(FLERR, arg[0], false, lmp);
  cut_hi = utils::numeric(FLERR, arg[1], false, lmp);
  cut_hi_sq = cut_hi * cut_hi;
  cut_width = cut_hi - cut_lo;

  if (cut_lo < 0 || cut_hi < cut_lo) error->all(FLERR, "fix lambda/apip: Illegal cutoff values");

  if (force->pair->cutforce < cut_hi)
    error->all(FLERR, "fix lambda: cutoff of potential smaller than cutoff of switching region");

  return 2;
}

/* ---------------------------------------------------------------------- */

FixLambdaAPIP::~FixLambdaAPIP()
{
  // check nfix in case all fixes have already been deleted
  if (fixstore && modify->nfix) modify->delete_fix(fixstore->id);
  if (fixstore2 && modify->nfix) modify->delete_fix(fixstore2->id);
  fixstore = fixstore2 = nullptr;

  memory->destroy(peratom_stats);

  delete[] group_name_simple;
  delete[] group_name_complex;
  delete[] group_name_ignore_lambda_input;
}

/* ---------------------------------------------------------------------- */

int FixLambdaAPIP::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= PRE_FORCE;    // for setup_pre_force only
  if (peratom_flag) mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLambdaAPIP::init()
{
  if (force->pair == nullptr) error->all(FLERR, "Fix lambda requires a pair style be defined");

  // only one fix lambda
  int count = 0;
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style, "lambda/apip") == 0) count++;
  }
  if (count > 1) error->all(FLERR, "More than one fix lambda.");

  // warn if there is no fix lambda_thermostat/apip
  if (comm->me == 0 && modify->get_fix_by_style("lambda_thermostat/apip").size() == 0)
    error->warning(FLERR,
                   "The energy is not conserved when lambda changes as fix lambda_thermostat/apip "
                   "is not used.");

  Pair *pair_tmp;
  // lambda_input
  pair_tmp = force->pair_match("lambda/input/", 0);
  if (!pair_tmp) error->all(FLERR, "fix lambda requires a `pair lambda_input`");
  pair_lambda_input = (PairLambdaInputAPIP *) pair_tmp;
  // lambda/zone
  pair_tmp = force->pair_match("lambda/zone/apip", 1);
  if (!pair_tmp) error->all(FLERR, "fix lambda requires a `pair lambda`");
  pair_lambda_zone = (PairLambdaZoneAPIP *) pair_tmp;

  if (force->pair->cutforce < cut_hi)
    error->all(FLERR, "fix lambda: cutoff of potential smaller than cutoff of switching region");

  if (strcmp(atom->atom_style, "apip")) error->all(FLERR, "fix lambda requires atom style apip");

  // check that groups have not been deleted
  if (group_name_simple) {
    int tmp = group->find(group_name_simple);
    if (tmp == -1) error->all(FLERR, "fix lambda: group {} does not exist", group_name_simple);
    group_bit_simple = group->bitmask[tmp];
  }
  if (group_name_complex) {
    int tmp = group->find(group_name_complex);
    if (tmp == -1) error->all(FLERR, "fix lambda: group {} does not exist", group_name_complex);
    group_bit_complex = group->bitmask[tmp];
  }
  if (group_name_ignore_lambda_input) {
    int tmp = group->find(group_name_ignore_lambda_input);
    if (tmp == -1)
      error->all(FLERR, "fix lambda: group {} does not exist", group_name_ignore_lambda_input);
    group_bit_ignore_lambda_input = group->bitmask[tmp];
  }
}

/**
  * allocate per-particle storage for past cs values via FixStoreAtom
  * fix could already be allocated if fix lambda is re-specified
  */

void FixLambdaAPIP::post_constructor()
{
  std::string cmd, cmd2;
  cmd = id;
  cmd2 = id;
  cmd += "LAMBDA_INPUT_HISTORY";
  cmd2 += "LAMBDA_HISTORY";

  // delete existing fix store if existing
  fixstore = dynamic_cast<FixStoreAtom *>(modify->get_fix_by_id(cmd));
  fixstore2 = dynamic_cast<FixStoreAtom *>(modify->get_fix_by_id(cmd2));
  // check nfix in case all fixes have already been deleted
  if (fixstore && modify->nfix) modify->delete_fix(fixstore->id);
  if (fixstore2 && modify->nfix) modify->delete_fix(fixstore2->id);
  fixstore = nullptr;

  // create new FixStoreAtom
  // store history_length of last values and the sum over all values
  char history_length_str[40], history2_length_str[40];
  sprintf(history_length_str, "%d", history_length + 2);      // lambda_input
  sprintf(history2_length_str, "%d", history2_length + 1);    // lambda

  // arguments of peratom:
  // first: 1 -> store in restart file
  // second: number of doubles to store per atom
  cmd += " all STORE/ATOM ";
  cmd2 += " all STORE/ATOM ";
  cmd += history_length_str;      // n1
  cmd2 += history2_length_str;    // n1
  cmd += " 0 0 1";                // n2 gflag rflag
  cmd2 += " 0 0 1";               // n2 gflag rflag
  fixstore = dynamic_cast<FixStoreAtom *>(modify->add_fix(cmd));
  fixstore2 = dynamic_cast<FixStoreAtom *>(modify->add_fix(cmd2));

  // carry weights with atoms during normal atom migration
  fixstore->disable = 0;
  fixstore2->disable = 0;
}

/**
  * Calculate lambda for initial atoms if required.
  * If required, this includes communication.
  */

void FixLambdaAPIP::setup_pre_force(int /*vflag*/)
{
  // lambda, lambda_input, lambda_input_ta and lambda_const are written to restart files.

  // Calculate lambda_input in pair style.
  pair_lambda_input->calculate_lambda_input();
  // Update lambda_input_history with lambda_input.
  update_lambda_input_history();
  // calculate and communicate lambda_input_ta to neighbours
  communicate_lambda_input_ta();
  // calculate lambda max with lambda_input_ta of own and ngh atoms
  pair_lambda_zone->calculate_lambda();
  // update lambda_history with calculated lambda_input_ta, set lambda, set lambda_input_ta to lambda(own ta lambda_input)
  post_integrate();
  // set initial lambda_const if required, communicate lambda and lambda_const to neighbours
  comm_forward_lambda();
  // pair_lambda_zone->calculate_lambda is again called as default setup force calculation
  // -> communicate lambda_input_ta again to have an appropriate input for this "force" calculation
  // increase invoked update to prevent double values
  communicate_lambda_input_ta();

  // just to be sure
  write_peratom_stats();
}

/**
  * The new lambda is stored for own atoms in lambda_input_ta since the last force calculation.
  * Update lambda and lambda_const with the running average including the new lambda.
  */

void FixLambdaAPIP::post_integrate()
{
  double *lambda, *lambda_const, *lambda_input_ta;

  lambda = atom->apip_lambda;
  lambda_const = atom->apip_lambda_const;
  lambda_input_ta = atom->apip_lambda_input_ta;

  update_lambda_history();    // update running average of lambda with lambda_input_ta
  get_lambda_average();       // and copy running average to lambda_input_ta

  // use lambda_input_ta to set lambda max
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    lambda_const[i] = lambda[i];
    lambda[i] = lambda_input_ta[i];
    // prevent useless fluctuations in the switching zone just due to atomic fluctuations
    // the values 0 and 1 are always permitted to prevent atoms from keeping a lambda of epsilon forever
    if (fabs(lambda[i] - lambda_const[i]) < min_delta_lambda && lambda[i] != 1 && lambda[i] != 0)
      lambda[i] = lambda_const[i];
  }

  // set lambda_input_ta to own lambda for new lambda calculation with pair_lambda_zone_apip.cpp
  calculate_lambda_input_ta();
}

/**
  * write stats at end of step
  */

void FixLambdaAPIP::end_of_step()
{
  write_peratom_stats();
}

/**
  * The new calculated lambda is stored in lambda.
  * Exchange lambda and lambda_const with neighbours.
  * Use only in setup_pre_force.
  */

void FixLambdaAPIP::comm_forward_lambda()
{
  if (history2_used == 1) {
    // There is only one one calculated lambda.
    // -> This is the first lambda calculation.
    // -> set lambda_const to lambda since there is no previous lambda
    double *lambda, *lambda_const;
    int nlocal;

    lambda = atom->apip_lambda;
    lambda_const = atom->apip_lambda_const;
    nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++) { lambda_const[i] = lambda[i]; }
  }

  comm_forward_flag = FORWARD_MAX;
  comm->forward_comm(this);
}

/**
  * write per atom stats
  */

void FixLambdaAPIP::write_peratom_stats()
{
  if (!peratom_flag) return;

  int i, j, nlocal;
  double **lambda_input_history, **lambda_history;

  nlocal = atom->nlocal;
  lambda_input_history = fixstore->astore;
  lambda_history = fixstore2->astore;

  // grow stats array if required
  if (atom->nmax > nmax_stats) {
    memory->destroy(peratom_stats);
    nmax_stats = atom->nmax;
    memory->create(peratom_stats, nmax_stats, size_peratom_cols, "lambda:peratom_stats");
    array_atom = peratom_stats;
  }

  for (i = 0; i < nlocal; i++) {
    peratom_stats[i][0] = lambda_input_history[i][history_last];    // lambda_input now
    peratom_stats[i][1] =
        lambda_input_history[i][history_length] / history_used;    // lambda_input averaged
    peratom_stats[i][2] =
        lambda_input_history[i][history_length + 1];           // lambda of own ta lambda_input
    peratom_stats[i][3] = lambda_history[i][history2_last];    // lambda max now
    peratom_stats[i][4] =
        lambda_history[i][history2_length] / history2_used;    // lambda max averaged
  }

  if (dump_history_flag) {
    for (i = 0; i < nlocal; i++) {
      // include lambda_input sum -> <= history_length
      for (j = 0; j <= history_length; j++) {
        peratom_stats[i][j + 5] = lambda_input_history[i][j];    // lambda_input history
      }
      for (j = 0; j <= history2_length; j++) {
        peratom_stats[i][j + 6 + history_length] = lambda_history[i][j];    // lambda_input history
      }
    }
  }
}

/**
  * Update running average of lambda_input.
  */

void FixLambdaAPIP::update_lambda_input_history()
{
  if (invoked_history_update == update->ntimestep) return;
  invoked_history_update = update->ntimestep;

  double *lambda_input, **lambda_input_history;
  int *mask;
  int nlocal;

  lambda_input = atom->apip_lambda_input;
  mask = atom->mask;
  lambda_input_history = fixstore->astore;
  nlocal = atom->nlocal;

  // update stats about written values
  history_last = (history_last + 1) % history_length;
  history_used = std::min(history_used + 1, history_length);

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) {
      lambda_input_history[i][history_length + 1] = lambda_non_group;
      continue;
    }

    // lambda_input_history[particle_number][history_number]
    // history_number:
    //     0 - history_length - 1 : single lambda_input of past time step
    //     history_length         : sum of all stored past lambda_input values
    //     history_length + 1     : lambda(time averaged lambda_input)

    // subtract the lambda_input to be overwritten from sum
    lambda_input_history[i][history_length] -= lambda_input_history[i][history_last];
    // store new lambda_input value
    lambda_input_history[i][history_last] = lambda_input[i];
    // add the new lambda_input to sum
    lambda_input_history[i][history_length] += lambda_input_history[i][history_last];

    // calculate lambda of new time average

    if (group_name_complex && (mask[i] & group_bit_complex)) {
      // hard code complex with highest priority
      lambda_input_history[i][history_length + 1] = 0;
    } else if (group_name_simple && (mask[i] & group_bit_simple)) {
      // hard code simple with second highest priority
      lambda_input_history[i][history_length + 1] = 1;
    } else if ((mask[i] & groupbit) &&
               !(group_name_ignore_lambda_input && (mask[i] & group_bit_ignore_lambda_input))) {
      // calculate lambda based on lambda_input
      lambda_input_history[i][history_length + 1] =
          switching_function_poly(lambda_input_history[i][history_length] / history_used);
    } else {
      lambda_input_history[i][history_length + 1] = lambda_non_group;
    }
  }
}

/**
  * Update running average of lambda with lambda_input_ta.
  */

void FixLambdaAPIP::update_lambda_history()
{
  if (invoked_history2_update == update->ntimestep) return;
  invoked_history2_update = update->ntimestep;

  double *lambda_input_ta, **lambda_history;
  int *mask;
  int nlocal;

  lambda_input_ta = atom->apip_lambda_input_ta;
  mask = atom->mask;
  lambda_history = fixstore2->astore;
  nlocal = atom->nlocal;

  // update stats about written values

  history2_last = (history2_last + 1) % history2_length;
  history2_used = std::min(history2_used + 1, history2_length);

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    // lambda_history[particle_number][history2_number]
    // history2_number:
    //     0 - history2_length - 1 : lambda of past time step
    //     history2_length         : sum of all stored past lambda values

    // subtract the lambda to be overwritten from sum
    lambda_history[i][history2_length] -= lambda_history[i][history2_last];
    // store new lambda_input value
    lambda_history[i][history2_last] = lambda_input_ta[i];
    // add the new lambda_input to sum
    lambda_history[i][history2_length] += lambda_history[i][history2_last];
  }
}

/**
  * Copy running average of lambda to lambda_input_ta.
  */

void FixLambdaAPIP::get_lambda_average()
{
  double *lambda_input_ta, **lambda_history, avg;
  int *mask, nlocal;

  mask = atom->mask;
  nlocal = atom->nlocal;
  lambda_history = fixstore2->astore;
  lambda_input_ta = atom->apip_lambda_input_ta;

  // recalculate history sum to limit floating point issues since only changes of the sum are tracked
  if (history2_last == history2_length - 1) {
    double sum;
    for (int i = 0; i < nlocal; i++) {
      sum = 0;
      for (int j = 0; j < history2_length; j++) sum += lambda_history[i][j];
      lambda_history[i][history2_length] = sum;
    }
  }

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      // avg is a division and not exactly 1 or 0
      // exact values are required otherwise many useless complex calculations follow since avg = 1 - epsilon < 1
      // -> hard code zero and one if required
      avg = lambda_history[i][history2_length] / history2_used;
      if (avg > 0.9999)
        lambda_input_ta[i] = 1;    // simple
      else if (avg < 0.0001)
        lambda_input_ta[i] = 0;    // complex
      else
        lambda_input_ta[i] = avg;    // switching
    }
  }
}

/**
  * calculate lambda_input_ta for own atoms
  */

void FixLambdaAPIP::calculate_lambda_input_ta()
{
  int i, nlocal;
  double *lambda_input_ta = atom->apip_lambda_input_ta;

  nlocal = atom->nlocal;

  // store time averaged lambda_input for own atoms
  double **lambda_input_history = fixstore->astore;
  for (i = 0; i < nlocal; i++) lambda_input_ta[i] = lambda_input_history[i][history_length + 1];
}

/* ----------------------------------------------------------------------
  Compute temporary lambda for owned atoms based on lambda_input_history of owned atoms.
  Save this temporary lambda as lambda_input_ta and send it to neighbours.
------------------------------------------------------------------------- */

void FixLambdaAPIP::communicate_lambda_input_ta()
{
  calculate_lambda_input_ta();

  // lambda_input_ta is known only for own atoms
  // -> exchange lambda_input_ta
  comm_forward_flag = FORWARD_TA;
  comm->forward_comm(this);
}

// helper function
// similar to cutoff_func_poly in ace_radial.cpp
// compare Phys Rev Mat 6, 013804 (2022) APPENDIX C: RADIAL AND CUTOFF FUNCTIONS 2. Cutoff function
// the first two derivatives of the switching function lambda vanishes at the boundaries of the switching region
double FixLambdaAPIP::switching_function_poly(double input)
{
  // calculate lambda
  if (input <= threshold_lo) {
    return 1;
  } else if (input >= threshold_hi) {
    return 0;
  } else {
    double deltatmp = 1 - 2 * (1 + (input - threshold_hi) / (threshold_width));
    return 0.5 + 7.5 / 2. * (deltatmp / 4. - pow(deltatmp, 3) / 6. + pow(deltatmp, 5) / 20.);
  }
}

/**
  * Send lambda to neighbours.
  */

int FixLambdaAPIP::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;
  double *lambda_input_ta = atom->apip_lambda_input_ta;
  m = 0;

  if (comm_forward_flag == FORWARD_TA) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = lambda_input_ta[j];
    }
  } else if (comm_forward_flag == FORWARD_MAX) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = atom->apip_lambda[j];
      buf[m++] = atom->apip_lambda_const[j];
    }
  }

  return m;
}

/**
  * Recv lambda from neighbours.
  */

void FixLambdaAPIP::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;
  double *lambda_input_ta = atom->apip_lambda_input_ta;

  m = 0;
  last = first + n;
  if (comm_forward_flag == FORWARD_TA) {
    for (i = first; i < last; i++) { lambda_input_ta[i] = buf[m++]; }
  } else if (comm_forward_flag == FORWARD_MAX) {
    for (i = first; i < last; i++) {
      atom->apip_lambda[i] = buf[m++];
      atom->apip_lambda_const[i] = buf[m++];
    }
  }
}

/**
   * store scalar history information
   */

void FixLambdaAPIP::write_restart(FILE *fp)
{
  bigint timesteps_since_invoked_history_update = update->ntimestep - invoked_history_update;
  bigint timesteps_since_invoked_history2_update = update->ntimestep - invoked_history2_update;

  int n = 0;
  double list[8];
  list[n++] = history_length;
  list[n++] = history_used;
  list[n++] = history_last;
  list[n++] = timesteps_since_invoked_history_update;
  list[n++] = history2_length;
  list[n++] = history2_used;
  list[n++] = history2_last;
  list[n++] = timesteps_since_invoked_history2_update;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(list, sizeof(double), n, fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixLambdaAPIP::restart(char *buf)
{
  int history_length_br, history2_length_br;

  int n = 0;
  auto *list = (double *) buf;

  history_length_br = static_cast<int>(list[n++]);
  history_used = static_cast<int>(list[n++]);
  history_last = static_cast<int>(list[n++]);
  invoked_history_update = update->ntimestep - (static_cast<int>(list[n++]));
  history2_length_br = static_cast<int>(list[n++]);
  history2_used = static_cast<int>(list[n++]);
  history2_last = static_cast<int>(list[n++]);
  invoked_history2_update = update->ntimestep - (static_cast<int>(list[n++]));

  // simple comparisons first
  if (history_length != history_length_br)
    error->all(FLERR, "fix lambda: history_length = {} != {} = history_length_before_restart",
               history_length, history_length_br);
  if (history2_length != history2_length_br)
    error->all(FLERR, "fix lambda: history2_length = {} != {} = history2_length_before_restart",
               history2_length, history2_length_br);
}

/**
  * extract lambda(time averaged lambda_input) and lambda_input_history_len
  */

void *FixLambdaAPIP::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "fix_lambda:lambda_input_history") == 0 && fixstore) { return fixstore->astore; }
  dim = 0;
  if (strcmp(str, "fix_lambda:lambda_input_history_len") == 0) { return &history_length; }
  if (strcmp(str, "fix_lambda:cut_lo") == 0) { return &cut_lo; }
  if (strcmp(str, "fix_lambda:cut_hi") == 0) { return &cut_hi; }
  if (strcmp(str, "fix_lambda:cut_hi_sq") == 0) { return &cut_hi_sq; }
  if (strcmp(str, "fix_lambda:cut_width") == 0) { return &cut_width; }
  if (strcmp(str, "fix_lambda:lambda_non_group") == 0) { return &lambda_non_group; }
  return nullptr;
}
