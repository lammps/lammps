// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <fstream>
#include <iostream>
#include <iomanip>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarvalue.h"
#include "colvarbias_restraint.h"



colvarbias_restraint::colvarbias_restraint(char const *key)
  : colvarbias_ti(key)
  // intermediate in virtual inheritance, no need to initialize base class colvarbias
{}


int colvarbias_restraint::init(std::string const &conf)
{
  state_keyword = "restraint";

  int err = colvarbias::init(conf);
  if (err != COLVARS_OK) {
    return err;
  }
  enable(f_cvb_apply_force);

  colvarbias_ti::init(conf);

  if (cvm::debug())
    cvmodule->log("Initializing a new restraint bias.\n");

  return COLVARS_OK;
}


int colvarbias_restraint::update()
{
  // Update base class (bias_energy and colvar_forces are zeroed there)
  colvarbias::update();

  // Force and energy calculation
  for (size_t i = 0; i < num_variables(); i++) {
    bias_energy += restraint_potential(i);
    colvar_forces[i].type(variables(i)->value());
    colvar_forces[i].is_derivative();
    colvar_forces[i] = restraint_force(i);
  }

  if (cvm::debug())
    cvmodule->log("Done updating the restraint bias \""+this->name+"\".\n");

  if (cvm::debug())
    cvmodule->log("Current forces for the restraint bias \""+
             this->name+"\": "+cvm::to_str(colvar_forces)+".\n");

  return COLVARS_OK;
}


colvarbias_restraint::~colvarbias_restraint()
{
}


std::string const colvarbias_restraint::get_state_params() const
{
  return colvarbias::get_state_params();
}


int colvarbias_restraint::set_state_params(std::string const &conf)
{
  return colvarbias::set_state_params(conf);
}


std::ostream & colvarbias_restraint::write_traj_label(std::ostream &os)
{
  return colvarbias::write_traj_label(os);
}


std::ostream & colvarbias_restraint::write_traj(std::ostream &os)
{
  return colvarbias::write_traj(os);
}



colvarbias_restraint_centers::colvarbias_restraint_centers(char const *key)
  : colvarbias_ti(key)
{
}


int colvarbias_restraint_centers::init(std::string const &conf)
{
  size_t i;

  bool null_centers = (colvar_centers.size() == 0);
  if (null_centers) {
    // try to initialize the restraint centers for the first time
    colvar_centers.resize(num_variables());
    for (i = 0; i < num_variables(); i++) {
      colvar_centers[i].type(variables(i)->value());
      colvar_centers[i].reset();
    }
  }

  if (get_keyval(conf, "centers", colvar_centers, colvar_centers)) {
    for (i = 0; i < num_variables(); i++) {
      if (cvm::debug()) {
        cvmodule->log("colvarbias_restraint: parsing initial centers, i = "+cvm::to_str(i)+".\n");
      }
      colvar_centers[i].apply_constraints();
    }
    null_centers = false;
  }

  if (null_centers) {
    colvar_centers.clear();
    cvmodule->error("Error: must define the initial centers of the restraints.\n", COLVARS_INPUT_ERROR);
    return COLVARS_INPUT_ERROR;
  }

  if (colvar_centers.size() != num_variables()) {
    cvmodule->error("Error: number of centers does not match "
               "that of collective variables.\n", COLVARS_INPUT_ERROR);
    return COLVARS_INPUT_ERROR;
  }

  return COLVARS_OK;
}


int colvarbias_restraint_centers::change_configuration(std::string const &conf)
{
  if (get_keyval(conf, "centers", colvar_centers, colvar_centers)) {
    for (size_t i = 0; i < num_variables(); i++) {
      colvar_centers[i].type(variables(i)->value());
      colvar_centers[i].apply_constraints();
    }
  }
  return COLVARS_OK;
}



colvarbias_restraint_k::colvarbias_restraint_k(char const *key)
  : colvarbias_ti(key), colvarbias_restraint(key)
{
  force_k = -1.0;
  check_positive_k = true;
}


int colvarbias_restraint_k::init(std::string const &conf)
{
  get_keyval(conf, "forceConstant", force_k, (force_k > 0.0 ? force_k : 1.0));
  if (check_positive_k && (force_k < 0.0)) {
    cvmodule->error("Error: undefined or invalid force constant.\n", COLVARS_INPUT_ERROR);
    return COLVARS_INPUT_ERROR;
  }
  return COLVARS_OK;
}


int colvarbias_restraint_k::change_configuration(std::string const &conf)
{
  get_keyval(conf, "forceConstant", force_k, force_k);
  return COLVARS_OK;
}


int colvarbias_restraint_moving::init(std::string const &conf)
{
  if (b_chg_centers && b_chg_force_k) {
    cvmodule->error("Error: cannot specify both targetCenters and targetForceConstant.\n",
               COLVARS_INPUT_ERROR);
    return COLVARS_INPUT_ERROR;
  }
  if (b_chg_centers && b_chg_walls) {
    cvmodule->error("Error: cannot specify both targetCenters and targetUpperWalls/targetLowerWalls.\n",
               COLVARS_INPUT_ERROR);
    return COLVARS_INPUT_ERROR;
  }

  if (b_chg_centers || b_chg_force_k || b_chg_walls) {

    first_step = cvmodule->step_absolute();

    cvmodule->log("Initial step for restraint change: " + cvm::to_str(first_step) + "\n");

    get_keyval(conf, "targetNumSteps", target_nsteps, target_nsteps);
    if (!target_nsteps) {
      return cvmodule->error("Error: targetNumSteps must be non-zero.\n", COLVARS_INPUT_ERROR);
    }

    get_keyval(conf, "targetEquilSteps", target_equil_steps, target_equil_steps);

    if ( target_nsteps && target_nsteps % time_step_factor ) {
      return cvmodule->error("Error: targetNumSteps must be a multiple of timeStepFactor.\n", COLVARS_INPUT_ERROR);
    }

    if (get_keyval(conf, "targetNumStages", target_nstages, target_nstages) &&
        lambda_schedule.size()) {
      return cvmodule->error("Error: targetNumStages and lambdaSchedule are incompatible.\n", COLVARS_INPUT_ERROR);
    }

    get_keyval_feature(this, conf, "outputAccumulatedWork",
                       f_cvb_output_acc_work,
                       is_enabled(f_cvb_output_acc_work));
    if (is_enabled(f_cvb_output_acc_work) && (target_nstages > 0)) {
      return cvmodule->error("Error: outputAccumulatedWork and targetNumStages "
                        "are incompatible.\n", COLVARS_INPUT_ERROR);
    }
  }

  return COLVARS_OK;
}


int colvarbias_restraint_moving::update() {
  if (!cvmodule->proxy->simulation_running() ||
      (!b_chg_centers && !b_chg_force_k && !b_chg_walls)) {
    // Nothing to update
    return COLVARS_OK;
  }

  if (target_nstages) {
    // Accumulate free energy derivative at every step except 0
    // (there are n+1 calls to update() in an n-step simulation)
    if (cvmodule->step_absolute() - first_step > 0) {
      if (b_chg_force_k) dA_dlambda += dU_dlambda_k();
      if (b_chg_centers) dA_dlambda += dU_dlambda_centers();
      if (b_chg_walls) dA_dlambda += dU_dlambda_walls();

      if ((cvmodule->step_absolute() - first_step) % target_nsteps == 0) {

        cvm::real lambda = current_lambda();
        dA_dlambda /= cvm::real(target_nsteps - target_equil_steps);
        std::string msg = "Restraint \"" + this->name + "\" end of stage " + cvm::to_str(stage) +
                          " at step " + cvm::to_str(cvmodule->step_absolute()) +
                          ", lambda= " + cvm::to_str(lambda) +
                          " : dA/dlambda= " + cvm::to_str(dA_dlambda);
        cvmodule->log(msg);
        dA_dlambda = 0.0;

        update_stage();
        lambda = current_lambda();
        if (b_chg_force_k) update_k(lambda);
        if (b_chg_centers) update_centers(lambda);
        if (b_chg_walls) update_walls(lambda);
      }
    }
  } else if (target_nsteps && (cvmodule->step_absolute() - first_step <= target_nsteps)) {
    // Continuous update (slow growth)
    cvm::real lambda = current_lambda();
    if (b_chg_force_k) update_k(lambda);
    if (b_chg_centers) update_centers(lambda);
    if (b_chg_walls) update_walls(lambda);
  }

  return COLVARS_OK;
}


cvm::real colvarbias_restraint_moving::current_lambda() const {
  cvm::real lambda = 0.;
  if (target_nstages) {
    if (lambda_schedule.size()) {
      lambda = lambda_schedule[stage];
    } else {
      lambda = cvm::real(stage) / cvm::real(target_nstages);
      if (b_decoupling) lambda = 1.0 - lambda;
    }
  } else if (target_nsteps) {
    lambda = cvm::real(cvmodule->step_absolute() - first_step) / cvm::real(target_nsteps);
    if (lambda > 1.0) lambda = 1.0;
    if (b_decoupling) lambda = 1.0 - lambda;
  }
  return lambda;
}


std::string const colvarbias_restraint_moving::get_state_params() const
{
  std::ostringstream os;
  os.setf(std::ios::scientific, std::ios::floatfield);
  if (b_chg_centers || b_chg_force_k) {
    os << "firstStep " << std::setw(cvmodule->it_width) << first_step << "\n";
    if (target_nstages) {
      os << "stage " << std::setw(cvmodule->it_width) << stage << "\n";
    }
  }
  return os.str();
}


int colvarbias_restraint_moving::set_state_params(std::string const &conf)
{
  if (b_chg_centers || b_chg_force_k) {
    auto first_step_flags = colvarparse::parse_restart;
    if (cvmodule->restart_version_number() > 20230906) {
      // Only require the first step when the code could produce it
      first_step_flags = colvarparse::parse_restart | colvarparse::parse_required;
    }
    get_keyval(conf, "firstStep", first_step, first_step, first_step_flags);
    if (target_nstages) {
      get_keyval(conf, "stage", stage, stage,
                 colvarparse::parse_restart | colvarparse::parse_required);
    }
  }
  return COLVARS_OK;
}


colvarbias_restraint_centers_moving::colvarbias_restraint_centers_moving(char const *key)
  : colvarbias_ti(key),
    colvarbias_restraint(key),
    colvarbias_restraint_centers(key),
    colvarbias_restraint_moving(key)
{
  b_chg_centers = false;
  b_output_centers = false;
}


int colvarbias_restraint_centers_moving::init(std::string const &conf)
{
  colvarbias_restraint_centers::init(conf);

  if (cvm::debug()) {
    cvmodule->log("colvarbias_restraint: parsing target centers.\n");
  }

  size_t i;
  if (get_keyval(conf, "targetCenters", target_centers, colvar_centers)) {
    if (target_centers.size() != num_variables()) {
      cvmodule->error("Error: number of target centers does not match "
                 "that of collective variables.\n", COLVARS_INPUT_ERROR);
    }
    b_chg_centers = true;
    for (i = 0; i < target_centers.size(); i++) {
      target_centers[i].apply_constraints();
      centers_incr.push_back(colvar_centers[i]);
      centers_incr[i].reset();
    }
  }

  if (b_chg_centers) {
    // parse moving schedule options
    colvarbias_restraint_moving::init(conf);
    if (initial_centers.size() == 0) {
      // One-time init
      initial_centers = colvar_centers;
    }
    // Call to check that the definition is correct
    for (i = 0; i < num_variables(); i++) {
      colvarvalue const midpoint =
        colvarvalue::interpolate(initial_centers[i],
                                 target_centers[i],
                                 0.5);
    }

  } else {
    target_centers.clear();
  }

  // Output restraint centers even when they do not change; some NAMD REUS
  // scripts expect this behavior
  get_keyval(conf, "outputCenters", b_output_centers, b_output_centers);

  return COLVARS_OK;
}


void colvarbias_restraint_centers_moving::update_centers(cvm::real lambda)
{
  if (cvm::debug()) {
    cvmodule->log("Updating centers for the restraint bias \""+
             this->name+"\": "+cvm::to_str(colvar_centers)+".\n");
  }
  size_t i;
  for (i = 0; i < num_variables(); i++) {
    colvarvalue const c_new = colvarvalue::interpolate(initial_centers[i],
                                                       target_centers[i],
                                                       lambda);
    centers_incr[i] = 0.5 * c_new.dist2_grad(colvar_centers[i]);
    colvar_centers[i] = c_new;
    variables(i)->wrap(colvar_centers[i]);
  }
  if (cvm::debug()) {
    cvmodule->log("New centers for the restraint bias \""+
             this->name+"\": "+cvm::to_str(colvar_centers)+".\n");
  }

  if (!target_nstages && (cvmodule->step_absolute() - first_step > target_nsteps)) {
    for (size_t i = 0; i < num_variables(); i++) {
      centers_incr[i].reset();
    }
  }
  return;
}


cvm::real colvarbias_restraint_centers_moving::dU_dlambda_centers() const
{
  cvm::real grad = 0.0;
  for (size_t i = 0; i < colvar_centers.size(); ++i) {
    grad += restraint_force(i) * (target_centers[i] - initial_centers[i]);
  }
  return grad;
}


int colvarbias_restraint_centers_moving::update_acc_work()
{
  if (!cvmodule->proxy->simulation_running()) {
    return COLVARS_OK;
  }
  if (b_chg_centers && is_enabled(f_cvb_output_acc_work)) {
    if ((cvmodule->step_relative() > 0) &&
        (cvmodule->step_absolute() - first_step <= target_nsteps)) {
      for (size_t i = 0; i < num_variables(); i++) {
        // project forces on the calculated increments at this step
        acc_work += colvar_forces[i] * centers_incr[i];
      }
    }
  }
  return COLVARS_OK;
}


std::string const colvarbias_restraint_centers_moving::get_state_params() const
{
  std::ostringstream os;
  os.setf(std::ios::scientific, std::ios::floatfield);

  if (b_chg_centers) {
    size_t i;
    os << "centers ";
    for (i = 0; i < num_variables(); i++) {
      os << " "
         << std::setprecision(cvmodule->cv_prec) << std::setw(cvmodule->cv_width)
         << colvar_centers[i];
    }
    os << "\n";

    if (is_enabled(f_cvb_output_acc_work)) {
      os << "accumulatedWork "
         << std::setprecision(cvmodule->en_prec) << std::setw(cvmodule->en_width)
         << acc_work << "\n";
    }
  }

  return os.str();
}


int colvarbias_restraint_centers_moving::set_state_params(std::string const &conf)
{
  colvarbias_restraint::set_state_params(conf);

  if (b_chg_centers) {
    get_keyval(conf, "centers", colvar_centers, colvar_centers,
               colvarparse::parse_restart | colvarparse::parse_required);
  }

  if (is_enabled(f_cvb_output_acc_work)) {
    get_keyval(conf, "accumulatedWork", acc_work, acc_work,
               colvarparse::parse_restart | colvarparse::parse_required);
  }

  return COLVARS_OK;
}


std::ostream & colvarbias_restraint_centers_moving::write_traj_label(std::ostream &os)
{
  if (b_output_centers) {
    for (size_t i = 0; i < num_variables(); i++) {
      size_t const this_cv_width = (variables(i)->value()).output_width(cvmodule->cv_width);
      os << " x0_"
         << cvmodule->wrap_string(variables(i)->name, this_cv_width-3);
    }
  }

  if (b_chg_centers && is_enabled(f_cvb_output_acc_work)) {
    os << " W_"
       << cvmodule->wrap_string(this->name, cvmodule->en_width-2);
  }

  return os;
}


std::ostream & colvarbias_restraint_centers_moving::write_traj(std::ostream &os)
{
  if (b_output_centers) {
    for (size_t i = 0; i < num_variables(); i++) {
      os << " "
         << std::setprecision(cvmodule->cv_prec) << std::setw(cvmodule->cv_width)
         << colvar_centers[i];
    }
  }

  if (b_chg_centers && is_enabled(f_cvb_output_acc_work)) {
    os << " "
       << std::setprecision(cvmodule->en_prec) << std::setw(cvmodule->en_width)
       << acc_work;
  }

  return os;
}



colvarbias_restraint_k_moving::colvarbias_restraint_k_moving(char const *key)
  : colvarbias_ti(key),
    colvarbias_restraint(key),
    colvarbias_restraint_k(key),
    colvarbias_restraint_moving(key)
{
  b_chg_force_k = false;
  b_decoupling = false;
  target_force_k = -1.0;
  starting_force_k = -1.0;
  lambda_exp = 1.0;
  force_k_incr = 0.0;
}


int colvarbias_restraint_k_moving::init(std::string const &conf)
{
  colvarbias_restraint_k::init(conf);

  get_keyval(conf, "decoupling", b_decoupling, b_decoupling);
  if (b_decoupling) {
    target_force_k = force_k;
    starting_force_k = 0.0;
    b_chg_force_k = true;
  }

  if (get_keyval(conf, "targetForceConstant", target_force_k, target_force_k)) {
    if (b_decoupling) {
      cvmodule->error("Error: targetForceConstant may not be specified together with decoupling.\n", COLVARS_INPUT_ERROR);
      return COLVARS_ERROR;
    }
    b_chg_force_k = true;
    starting_force_k = force_k;
  }

  if (!b_chg_force_k) {
    return COLVARS_OK;
  }

  // parse moving restraint options
  colvarbias_restraint_moving::init(conf);

  if (get_keyval(conf, "lambdaSchedule", lambda_schedule, lambda_schedule) &&
      target_nstages > 0) {
    cvmodule->error("Error: targetNumStages and lambdaSchedule are incompatible.\n", COLVARS_INPUT_ERROR);
    return cvmodule->get_error();
  }

  if (lambda_schedule.size()) {
    // There is one more lambda-point than stages
    target_nstages = lambda_schedule.size() - 1;
  }

  if ((get_keyval(conf, "targetForceExponent", lambda_exp, lambda_exp, parse_deprecated)
    || get_keyval(conf, "lambdaExponent", lambda_exp, lambda_exp))
    && !b_chg_force_k) {
    cvmodule->error("Error: cannot set lambdaExponent unless a changing force constant is active.\n", COLVARS_INPUT_ERROR);
  }
  if (lambda_exp < 1.0) {
    cvmodule->log("Warning: for all practical purposes, lambdaExponent should be 1.0 or greater.\n");
  }

  return COLVARS_OK;
}


void colvarbias_restraint_k_moving::update_k(cvm::real lambda) {
  cvm::real const force_k_old = force_k;

  force_k = starting_force_k + (target_force_k - starting_force_k) * cvm::pow(lambda, lambda_exp);
  force_k_incr = force_k - force_k_old;
  if (!target_nstages && (cvmodule->step_absolute() > first_step + target_nsteps)) {
    force_k_incr = 0.0;
  }
  cvmodule->log("Updated force constant for the restraint bias \""+
           this->name+"\": "+cvm::to_str(force_k)+".\n");
}


cvm::real colvarbias_restraint_k_moving::dU_dlambda_k() const {
  if (force_k == 0.0) return 0.0;
  cvm::real dU_dk = 0.0;
  for (size_t i = 0; i < num_variables(); i++) {
    dU_dk += d_restraint_potential_dk(i);
  }
  return dU_dk * (target_force_k - starting_force_k) * lambda_exp * cvm::pow(current_lambda(), lambda_exp - 1.0);
}


int colvarbias_restraint_k_moving::update_acc_work()
{
  if (!cvmodule->proxy->simulation_running()) {
    return COLVARS_OK;
  }
  if (b_chg_force_k && is_enabled(f_cvb_output_acc_work)) {
    if ((cvmodule->step_relative() > 0) &&
        (cvmodule->step_absolute() - first_step <= target_nsteps)) {
      cvm::real dU_dk = 0.0;
      for (size_t i = 0; i < num_variables(); i++) {
        dU_dk += d_restraint_potential_dk(i);
      }
      acc_work += dU_dk * force_k_incr;
    }
  }
  return COLVARS_OK;
}


std::string const colvarbias_restraint_k_moving::get_state_params() const
{
  std::ostringstream os;
  os.setf(std::ios::scientific, std::ios::floatfield);
  if (b_chg_force_k) {
    os << "forceConstant "
       << std::setprecision(cvmodule->en_prec)
       << std::setw(cvmodule->en_width) << force_k << "\n";

    if (is_enabled(f_cvb_output_acc_work)) {
      os << "accumulatedWork "
         << std::setprecision(cvmodule->en_prec) << std::setw(cvmodule->en_width)
         << acc_work << "\n";
    }
  }
  return os.str();
}


int colvarbias_restraint_k_moving::set_state_params(std::string const &conf)
{
  colvarbias_restraint::set_state_params(conf);

  if (b_chg_force_k) {
    get_keyval(conf, "forceConstant", force_k, force_k,
               colvarparse::parse_restart | colvarparse::parse_required);
  }

  if (is_enabled(f_cvb_output_acc_work)) {
    get_keyval(conf, "accumulatedWork", acc_work, acc_work,
               colvarparse::parse_restart | colvarparse::parse_required);
  }

  return COLVARS_OK;
}


std::ostream & colvarbias_restraint_k_moving::write_traj_label(std::ostream &os)
{
  if (b_chg_force_k && is_enabled(f_cvb_output_acc_work)) {
    os << " W_"
       << cvmodule->wrap_string(this->name, cvmodule->en_width-2);
  }
  return os;
}


std::ostream & colvarbias_restraint_k_moving::write_traj(std::ostream &os)
{
  if (b_chg_force_k && is_enabled(f_cvb_output_acc_work)) {
    os << " "
       << std::setprecision(cvmodule->en_prec) << std::setw(cvmodule->en_width)
       << acc_work;
  }
  return os;
}



colvarbias_restraint_harmonic::colvarbias_restraint_harmonic(colvarmodule *cvmodule_in, char const *key)
  : colvarbias(cvmodule_in, key),
    colvarbias_ti(key),
    colvarbias_restraint(key),
    colvarbias_restraint_centers(key),
    colvarbias_restraint_moving(key),
    colvarbias_restraint_k(key),
    colvarbias_restraint_centers_moving(key),
    colvarbias_restraint_k_moving(key)
{
}


int colvarbias_restraint_harmonic::init(std::string const &conf)
{
  colvarbias_restraint::init(conf);
  colvarbias_restraint_moving::init(conf);
  colvarbias_restraint_centers_moving::init(conf);
  colvarbias_restraint_k_moving::init(conf);

  cvmodule->cite_feature("Harmonic colvar bias implementation");

  for (size_t i = 0; i < num_variables(); i++) {
    cvm::real const w = variables(i)->width;
    cvmodule->log("The force constant for colvar \""+variables(i)->name+
             "\" will be rescaled to "+
             cvm::to_str(force_k/(w*w))+
             " according to the specified width ("+cvm::to_str(w)+").\n");
  }

  return COLVARS_OK;
}


int colvarbias_restraint_harmonic::update()
{
  int error_code = COLVARS_OK;

  // update the TI estimator (if defined)
  error_code |= colvarbias_ti::update();

  // update parameters (centers or force constant)
  error_code |= colvarbias_restraint_moving::update();

  // update restraint energy and forces
  error_code |= colvarbias_restraint::update();

  // update accumulated work using the current forces
  error_code |= colvarbias_restraint_centers_moving::update_acc_work();
  error_code |= colvarbias_restraint_k_moving::update_acc_work();

  return error_code;
}


cvm::real colvarbias_restraint_harmonic::restraint_potential(size_t i) const
{
  return 0.5 * force_k / (variables(i)->width * variables(i)->width) *
    variables(i)->dist2(variables(i)->value(), colvar_centers[i]);
}


colvarvalue const colvarbias_restraint_harmonic::restraint_force(size_t i) const
{
  return -0.5 * force_k / (variables(i)->width * variables(i)->width) *
    variables(i)->dist2_lgrad(variables(i)->value(), colvar_centers[i]);
}


cvm::real colvarbias_restraint_harmonic::d_restraint_potential_dk(size_t i) const
{
  return 0.5 / (variables(i)->width * variables(i)->width) *
    variables(i)->dist2(variables(i)->value(), colvar_centers[i]);
}


std::string const colvarbias_restraint_harmonic::get_state_params() const
{
  return colvarbias_restraint::get_state_params() +
    colvarbias_restraint_moving::get_state_params() +
    colvarbias_restraint_centers_moving::get_state_params() +
    colvarbias_restraint_k_moving::get_state_params();
}


int colvarbias_restraint_harmonic::set_state_params(std::string const &conf)
{
  int error_code = COLVARS_OK;
  error_code |= colvarbias_restraint::set_state_params(conf);
  error_code |= colvarbias_restraint_moving::set_state_params(conf);
  error_code |= colvarbias_restraint_centers_moving::set_state_params(conf);
  error_code |= colvarbias_restraint_k_moving::set_state_params(conf);
  return error_code;
}


std::ostream & colvarbias_restraint_harmonic::write_traj_label(std::ostream &os)
{
  colvarbias_restraint::write_traj_label(os);
  colvarbias_restraint_centers_moving::write_traj_label(os);
  colvarbias_restraint_k_moving::write_traj_label(os);
  return os;
}


std::ostream & colvarbias_restraint_harmonic::write_traj(std::ostream &os)
{
  colvarbias_restraint::write_traj(os);
  colvarbias_restraint_centers_moving::write_traj(os);
  colvarbias_restraint_k_moving::write_traj(os);
  return os;
}


int colvarbias_restraint_harmonic::change_configuration(std::string const &conf)
{
  cvmodule->log("Warning: Updating configuration mid-run for an existing restraint is unsupported, and will raise an error in a future release.  The supported method is deleting the restraint and re-creating, which will also do a better error check.");
  return colvarbias_restraint_centers::change_configuration(conf) |
    colvarbias_restraint_k::change_configuration(conf);
}


cvm::real colvarbias_restraint_harmonic::energy_difference(std::string const &conf)
{
  cvm::real const old_bias_energy = bias_energy;
  cvm::real const old_force_k = force_k;
  std::vector<colvarvalue> const old_centers = colvar_centers;

  change_configuration(conf);
  update();

  cvm::real const result = (bias_energy - old_bias_energy);

  bias_energy = old_bias_energy;
  force_k = old_force_k;
  colvar_centers = old_centers;

  return result;
}



colvarbias_restraint_harmonic_walls::colvarbias_restraint_harmonic_walls(colvarmodule *cvmodule_in, char const *key)
  : colvarbias(cvmodule_in, key),
    colvarbias_ti(key),
    colvarbias_restraint(key),
    colvarbias_restraint_k(key),
    colvarbias_restraint_moving(key),
    // colvarbias_restraint_centers_moving(key)
    colvarbias_restraint_k_moving(key)
{
  lower_wall_k = -1.0;
  upper_wall_k = -1.0;
  // This bias implements the bias_actual_colvars feature (most others do not)
  provide(f_cvb_bypass_ext_lagrangian);
  set_enabled(f_cvb_bypass_ext_lagrangian); // Defaults to enabled
}


int colvarbias_restraint_harmonic_walls::init(std::string const &conf)
{


get_keyval(conf, "targetLowerWalls", target_lower_walls, target_lower_walls);
get_keyval(conf, "targetUpperWalls", target_upper_walls, target_upper_walls);

  if ((target_lower_walls.size() > 0) || (target_upper_walls.size() > 0)) {
    b_chg_walls=true;
  }

  colvarbias_restraint::init(conf);
  colvarbias_restraint_moving::init(conf);
  colvarbias_restraint_k_moving::init(conf);

  cvmodule->cite_feature("harmonicWalls colvar bias implementation");

  enable(f_cvb_scalar_variables);

  size_t i;

  bool b_null_lower_walls = false;
  if (lower_walls.size() == 0) {
    b_null_lower_walls = true;
    lower_walls.resize(num_variables());
    for (i = 0; i < num_variables(); i++) {
      lower_walls[i].type(variables(i)->value());
      lower_walls[i].reset();
    }
  }
  if (!get_keyval(conf, "lowerWalls", lower_walls, lower_walls) &&
      b_null_lower_walls) {
    cvmodule->log("Lower walls were not provided.\n");
    lower_walls.clear();
  }

  bool b_null_upper_walls = false;
  if (upper_walls.size() == 0) {
    b_null_upper_walls = true;
    upper_walls.resize(num_variables());
    for (i = 0; i < num_variables(); i++) {
      upper_walls[i].type(variables(i)->value());
      upper_walls[i].reset();
    }
  }
  if (!get_keyval(conf, "upperWalls", upper_walls, upper_walls) &&
      b_null_upper_walls) {
    cvmodule->log("Upper walls were not provided.\n");
    upper_walls.clear();
  }

  if ((lower_walls.size() == 0) && (upper_walls.size() == 0)) {
    return cvmodule->error("Error: no walls provided.\n", COLVARS_INPUT_ERROR);
  }

  if (lower_walls.size() > 0) {
    get_keyval(conf, "lowerWallConstant", lower_wall_k,
               (lower_wall_k > 0.0) ? lower_wall_k : force_k);
  }
  if (upper_walls.size() > 0) {
    get_keyval(conf, "upperWallConstant", upper_wall_k,
               (upper_wall_k > 0.0) ? upper_wall_k : force_k);
  }

  if ((lower_walls.size() == 0) || (upper_walls.size() == 0)) {
    for (i = 0; i < num_variables(); i++) {
      if (variables(i)->is_enabled(f_cv_periodic)) {
        return cvmodule->error("Error: at least one variable is periodic, "
                          "both walls must be provided.\n", COLVARS_INPUT_ERROR);
      }
    }
  }

bool b_null_target_lower_walls = false;
  if (target_lower_walls.size() == 0) {
    b_null_target_lower_walls = true;
    target_lower_walls.resize(num_variables());
    for (i = 0; i < num_variables(); i++) {
      target_lower_walls[i].type(variables(i)->value());
      target_lower_walls[i].reset();
    }
  }
  if (!get_keyval(conf, "targetLowerWalls", target_lower_walls, target_lower_walls) &&
      b_null_target_lower_walls) {
    cvmodule->log("Target lower walls were not provided.\n");
    target_lower_walls.clear();
  }

  bool b_null_target_upper_walls = false;
  if (target_upper_walls.size() == 0) {
    b_null_target_upper_walls = true;
    target_upper_walls.resize(num_variables());
    for (i = 0; i < num_variables(); i++) {
      target_upper_walls[i].type(variables(i)->value());
      target_upper_walls[i].reset();
    }
  }
  if (!get_keyval(conf, "targetUpperWalls", target_upper_walls, target_upper_walls) &&
      b_null_target_upper_walls) {
    cvmodule->log("Target upper walls were not provided.\n");
    target_upper_walls.clear();
  }

  if (target_upper_walls.size() > upper_walls.size()) {
    cvmodule->error("At least 1 target wall was provided with no initial wall!.\n",
      COLVARS_INPUT_ERROR);
  }

    if (target_lower_walls.size() > lower_walls.size()) {
    cvmodule->error("At least 1 target wall was provided with no initial wall!.\n",
      COLVARS_INPUT_ERROR);
  }

  if ((lower_walls.size() > 0) && (upper_walls.size() > 0)) {
    for (i = 0; i < num_variables(); i++) {
      if (lower_walls[i] >= upper_walls[i]) {
        return cvmodule->error("Error: one upper wall, "+
                          cvm::to_str(upper_walls[i])+
                          ", is not higher than the lower wall, "+
                          cvm::to_str(lower_walls[i])+".\n",
                          COLVARS_INPUT_ERROR);
      }
      if (variables(i)->dist2(lower_walls[i], upper_walls[i]) < 1.0e-12) {
        return cvmodule->error("Error: lower wall and upper wall are equal "
                          "in the domain of the variable \""+
                          variables(i)->name+"\".\n", COLVARS_INPUT_ERROR);
      }
    }
    if ((target_lower_walls.size() > 0) && (target_upper_walls.size() > 0)){
      for (i = 0; i < num_variables(); i++) {
        if (target_lower_walls[i] >= target_upper_walls[i]) {
          return cvmodule->error("Error: one target upper wall, "+
                            cvm::to_str(target_upper_walls[i])+
                            ", will be lower then the target lower wall, "+
                            cvm::to_str(target_lower_walls[i])+".\n",
                            COLVARS_INPUT_ERROR);
      }
    }
  }
    if ((target_lower_walls.size() > 0) && !(target_upper_walls.size() > 0)){
      for (i = 0; i < num_variables(); i++) {
        if (target_lower_walls[i] >= upper_walls[i]) {
          return cvmodule->error("Error: one upper wall, "+
                            cvm::to_str(upper_walls[i])+
                            ", will be surpressed by the target lower wall, "+
                            cvm::to_str(target_lower_walls[i])+".\n",
                            COLVARS_INPUT_ERROR);
      }
    }
  }

    if (!(target_lower_walls.size() > 0) && (target_upper_walls.size() > 0)){
      for (i = 0; i < num_variables(); i++) {
        if (lower_walls[i] >= target_upper_walls[i]) {
          return cvmodule->error("Error: one target upper wall, "+
                            cvm::to_str(target_upper_walls[i])+
                            ", will surpress the lower wall, "+
                            cvm::to_str(lower_walls[i])+".\n",
                            COLVARS_INPUT_ERROR);
      }
    }
  }

    if (lower_wall_k * upper_wall_k == 0.0 && !b_chg_force_k) {
      cvmodule->error("Error: lowerWallConstant and upperWallConstant, "
                 "when defined, must both be positive.\n",
                 COLVARS_INPUT_ERROR);
      return COLVARS_INPUT_ERROR;
    }
    force_k = cvm::sqrt(lower_wall_k * upper_wall_k);
    // transform the two constants to relative values using geometric mean as ref
    // to preserve force_k if provided as single parameter
    // (allow changing both via force_k)
    lower_wall_k /= force_k;
    upper_wall_k /= force_k;
  } else {
    // If only one wall is defined, need to rescale as well
    if (lower_walls.size() > 0) {
      force_k = lower_wall_k;
      lower_wall_k = 1.0;
    }
    if (upper_walls.size() > 0) {
      force_k = upper_wall_k;
      upper_wall_k = 1.0;
    }
  }

  // Initialize starting value of the force constant (in case it's changing)
  starting_force_k = (b_decoupling ? 0.0 : force_k);

  if (lower_walls.size() > 0) {
    for (i = 0; i < num_variables(); i++) {
      cvm::real const w = variables(i)->width;
      cvmodule->log("The lower wall force constant for colvar \""+
               variables(i)->name+"\" will be rescaled to "+
               cvm::to_str(lower_wall_k * force_k / (w*w))+
               " according to the specified width ("+cvm::to_str(w)+").\n");
    }
  }

  if (upper_walls.size() > 0) {
    for (i = 0; i < num_variables(); i++) {
      cvm::real const w = variables(i)->width;
      cvmodule->log("The upper wall force constant for colvar \""+
               variables(i)->name+"\" will be rescaled to "+
               cvm::to_str(upper_wall_k * force_k / (w*w))+
               " according to the specified width ("+cvm::to_str(w)+").\n");
    }
  }

//Initializing starting value for moving walls
if (target_upper_walls.size() > 0) {
  if (initial_upper_walls.size() == 0) {
    // One-time init
    initial_upper_walls = upper_walls;
    upper_walls_incr.resize(target_upper_walls.size());
  }
}

if (target_lower_walls.size() > 0) {
  if (initial_lower_walls.size() == 0) {
    // One-time init
    initial_lower_walls = lower_walls;
    lower_walls_incr.resize(target_lower_walls.size());
  }
}


  return COLVARS_OK;
}

void colvarbias_restraint_harmonic_walls::update_walls(cvm::real lambda)
{
  if (target_upper_walls.size() > 0){

    if (cvm::debug()) {
      cvmodule->log("Updating upper walls for the restraint bias \""+
              this->name+"\": "+cvm::to_str(lower_walls)+".\n");
    }

    size_t i;
    for (i = 0; i < num_variables(); i++) {
      colvarvalue const c_new = colvarvalue::interpolate(initial_upper_walls[i],
                                                        target_upper_walls[i],
                                                        lambda);

      upper_walls_incr[i] = 0.5 * c_new.dist2_grad(upper_walls[i]);
      upper_walls[i] = c_new;
      variables(i)->wrap(upper_walls[i]);
    }
    if (cvm::debug()) {
      cvmodule->log("New upper walls for the restraint bias \""+
              this->name+"\": "+cvm::to_str(upper_walls)+".\n");
    }
  }

  if (target_lower_walls.size() > 0){
    if (cvm::debug()) {
      cvmodule->log("Updating lower walls for the restraint bias \""+
              this->name+"\": "+cvm::to_str(lower_walls)+".\n");
    }

    size_t i;
    for (i = 0; i < num_variables(); i++) {
      colvarvalue const c_new = colvarvalue::interpolate(initial_lower_walls[i],
                                                        target_lower_walls[i],
                                                        lambda);
      lower_walls_incr[i] = 0.5 * c_new.dist2_grad(lower_walls[i]);
      lower_walls[i] = c_new;
      variables(i)->wrap(lower_walls[i]);
    }
    if (cvm::debug()) {
      cvmodule->log("New lower walls for the restraint bias \""+
              this->name+"\": "+cvm::to_str(lower_walls)+".\n");
    }
  }
  if (!target_nstages && (cvmodule->step_absolute() - first_step > target_nsteps)) {
    for (size_t i = 0; i < num_variables(); i++) {
        if(lower_walls_incr.size() > 0) {
          lower_walls_incr[i].reset();
        }
        if(upper_walls_incr.size() > 0) {
          upper_walls_incr[i].reset();
        }
    }
  }
  return;
}


cvm::real colvarbias_restraint_harmonic_walls::dU_dlambda_walls() const {
  cvm::real dU_dwall = 0.0;
  for (size_t i = 0; i < num_variables(); i++) {
    cvm::real const dist = colvar_distance(i);
    if (dist > 0.0) {
        dU_dwall += force_k * upper_wall_k * dist/(variables(i)->width);
    } else{
        dU_dwall += force_k * lower_wall_k * dist/(variables(i)->width);
    }
  }
  return dU_dwall;
}


int colvarbias_restraint_harmonic_walls::update()
{
  int error_code = COLVARS_OK;

  error_code |= colvarbias_ti::update();

  error_code |= colvarbias_restraint_moving::update();

  error_code |= colvarbias_restraint::update();

  error_code |= colvarbias_restraint_k_moving::update_acc_work();
  error_code |= colvarbias_restraint_harmonic_walls::update_acc_work();

  return error_code;
}


int colvarbias_restraint_harmonic_walls::update_acc_work()
{
  if (!cvmodule->proxy->simulation_running()) {
    return COLVARS_OK;
  }
  if (b_chg_walls && is_enabled(f_cvb_output_acc_work)) {
    if ((cvmodule->step_relative() > 0) &&
        (cvmodule->step_absolute() - first_step <= target_nsteps)) {
      for (size_t i = 0; i < num_variables(); i++) {
        cvm::real const dist = colvar_distance(i);
        if (dist > 0.0) {
          acc_work += colvar_forces[i] * upper_walls_incr[i];
        } else{
          acc_work += colvar_forces[i] * lower_walls_incr[i];
        }
      }
    }
  }
  return COLVARS_OK;
}


cvm::real colvarbias_restraint_harmonic_walls::colvar_distance(size_t i) const
{
  colvar *cv = variables(i);

  colvarvalue const &cvv = is_enabled(f_cvb_bypass_ext_lagrangian) ?
    variables(i)->actual_value() :
    variables(i)->value();

  // For a periodic colvar, both walls may be applicable at the same time
  // in which case we pick the closer one

  if (cv->is_enabled(f_cv_periodic)) {
    cvm::real const lower_wall_dist2 = cv->dist2(cvv, lower_walls[i]);
    cvm::real const upper_wall_dist2 = cv->dist2(cvv, upper_walls[i]);
    if (lower_wall_dist2 < upper_wall_dist2) {
      cvm::real const grad = cv->dist2_lgrad(cvv, lower_walls[i]);
      if (grad < 0.0) { return 0.5 * grad; }
    } else {
      cvm::real const grad = cv->dist2_lgrad(cvv, upper_walls[i]);
      if (grad > 0.0) { return 0.5 * grad; }
    }
    return 0.0;
  }

  if (lower_walls.size() > 0) {
    cvm::real const grad = cv->dist2_lgrad(cvv, lower_walls[i]);
    if (grad < 0.0) { return 0.5 * grad; }
  }
  if (upper_walls.size() > 0) {
    cvm::real const grad = cv->dist2_lgrad(cvv, upper_walls[i]);
    if (grad > 0.0) { return 0.5 * grad; }
  }
  return 0.0;
}


cvm::real colvarbias_restraint_harmonic_walls::restraint_potential(size_t i) const
{
  cvm::real const dist = colvar_distance(i);
  cvm::real const scale = dist > 0.0 ? upper_wall_k : lower_wall_k;
  return 0.5 * force_k * scale / (variables(i)->width * variables(i)->width) *
    dist * dist;
}


colvarvalue const colvarbias_restraint_harmonic_walls::restraint_force(size_t i) const
{
  cvm::real const dist = colvar_distance(i);
  cvm::real const scale = dist > 0.0 ? upper_wall_k : lower_wall_k;
  return - force_k * scale / (variables(i)->width * variables(i)->width) * dist;
}


cvm::real colvarbias_restraint_harmonic_walls::d_restraint_potential_dk(size_t i) const
{
  cvm::real const dist = colvar_distance(i);
  cvm::real const scale = dist > 0.0 ? upper_wall_k : lower_wall_k;
  return 0.5 * scale / (variables(i)->width * variables(i)->width) *
    dist * dist;
}


std::string const colvarbias_restraint_harmonic_walls::get_state_params() const
{
  return colvarbias_restraint::get_state_params() +
    colvarbias_restraint_moving::get_state_params() +
    colvarbias_restraint_k_moving::get_state_params();
}


int colvarbias_restraint_harmonic_walls::set_state_params(std::string const &conf)
{
  int error_code = COLVARS_OK;
  error_code |= colvarbias_restraint::set_state_params(conf);
  error_code |= colvarbias_restraint_moving::set_state_params(conf);
  error_code |= colvarbias_restraint_k_moving::set_state_params(conf);
  return error_code;
}

int colvarbias_restraint_harmonic_walls::change_configuration(std::string const &conf)
{
  return cvmodule->error("Updating configuration mid-run is not supported for harmonicWalls.",
                    COLVARS_NOT_IMPLEMENTED);
}

std::ostream & colvarbias_restraint_harmonic_walls::write_traj_label(std::ostream &os)
{
  colvarbias_restraint::write_traj_label(os);
  colvarbias_restraint_k_moving::write_traj_label(os);
  return os;
}


std::ostream & colvarbias_restraint_harmonic_walls::write_traj(std::ostream &os)
{
  colvarbias_restraint::write_traj(os);
  colvarbias_restraint_k_moving::write_traj(os);
  return os;
}



colvarbias_restraint_linear::colvarbias_restraint_linear(colvarmodule *cvmodule_in, char const *key)
  : colvarbias(cvmodule_in, key),
    colvarbias_ti(key),
    colvarbias_restraint(key),
    colvarbias_restraint_centers(key),
    colvarbias_restraint_moving(key),
    colvarbias_restraint_k(key),
    colvarbias_restraint_centers_moving(key),
    colvarbias_restraint_k_moving(key)
{
  check_positive_k = false;
}


int colvarbias_restraint_linear::init(std::string const &conf)
{
  colvarbias_restraint::init(conf);
  colvarbias_restraint_moving::init(conf);
  colvarbias_restraint_centers_moving::init(conf);
  colvarbias_restraint_k_moving::init(conf);

  cvmodule->cite_feature("harmonicWalls colvar bias implementation");

  for (size_t i = 0; i < num_variables(); i++) {
    if (variables(i)->is_enabled(f_cv_periodic)) {
      cvmodule->error("Error: linear biases cannot be applied to periodic variables.\n",
                 COLVARS_INPUT_ERROR);
      return COLVARS_INPUT_ERROR;
    }
    cvm::real const w = variables(i)->width;
    cvmodule->log("The force constant for colvar \""+variables(i)->name+
             "\" will be rescaled to "+
             cvm::to_str(force_k / w)+
             " according to the specified width ("+cvm::to_str(w)+").\n");
  }

  return COLVARS_OK;
}


int colvarbias_restraint_linear::update()
{
  int error_code = COLVARS_OK;

  // update the TI estimator (if defined)
  error_code |= colvarbias_ti::update();

  // update parameters (centers or force constant)
  error_code |= colvarbias_restraint_moving::update();

  // update restraint energy and forces
  error_code |= colvarbias_restraint::update();

  // update accumulated work using the current forces
  error_code |= colvarbias_restraint_centers_moving::update_acc_work();
  error_code |= colvarbias_restraint_k_moving::update_acc_work();

  return error_code;
}


int colvarbias_restraint_linear::change_configuration(std::string const &conf)
{
  return cvmodule->error("Updating configuration mid-run is not supported for a linear bias.",
                    COLVARS_NOT_IMPLEMENTED);
}


cvm::real colvarbias_restraint_linear::energy_difference(std::string const &conf)
{
  cvm::real const old_bias_energy = bias_energy;
  cvm::real const old_force_k = force_k;

  change_configuration(conf);
  update();

  cvm::real const result = (bias_energy - old_bias_energy);

  bias_energy = old_bias_energy;
  force_k = old_force_k;

  return result;
}


cvm::real colvarbias_restraint_linear::restraint_potential(size_t i) const
{
  return force_k / variables(i)->width * (variables(i)->value() -
                                          colvar_centers[i]).sum();
}


colvarvalue const colvarbias_restraint_linear::restraint_force(size_t i) const
{
  colvarvalue dummy(variables(i)->value());
  dummy.set_ones();
  return -1.0 * force_k / variables(i)->width * dummy;
}


cvm::real colvarbias_restraint_linear::d_restraint_potential_dk(size_t i) const
{
  return 1.0 / variables(i)->width * (variables(i)->value() -
                                      colvar_centers[i]).sum();
}


std::string const colvarbias_restraint_linear::get_state_params() const
{
  return colvarbias_restraint::get_state_params() +
    colvarbias_restraint_moving::get_state_params() +
    colvarbias_restraint_centers_moving::get_state_params() +
    colvarbias_restraint_k_moving::get_state_params();
}


int colvarbias_restraint_linear::set_state_params(std::string const &conf)
{
  int error_code = COLVARS_OK;
  error_code |= colvarbias_restraint::set_state_params(conf);
  error_code |= colvarbias_restraint_moving::set_state_params(conf);
  error_code |= colvarbias_restraint_centers_moving::set_state_params(conf);
  error_code |= colvarbias_restraint_k_moving::set_state_params(conf);
  return error_code;
}


std::ostream & colvarbias_restraint_linear::write_traj_label(std::ostream &os)
{
  colvarbias_restraint::write_traj_label(os);
  colvarbias_restraint_centers_moving::write_traj_label(os);
  colvarbias_restraint_k_moving::write_traj_label(os);
  return os;
}


std::ostream & colvarbias_restraint_linear::write_traj(std::ostream &os)
{
  colvarbias_restraint::write_traj(os);
  colvarbias_restraint_centers_moving::write_traj(os);
  colvarbias_restraint_k_moving::write_traj(os);
  return os;
}



colvarbias_restraint_histogram::colvarbias_restraint_histogram(colvarmodule *cvmodule_in, char const *key)
  : colvarbias(cvmodule_in, key)
{
  lower_boundary = 0.0;
  upper_boundary = 0.0;
  width = 0.0;
  gaussian_width = 0.0;
}


int colvarbias_restraint_histogram::init(std::string const &conf)
{
  int error_code = COLVARS_OK;

  int err = colvarbias::init(conf);
  if (err != COLVARS_OK) {
    return err;
  }
  enable(f_cvb_apply_force);

  cvmodule->cite_feature("histogramRestraint colvar bias implementation");

  get_keyval(conf, "lowerBoundary", lower_boundary, lower_boundary);
  get_keyval(conf, "upperBoundary", upper_boundary, upper_boundary);
  get_keyval(conf, "width", width, width);

  if (width <= 0.0) {
    error_code |= cvmodule->error("Error: \"width\" must be positive.\n",
                             COLVARS_INPUT_ERROR);
  }

  get_keyval(conf, "gaussianWidth", gaussian_width, 2.0 * width, colvarparse::parse_silent);
  get_keyval(conf, "gaussianSigma", gaussian_width, 2.0 * width);

  if (lower_boundary >= upper_boundary) {
    error_code |= cvmodule->error("Error: the upper boundary, "+
                             cvm::to_str(upper_boundary)+
                             ", is not higher than the lower boundary, "+
                             cvm::to_str(lower_boundary)+".\n",
                             COLVARS_INPUT_ERROR);
  }

  cvm::real const nbins = (upper_boundary - lower_boundary) / width;
  int const nbins_round = (int)(nbins);

  if (cvm::fabs(nbins - cvm::real(nbins_round)) > 1.0E-10) {
    cvmodule->log("Warning: grid interval ("+
             cvm::to_str(lower_boundary, cvmodule->cv_width, cvmodule->cv_prec)+" - "+
             cvm::to_str(upper_boundary, cvmodule->cv_width, cvmodule->cv_prec)+
             ") is not commensurate to its bin width ("+
             cvm::to_str(width, cvmodule->cv_width, cvmodule->cv_prec)+").\n");
  }

  p.resize(nbins_round);
  ref_p.resize(nbins_round);
  p_diff.resize(nbins_round);

  bool const inline_ref_p =
    get_keyval(conf, "refHistogram", ref_p.data_array(), ref_p.data_array());
  std::string ref_p_file;
  get_keyval(conf, "refHistogramFile", ref_p_file, std::string(""));
  if (ref_p_file.size()) {
    if (inline_ref_p) {
      error_code |= cvmodule->error("Error: cannot specify both refHistogram and refHistogramFile at the same time.\n",
                               COLVARS_INPUT_ERROR);
    } else {

      std::istream &is =
        cvmodule->proxy->input_stream(ref_p_file,
                                         "reference histogram file");

      std::string data_s = "";
      std::string line;
      while (getline_nocomments(is, line)) {
        data_s.append(line+"\n");
      }
      if (data_s.size() == 0) {
        error_code |= cvmodule->error("Error: file \""+ref_p_file+
                                 "\" empty or unreadable.\n",
                                 COLVARS_FILE_ERROR);
      }
      error_code |= cvmodule->proxy->close_input_stream(ref_p_file);

      cvm::vector1d<cvm::real> data;
      if (data.from_simple_string(data_s) != 0) {
        error_code |= cvmodule->error("Error: could not read histogram from file \""+
                                 ref_p_file+"\".\n");
      }
      if (data.size() == 2*ref_p.size()) {
        // file contains both x and p(x)
        size_t i;
        for (i = 0; i < ref_p.size(); i++) {
          ref_p[i] = data[2*i+1];
        }
      } else if (data.size() == ref_p.size()) {
        ref_p = data;
      } else {
        error_code |= cvmodule->error("Error: file \""+ref_p_file+
                                 "\" contains a histogram of different length.\n",
                                 COLVARS_INPUT_ERROR);
      }
    }
  }

  cvm::real const ref_integral = ref_p.sum() * width;
  if (cvm::fabs(ref_integral - 1.0) > 1.0e-03) {
    cvmodule->log("Reference distribution not normalized, normalizing to unity.\n");
    ref_p /= ref_integral;
  }

  get_keyval(conf, "writeHistogram", b_write_histogram, false);
  get_keyval(conf, "forceConstant", force_k, 1.0);

  return error_code;
}


colvarbias_restraint_histogram::~colvarbias_restraint_histogram()
{
  p.clear();
  ref_p.clear();
  p_diff.clear();
}


int colvarbias_restraint_histogram::update()
{
  if (cvm::debug())
    cvmodule->log("Updating the histogram restraint bias \""+this->name+"\".\n");

  size_t vector_size = 0;
  size_t icv;
  for (icv = 0; icv < num_variables(); icv++) {
    vector_size += variables(icv)->value().size();
  }

  cvm::real const norm = 1.0/(cvm::sqrt(2.0*PI)*gaussian_width*vector_size);

  // calculate the histogram
  p.reset();
  for (icv = 0; icv < num_variables(); icv++) {
    colvarvalue const &cv = variables(icv)->value();
    if (cv.type() == colvarvalue::type_scalar) {
      cvm::real const cv_value = cv.real_value;
      size_t igrid;
      for (igrid = 0; igrid < p.size(); igrid++) {
        cvm::real const x_grid = (lower_boundary + (igrid+0.5)*width);
        p[igrid] += norm * cvm::exp(-1.0 * (x_grid - cv_value) * (x_grid - cv_value) /
                                    (2.0 * gaussian_width * gaussian_width));
      }
    } else if (cv.type() == colvarvalue::type_vector) {
      size_t idim;
      for (idim = 0; idim < cv.vector1d_value.size(); idim++) {
        cvm::real const cv_value = cv.vector1d_value[idim];
        size_t igrid;
        for (igrid = 0; igrid < p.size(); igrid++) {
          cvm::real const x_grid = (lower_boundary + (igrid+0.5)*width);
          p[igrid] += norm * cvm::exp(-1.0 * (x_grid - cv_value) * (x_grid - cv_value) /
                                      (2.0 * gaussian_width * gaussian_width));
        }
      }
    } else {
      cvmodule->error("Error: unsupported type for variable "+variables(icv)->name+".\n",
                 COLVARS_NOT_IMPLEMENTED);
      return COLVARS_NOT_IMPLEMENTED;
    }
  }

  cvm::real const force_k_cv = force_k * vector_size;

  // calculate the difference between current and reference
  p_diff = p - ref_p;
  bias_energy = 0.5 * force_k_cv * p_diff * p_diff;

  // calculate the forces
  for (icv = 0; icv < num_variables(); icv++) {
    colvarvalue const &cv = variables(icv)->value();
    colvarvalue &cv_force = colvar_forces[icv];
    cv_force.type(cv);
    cv_force.reset();

    if (cv.type() == colvarvalue::type_scalar) {
      cvm::real const cv_value = cv.real_value;
      cvm::real &force = cv_force.real_value;
      size_t igrid;
      for (igrid = 0; igrid < p.size(); igrid++) {
        cvm::real const x_grid = (lower_boundary + (igrid+0.5)*width);
        force += force_k_cv * p_diff[igrid] *
          norm * cvm::exp(-1.0 * (x_grid - cv_value) * (x_grid - cv_value) /
                          (2.0 * gaussian_width * gaussian_width)) *
          (-1.0 * (x_grid - cv_value) / (gaussian_width * gaussian_width));
      }
    } else if (cv.type() == colvarvalue::type_vector) {
      size_t idim;
      for (idim = 0; idim < cv.vector1d_value.size(); idim++) {
        cvm::real const cv_value = cv.vector1d_value[idim];
        cvm::real &force = cv_force.vector1d_value[idim];
        size_t igrid;
        for (igrid = 0; igrid < p.size(); igrid++) {
          cvm::real const x_grid = (lower_boundary + (igrid+0.5)*width);
          force += force_k_cv * p_diff[igrid] *
            norm * cvm::exp(-1.0 * (x_grid - cv_value) * (x_grid - cv_value) /
                            (2.0 * gaussian_width * gaussian_width)) *
            (-1.0 * (x_grid - cv_value) / (gaussian_width * gaussian_width));
        }
      }
    } else {
      // TODO
    }
  }

  return COLVARS_OK;
}


int colvarbias_restraint_histogram::write_output_files()
{
  if (b_write_histogram) {
    colvarproxy *proxy = cvmodule->proxy;
    std::string file_name(cvmodule->output_prefix()+"."+this->name+".hist.dat");
    std::ostream &os = proxy->output_stream(file_name,
                                            "histogram output file");
    os << "# " << cvmodule->wrap_string(variables(0)->name, cvmodule->cv_width)
       << "  " << "p(" << cvmodule->wrap_string(variables(0)->name, cvmodule->cv_width-3)
       << ")\n";

    os.setf(std::ios::fixed, std::ios::floatfield);

    size_t igrid;
    for (igrid = 0; igrid < p.size(); igrid++) {
      cvm::real const x_grid = (lower_boundary + (igrid+1)*width);
      os << "  "
         << std::setprecision(cvmodule->cv_prec)
         << std::setw(cvmodule->cv_width)
         << x_grid
         << "  "
         << std::setprecision(cvmodule->cv_prec)
         << std::setw(cvmodule->cv_width)
         << p[igrid] << "\n";
    }
    proxy->close_output_stream(file_name);
  }
  return COLVARS_OK;
}


std::ostream & colvarbias_restraint_histogram::write_traj_label(std::ostream &os)
{
  os << " ";
  if (b_output_energy) {
    os << " E_"
       << cvmodule->wrap_string(this->name, cvmodule->en_width-2);
  }
  return os;
}


std::ostream & colvarbias_restraint_histogram::write_traj(std::ostream &os)
{
  os << " ";
  if (b_output_energy) {
    os << " "
       << std::setprecision(cvmodule->en_prec) << std::setw(cvmodule->en_width)
       << bias_energy;
  }
  return os;
}
