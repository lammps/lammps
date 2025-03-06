// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.


#include "colvarbias_abmd.h"
#include "colvarproxy.h"
#include <iomanip>


colvarbias_abmd::colvarbias_abmd(char const *key)
  : colvarbias(key),
    colvarbias_ti(key)
{
}


int colvarbias_abmd::init(std::string const &conf)
{
  cvm::main()->cite_feature("ABMD bias");

  int err = colvarbias::init(conf);
  err |= colvarbias_ti::init(conf);
  if (err != COLVARS_OK) return err;

  enable(f_cvb_apply_force);

  if (num_variables() != 1) {
    return cvm::error("ABMD requires exactly one collective variable.\n", COLVARS_INPUT_ERROR);
  }

  if ( ! (variables(0))->is_enabled(f_cv_scalar) ) {
    return cvm::error("ABMD colvar must be scalar.\n", COLVARS_INPUT_ERROR);
  }

  get_keyval(conf, "forceConstant", k);
  get_keyval(conf, "decreasing", decreasing, decreasing);
  get_keyval(conf, "stoppingValue", stopping_val);

  return COLVARS_OK;
}


int colvarbias_abmd::update()
{
  if (!cvm::main()->proxy->simulation_running()) {
    return COLVARS_OK;
  }

  colvar const *cv = variables(0);
  cvm::real const val = cv->value().real_value;

  if (!ref_initialized) {
    ref_val = val;
    ref_initialized = true;
  }

  // Compute sign factor to unify increasing and decreasing cases below
  // less conditionals, more arithmetic
  cvm::real const sign = decreasing ? -1. : 1.;
  cvm::real const diff = (val - ref_val) * sign;

  if ( diff > 0. ) {
    colvar_forces[0] = 0.;
    bias_energy = 0.;
    if ( (ref_val-stopping_val) * sign <= 0. ) ref_val = val;
  } else {
    colvar_forces[0] = - sign * k * diff;
    bias_energy = 0.5 * k * diff * diff;;
  }
  return COLVARS_OK;
}


std::string const colvarbias_abmd::get_state_params() const
{
  std::ostringstream os;
  os.setf(std::ios::scientific, std::ios::floatfield);

  os << "    refValue "
      << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
      << ref_val << "\n";
  os << "    stoppingValue " << stopping_val << "\n";
  os << "    forceConstant " << k << "\n";
  os << "    decreasing " << (decreasing ? "on" : "off") << "\n";

  return (colvarbias::get_state_params() + os.str());
}


int colvarbias_abmd::set_state_params(std::string const &conf)
{
  int error_code = colvarbias::set_state_params(conf);

  if (error_code != COLVARS_OK) {
    return error_code;
  }

  get_keyval(conf, "refValue", ref_val, ref_val,
    colvarparse::parse_restart | colvarparse::parse_required);
  ref_initialized = true;

  get_keyval(conf, "forceConstant", k, k,
    colvarparse::parse_restart | colvarparse::parse_required);
  get_keyval(conf, "decreasing", decreasing, decreasing,
    colvarparse::parse_restart | colvarparse::parse_required);
  get_keyval(conf, "stoppingValue", stopping_val, stopping_val,
    colvarparse::parse_restart | colvarparse::parse_required);

  return COLVARS_OK;
}


std::ostream & colvarbias_abmd::write_traj_label(std::ostream &os)
{
  size_t const this_cv_width = (variables(0)->value()).output_width(cvm::cv_width);
  os << " ref_"
      << cvm::wrap_string(variables(0)->name, this_cv_width-4);

  return os;
}


std::ostream & colvarbias_abmd::write_traj(std::ostream &os)
{
  os << " "
      << std::setprecision(cvm::en_prec) << std::setw(cvm::en_width)
      << ref_val;

  return os;
}
