// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <iostream>

#include "colvarmodule.h"
#include "colvar.h"
#include "colvarbias_abf.h"
#include "colvars_memstream.h"

colvarbias_abf::colvarbias_abf(char const *key)
  : colvarbias(key),
    b_UI_estimator(false),
    b_CZAR_estimator(false),
    pabf_freq(0),
    system_force(NULL)
{
}


int colvarbias_abf::init(std::string const &conf)
{
  colvarproxy *proxy = cvm::main()->proxy;

  int err = colvarbias::init(conf);
  if (err != COLVARS_OK) {
    return err;
  }
  cvm::main()->cite_feature("ABF colvar bias implementation");

  enable(f_cvb_scalar_variables);
  enable(f_cvb_calc_pmf);

  if ((proxy->target_temperature() == 0.0) && proxy->simulation_running()) {
    cvm::log("WARNING: ABF should not be run without a thermostat or at 0 Kelvin!\n");
  }

  // ************* parsing general ABF options ***********************

  get_keyval_feature(this, conf, "applyBias", f_cvb_apply_force, true);
  if (!is_enabled(f_cvb_apply_force)){
    cvm::log("WARNING: ABF biases will *not* be applied!\n");
  }

  get_keyval(conf, "hideJacobian", hide_Jacobian, false);
  if (hide_Jacobian) {
    cvm::log("Jacobian (geometric) forces will be handled internally.\n");
  } else {
    cvm::log("Jacobian (geometric) forces will be included in reported free energy gradients.\n");
  }

  full_samples = 200;
  get_keyval(conf, "fullSamples", full_samples, full_samples);
  get_keyval(conf, "minSamples", min_samples, full_samples / 2);

  if (full_samples <= 1 ) {
    full_samples = 1;
    min_samples = 0;
  }
  if (min_samples >= full_samples) {
    return cvm::error("Error: minSamples must be lower than fullSamples\n");
  }

  get_keyval(conf, "inputPrefix",  input_prefix, std::vector<std::string>());

  history_last_step = -1;
  get_keyval(conf, "historyFreq", history_freq, 0);
  if (history_freq != 0) {
    if (output_freq == 0) {
      cvm::error("Error: historyFreq must be a multiple of outputFreq.\n",
                 COLVARS_INPUT_ERROR);
    } else {
      if ((history_freq % output_freq) != 0) {
        cvm::error("Error: historyFreq must be a multiple of outputFreq.\n",
                   COLVARS_INPUT_ERROR);
      }
    }
  }

  // shared ABF
  get_keyval(conf, "shared", shared_on, false);
  if (shared_on) {
    cvm::main()->cite_feature("Multiple-walker ABF implementation");
    if ((proxy->replica_enabled() != COLVARS_OK) ||
        (proxy->num_replicas() <= 1)) {
      return cvm::error("Error: shared ABF requires more than one replica.",
                        COLVARS_INPUT_ERROR);
    }
    cvm::log("shared ABF will be applied among "+
             cvm::to_str(proxy->num_replicas()) + " replicas.\n");

    // If shared_freq is not set, we default to output_freq
    get_keyval(conf, "sharedFreq", shared_freq, output_freq);
    if ( shared_freq && output_freq % shared_freq ) {
      return cvm::error("Error: outputFreq must be a multiple of sharedFreq.\n");
    }

    // Allocate these at init time if possible
    local_samples.reset(new colvar_grid_count(colvars));
    local_gradients.reset(new colvar_grid_gradient(colvars, local_samples));
    local_pmf.reset(new integrate_potential(colvars, local_gradients));
  }

  // ************* checking the associated colvars *******************

  if (num_variables() == 0) {
    return cvm::error("Error: no collective variables specified for the ABF bias.\n");
  }

  size_t i;
  for (i = 0; i < num_variables(); i++) {

    if (colvars[i]->value().type() != colvarvalue::type_scalar) {
      cvm::error("Error: ABF bias can only use scalar-type variables.\n");
    }
    colvars[i]->enable(f_cv_grid); // Could be a child dependency of a f_cvb_use_grids feature
    if (hide_Jacobian) {
      colvars[i]->enable(f_cv_hide_Jacobian);
    }

    // If any colvar is extended-system, we need to collect the extended
    // system gradient
    if (colvars[i]->is_enabled(f_cv_extended_Lagrangian))
      enable(f_cvb_extended);

    // Cannot mix and match coarse time steps with ABF because it gives
    // wrong total force averages - total force needs to be averaged over
    // every time step
    if (colvars[i]->get_time_step_factor() != time_step_factor) {
      cvm::error("Error: " + colvars[i]->description + " has a value of timeStepFactor ("
        + cvm::to_str(colvars[i]->get_time_step_factor()) + ") different from that of "
        + description + " (" + cvm::to_str(time_step_factor) + ").\n");
      return COLVARS_ERROR;
    }

    // Here we could check for orthogonality of the Cartesian coordinates
    // and make it just a warning if some parameter is set?
  }

  get_keyval(conf, "updateBias",  update_bias, true);
  if (update_bias) {
    enable(f_cvb_history_dependent);
    enable(f_cvb_get_total_force);
  } else {
    cvm::log("WARNING: ABF biases will *not* be updated!\n");
  }

  if (is_enabled(f_cvb_extended)) {
    cvm::main()->cite_feature("eABF implementation");
  } else {
    cvm::main()->cite_feature("Internal-forces free energy estimator");
  }

  if (get_keyval(conf, "maxForce", max_force)) {
    if (max_force.size() != num_variables()) {
      cvm::error("Error: Number of parameters to maxForce does not match number of colvars.");
    }
    for (i = 0; i < num_variables(); i++) {
      if (max_force[i] < 0.0) {
        cvm::error("Error: maxForce should be non-negative.");
        return COLVARS_ERROR;
      }
    }
    cap_force = true;
  } else {
    cap_force = false;
  }

  bin.assign(num_variables(), 0);
  force_bin.assign(num_variables(), 0);
  system_force = new cvm::real [num_variables()];

  // Construct empty grids based on the colvars
  if (cvm::debug()) {
    cvm::log("Allocating count and free energy gradient grids.\n");
  }

  samples.reset(new colvar_grid_count(colvars));
  gradients.reset(new colvar_grid_gradient(colvars, samples));

  gradients->full_samples = full_samples;
  gradients->min_samples = min_samples;

  // Data for eABF z-based estimator
  if (is_enabled(f_cvb_extended)) {
    get_keyval(conf, "CZARestimator", b_CZAR_estimator, true);
    if ( b_CZAR_estimator ) {
      cvm::main()->cite_feature("CZAR eABF estimator");
    }
    // CZAR output files for stratified eABF
    get_keyval(conf, "writeCZARwindowFile", b_czar_window_file, false,
               colvarparse::parse_silent);

    z_bin.assign(num_variables(), 0);
    z_samples.reset(new colvar_grid_count(colvars));
    z_samples->request_actual_value();
    z_gradients.reset(new colvar_grid_gradient(colvars, z_samples));
    z_gradients->request_actual_value();
    czar_gradients.reset(new colvar_grid_gradient(colvars));
  }

  get_keyval(conf, "integrate", b_integrate, num_variables() <= 3); // Integrate for output if d<=3
  if (b_integrate) {
    // For now, we integrate on-the-fly iff the grid is < 3D
    if ( num_variables() > 3 ) {
      cvm::error("Error: cannot integrate free energy in dimension > 3.\n");
      return COLVARS_ERROR;
    }
    pmf.reset(new integrate_potential(colvars, gradients));
    if (b_CZAR_estimator) {
      czar_pmf.reset(new integrate_potential(colvars, czar_gradients));
    }
    // Parameters for integrating initial (and final) gradient data
    get_keyval(conf, "integrateMaxIterations", integrate_iterations, 10000, colvarparse::parse_silent);
    get_keyval(conf, "integrateTol", integrate_tol, 1e-6, colvarparse::parse_silent);
    // Projected ABF, updating the integrated PMF on the fly
    get_keyval(conf, "pABFintegrateFreq", pabf_freq, 0, colvarparse::parse_silent);
    get_keyval(conf, "pABFintegrateMaxIterations", pabf_integrate_iterations, 100, colvarparse::parse_silent);
    get_keyval(conf, "pABFintegrateTol", pabf_integrate_tol, 1e-4, colvarparse::parse_silent);
  }

  if (b_CZAR_estimator && shared_on && cvm::main()->proxy->replica_index() == 0) {
    // The pointers below are used for outputting CZAR data
    // Allocate grids for collected global data, on replica 0 only
    global_z_samples.reset(new colvar_grid_count(colvars));
    global_z_gradients.reset(new colvar_grid_gradient(colvars, global_z_samples));
    global_czar_gradients.reset(new colvar_grid_gradient(colvars));
    global_czar_pmf.reset(new integrate_potential(colvars, global_czar_gradients));
  } else {
    // otherwise they are just aliases for the local CZAR grids
    global_z_samples = z_samples;
    global_z_gradients = z_gradients;
    global_czar_gradients = czar_gradients;
    global_czar_pmf = czar_pmf;
  }

  // For shared ABF, we store a second set of grids.
  // This used to be only if "shared" was defined,
  // but now we allow calling share externally (e.g. from Tcl).
  if (b_CZAR_estimator) {
    z_samples_in.reset(new colvar_grid_count(colvars));
    z_gradients_in.reset(new colvar_grid_gradient(colvars, z_samples_in));
  }
  last_samples.reset(new colvar_grid_count(colvars));
  last_gradients.reset(new colvar_grid_gradient(colvars, last_samples));
  // Any data collected after now is new for shared ABF purposes
  shared_last_step = cvm::step_absolute();

  // Read any custom input ABF data
  if ( input_prefix.size() > 0 ) {
    read_gradients_samples();
    // Update divergence to account for input data
    pmf->set_div();
  }

  // if extendedLangrangian is on, then call UI estimator
  if (is_enabled(f_cvb_extended)) {
    get_keyval(conf, "UIestimator", b_UI_estimator, false);

    if (b_UI_estimator) {
      if (shared_on) {
        cvm::error("Error: UI estimator is not available for multiple-walker (shared) ABF.\n");
        b_UI_estimator = false;
        return COLVARS_ERROR;
      }
      cvm::main()->cite_feature("Umbrella-integration eABF estimator");
      std::vector<double> UI_lowerboundary;
      std::vector<double> UI_upperboundary;
      std::vector<double> UI_width;
      std::vector<double> UI_krestr;

      bool UI_restart = (input_prefix.size() > 0);

      for (i = 0; i < num_variables(); i++) {
        UI_lowerboundary.push_back(colvars[i]->lower_boundary);
        UI_upperboundary.push_back(colvars[i]->upper_boundary);
        UI_width.push_back(colvars[i]->width);
        UI_krestr.push_back(colvars[i]->force_constant());
      }
      eabf_UI = UIestimator::UIestimator(UI_lowerboundary,
                                         UI_upperboundary,
                                         UI_width,
                                         UI_krestr,                // force constant in eABF
                                         output_prefix,              // the prefix of output files
                                         cvm::restart_out_freq,
                                         UI_restart,                    // whether restart from a .count and a .grad file
                                         input_prefix,   // the prefixes of input files
                                         proxy->target_temperature());
    }
  }

  cvm::log("Finished ABF setup.\n");
  return COLVARS_OK;
}

/// Destructor
colvarbias_abf::~colvarbias_abf()
{
  if (system_force) delete[] system_force;
}


/// Update the FE gradient, compute and apply biasing force

int colvarbias_abf::update()
{
  if (cvm::debug()) cvm::log("Updating ABF bias " + this->name);

  size_t i;
  for (i = 0; i < num_variables(); i++) {
    bin[i] = samples->current_bin_scalar(i);
  }


  // ***********************************************************
  // ******  ABF Part I: update the FE gradient estimate  ******
  // ***********************************************************


  if (cvm::proxy->total_forces_same_step()) {
    // e.g. in LAMMPS, total forces are current
    force_bin = bin;
  }

  if (can_accumulate_data() && is_enabled(f_cvb_history_dependent)) {

    if (cvm::step_relative() > 0 || cvm::proxy->total_forces_same_step()) {
      if (samples->index_ok(force_bin)) {
        // Only if requested and within bounds of the grid...

        // get total forces (lagging by 1 timestep) from colvars
        // and subtract previous ABF force if necessary
        update_system_force();

        gradients->acc_force(force_bin, system_force);
        if ( b_integrate ) {
          pmf->update_div_neighbors(force_bin);
        }
      }

      if ( z_gradients ) {
        for (i = 0; i < num_variables(); i++) {
          z_bin[i] = z_samples->current_bin_scalar(i);
        }
        if ( z_samples->index_ok(z_bin) ) {
          // If we are outside the range of z, the force has not been obtained above
          // the function is just an accessor, so cheap to call again anyway
          update_system_force();
          z_gradients->acc_force(z_bin, system_force);
        }
      }

      if ( pabf_freq && cvm::step_relative() % pabf_freq == 0 ) {
        cvm::real err;
        int iter = pmf->integrate(integrate_iterations, integrate_tol, err);
        if ( iter == integrate_iterations ) {
          cvm::log("Warning: PMF integration did not converge to " + cvm::to_str(integrate_tol)
            + " in " + cvm::to_str(integrate_iterations)
            + " steps. Residual error: " +  cvm::to_str(err));
        }
      }
    }
  }

  if (!(cvm::proxy->total_forces_same_step())) {
    // e.g. in NAMD, total forces will be available for next timestep
    // hence we store the current colvar bin
    force_bin = bin;
  }

  // Share data after force sample is collected for this time step
  // shared_on can be true with shared_freq 0 if we are sharing via script
  if (shared_on && shared_freq &&
      cvm::step_absolute() > shared_last_step &&  // time has passed since the last sharing timestep
                                                  // (avoid re-sharing at last and first ts of successive run statements)
      cvm::step_absolute() % shared_freq == 0) {
    // Share gradients and samples for shared ABF.
    replica_share();
  }

  // ******************************************************************
  // ******  ABF Part II: calculate and apply the biasing force  ******
  // ******************************************************************


  // Reset biasing forces from previous timestep
  for (i = 0; i < num_variables(); i++) {
    colvar_forces[i].reset();
  }

  // Compute and apply the new bias, if applicable
  if (is_enabled(f_cvb_apply_force) && samples->index_ok(bin)) {

    std::vector<cvm::real>  force(num_variables());
    calc_biasing_force(force);

    for (size_t i = 0; i < num_variables(); i++) {
      colvar_forces[i].real_value = force[i];
    }
  }


  // *********************************
  // ******  End of ABF proper  ******
  // *********************************


  // update the output prefix; TODO: move later to setup_output() function
  if (cvm::main()->num_biases_feature(colvardeps::f_cvb_calc_pmf) == 1) {
    // This is the only bias computing PMFs
    output_prefix = cvm::output_prefix();
  } else {
    output_prefix = cvm::output_prefix() + "." + this->name;
  }


  // update UI estimator every step
  if (b_UI_estimator)
  {
    std::vector<double> x(num_variables(),0);
    std::vector<double> y(num_variables(),0);
    for (i = 0; i < num_variables(); i++)
    {
      x[i] = colvars[i]->actual_value();
      y[i] = colvars[i]->value();
    }
    eabf_UI.update_output_filename(output_prefix);
    eabf_UI.update(cvm::step_absolute(), x, y);
  }

  /// Compute the bias energy
  int error_code = calc_energy(NULL);

  return error_code;
}


  // ****************************************
  // ******  Helper functions for ABF  ******
  // ****************************************


int colvarbias_abf::update_system_force()
{
  size_t i;
  // System force from atomic forces (or extended Lagrangian if applicable)

  for (i = 0; i < num_variables(); i++) {
    if (colvars[i]->is_enabled(f_cv_subtract_applied_force)) {
      // this colvar is already subtracting the ABF force
      system_force[i] = colvars[i]->total_force().real_value;
    } else {
      system_force[i] = colvars[i]->total_force().real_value
        - colvar_forces[i].real_value;
    }
  }
  return COLVARS_OK;
}


cvm::real colvarbias_abf::smoothing_factor(cvm::real weight)
{
  cvm::real fact = 1.0;
  if ( weight < full_samples ) {
    if ( weight < min_samples) {
      fact = 0.0;
    } else {
      fact = (weight - min_samples) / cvm::real(full_samples - min_samples);
    }
  }
  return fact;
}


int colvarbias_abf::calc_biasing_force(std::vector<cvm::real> &force)
{
  size_t i;

  // Pick between different types of biasing force
  if ( pabf_freq ) {
    // In projected ABF, the force is the PMF gradient estimate
    pmf->vector_gradient_finite_diff(bin, force);
    // Calculate ramp factor that ensures smooth introduction of the force
    const cvm::real count = samples->value(bin);
    const cvm::real fact = smoothing_factor(count);
    for (i = 0; i < num_variables(); i++) {
      force[i] *= fact;
    }
  } else {
    // Normal ABF or eABF: use accumulated gradient average
    gradients->vector_value_smoothed(bin, &force[0], true);
    if ( (num_variables() == 1) && gradients->periodic[0] ) {
      // Enforce a zero-mean bias on periodic, 1D coordinates
      // in other words: boundary condition is that the biasing potential is periodic
      // Only plain ABF needs this
      force[0] = force[0] - gradients->average();
    }
  }

  if (cap_force) {
    for (i = 0; i < num_variables(); i++) {
      if ( force[i] * force[i] > max_force[i] * max_force[i] ) {
        force[i] = (force[i] > 0 ? max_force[i] : -1.0 * max_force[i]);
      }
    }
  }

  return COLVARS_OK;
}


  // ************************************
  // ******  Shared ABF functions  ******
  // ************************************



int colvarbias_abf::replica_share() {

  colvarproxy *proxy = cvm::main()->proxy;

  if (proxy->replica_enabled() != COLVARS_OK) {
    cvm::error("Error: shared ABF: No replicas.\n");
    return COLVARS_ERROR;
  }
  // We must have stored the last_gradients and last_samples.
  if (shared_last_step < 0 ) {
    cvm::error("Error: shared ABF: Tried to apply shared ABF before any sampling had occurred.\n");
    return COLVARS_ERROR;
  }
  shared_on = true; // If called by a script, inform the rest of the code that we're sharing, eg. CZAR

  // Share gradients for shared ABF.
  cvm::log("shared ABF: Sharing gradient and samples among replicas at step "+cvm::to_str(cvm::step_absolute()) );

  if (!local_samples) {
    // We arrive here if sharing has just been enabled by a script
    // in which case local arrays have not been initialized yet
    local_samples.reset(new colvar_grid_count(colvars));
    local_gradients.reset(new colvar_grid_gradient(colvars, local_samples));
    local_pmf.reset(new integrate_potential(colvars, local_gradients));
  }
  // Calculate the delta gradient and count for the local replica
  last_gradients->delta_grid(*gradients);
  // Add the delta gradient and count to the accumulated local data
  local_gradients->add_grid(*last_gradients);

  last_samples->delta_grid(*samples);
  local_samples->add_grid(*last_samples);


  // Count of data items.
  size_t samples_n = samples->raw_data_num();
  size_t gradients_n = gradients->raw_data_num();

  size_t samp_start = gradients_n * sizeof(cvm::real);
  int msg_total = samples_n * sizeof(size_t) + samp_start;
  char* msg_data = new char[msg_total];

  if (cvm::main()->proxy->replica_index() == 0) {
    int p;
    // Replica 0 collects the delta gradient and count from the others.
    for (p = 1; p < proxy->num_replicas(); p++) {
      // Receive the deltas.
      if (proxy->replica_comm_recv(msg_data, msg_total, p) != msg_total) {
        cvm::error("Error getting shared ABF data from replica.");
        return COLVARS_ERROR;
      }

      // Map the deltas from the others into the grids.
      // Re-use last_gradients as temp array, erasing its contents each time
      last_gradients->raw_data_in((cvm::real*)(&msg_data[0]));
      // Combine the delta gradient and count of the other replicas
      // with Replica 0's current state (including its delta).
      gradients->add_grid(*last_gradients);

      last_samples->raw_data_in((size_t*)(&msg_data[samp_start]));
      samples->add_grid(*last_samples);
    }

    // Now we must send the combined gradient to the other replicas.
    gradients->raw_data_out((cvm::real*)(&msg_data[0]));
    samples->raw_data_out((size_t*)(&msg_data[samp_start]));

    for (p = 1; p < proxy->num_replicas(); p++) {
      if (proxy->replica_comm_send(msg_data, msg_total, p) != msg_total) {
        cvm::error("Error sending shared ABF data to replica.");
        return COLVARS_ERROR;
      }
    }

  } else {
    // All other replicas send their delta gradient and count.
    // Cast the raw char data to the gradient and samples.
    last_gradients->raw_data_out((cvm::real*)(&msg_data[0]));
    last_samples->raw_data_out((size_t*)(&msg_data[samp_start]));

    if (proxy->replica_comm_send(msg_data, msg_total, 0) != msg_total) {
      cvm::error("Error sending shared ABF data to replica.");
      return COLVARS_ERROR;
    }
    // We now receive the combined gradient from Replica 0.
    if (proxy->replica_comm_recv(msg_data, msg_total, 0) != msg_total) {
      cvm::error("Error getting shared ABF data from replica 0.");
      return COLVARS_ERROR;
    }
    // We sync to the combined gradient computed by Replica 0.
    gradients->raw_data_in((cvm::real*)(&msg_data[0]));
    samples->raw_data_in((size_t*)(&msg_data[samp_start]));
  }

  // Without a barrier it's possible that one replica starts
  // share 2 when other replicas haven't finished share 1.
  proxy->replica_comm_barrier();
  // Done syncing the replicas.
  delete[] msg_data;

  // Copy the current gradient and count values into last.
  last_gradients->copy_grid(*gradients);
  last_samples->copy_grid(*samples);
  shared_last_step = cvm::step_absolute();

  cvm::log("RMSD btw. local and global ABF gradients: " + cvm::to_str(gradients->grid_rmsd(*local_gradients)));

  if (b_integrate) {
    cvm::real err;

    // Update whole divergence field to account for newly shared gradients
    pmf->set_div();
    pmf->integrate(integrate_iterations, integrate_tol, err);
    pmf->set_zero_minimum();
    local_pmf->set_div();
    local_pmf->integrate(integrate_iterations, integrate_tol, err);
    local_pmf->set_zero_minimum();
    cvm::log("RMSD btw. local and global ABF FES: " + cvm::to_str(pmf->grid_rmsd(*local_pmf)));
  }
  return COLVARS_OK;
}


int colvarbias_abf::replica_share_CZAR() {
  colvarproxy *proxy = cvm::main()->proxy;

  cvm::log("shared eABF: Gathering CZAR gradient and samples from replicas at step "+cvm::to_str(cvm::step_absolute()) );

  // Count of data items.
  size_t samples_n = z_samples->raw_data_num();
  size_t gradients_n = z_gradients->raw_data_num();

  size_t samp_start = gradients_n*sizeof(cvm::real);
  int msg_total = samples_n*sizeof(size_t) + samp_start;
  char* msg_data = new char[msg_total];

  if (cvm::main()->proxy->replica_index() == 0) {
     if (!global_z_samples) {
      // We arrive here if sharing has just been enabled by a script
      // Allocate grids for collective data, on replica 0 only
      // overriding CZAR grids that are equal to local ones by default
      global_z_samples.reset(new colvar_grid_count(colvars));
      global_z_gradients.reset(new colvar_grid_gradient(colvars, global_z_samples));
      global_czar_gradients.reset(new colvar_grid_gradient(colvars));
      global_czar_pmf.reset(new integrate_potential(colvars, global_czar_gradients));
    }

    // Start with data from replica 0
    global_z_gradients->copy_grid(*z_gradients);
    global_z_samples->copy_grid(*z_samples);

    int p;
    // Replica 0 collects the gradient and count from the others.
    for (p = 1; p < proxy->num_replicas(); p++) {
      if (proxy->replica_comm_recv(msg_data, msg_total, p) != msg_total) {
        cvm::error("Error getting shared ABF data from replica.");
        return COLVARS_ERROR;
      }

      // Map the deltas from the others into the grids.
      // Re-use z_gradients_in, erasing its contents each time
      z_gradients_in->raw_data_in((cvm::real*)(&msg_data[0]));
      z_samples_in->raw_data_in((size_t*)(&msg_data[samp_start]));

      // Combine the new gradient and count of the other replicas
      // with Replica 0's current state
      global_z_gradients->add_grid(*z_gradients_in);
      global_z_samples->add_grid(*z_samples_in);
    }
  } else {
    // All other replicas send their current z gradient and z count.
    z_gradients->raw_data_out((cvm::real*)(&msg_data[0]));
    z_samples->raw_data_out((size_t*)(&msg_data[samp_start]));
    if (proxy->replica_comm_send(msg_data, msg_total, 0) != msg_total) {
      cvm::error("Error sending shared ABF data to replica.");
      return COLVARS_ERROR;
    }
  }

  // Without a barrier it's possible that one replica starts
  // share 2 when other replicas haven't finished share 1.
  proxy->replica_comm_barrier();
  // Done syncing the replicas.
  delete[] msg_data;

  return COLVARS_OK;
}


  // *****************************
  // ******  I/O functions  ******
  // *****************************


size_t colvarbias_abf::replica_share_freq() const
{
  return shared_freq;
}


template <class T> int colvarbias_abf::write_grid_to_file(T const *grid,
                                                          std::string const &filename,
                                                          bool close) {
  std::ostream &os = cvm::proxy->output_stream(filename, "multicolumn grid file");
  if (!os) {
    return cvm::error("Error opening file " + filename + " for writing.\n", COLVARS_ERROR | COLVARS_FILE_ERROR);
  }
  grid->write_multicol(os);
  if (close) {
    cvm::proxy->close_output_stream(filename);
  } else {
    // Insert empty line between frames in history files
    os << std::endl;
    cvm::proxy->flush_output_stream(filename);
  }

  // In dimension higher than 2, dx is easier to handle and visualize
  // but we cannot write multiple frames in a dx file now
  // (could be implemented as multiple dx files)
  if (num_variables() > 2 && close) {
    std::string  dx = filename + ".dx";
    std::ostream &dx_os = cvm::proxy->output_stream(dx, "OpenDX grid file");
    if (!dx_os)  {
      return cvm::error("Error opening file " + dx + " for writing.\n", COLVARS_ERROR | COLVARS_FILE_ERROR);
    }
    grid->write_opendx(dx_os);
    // if (close) {
      cvm::proxy->close_output_stream(dx);
    // }
    // else {
    //   // TODO, decide convention for multiple datasets in dx file
    //   *dx_os << std::endl;
    //   dx_os->flush();
    // }
  }
  return COLVARS_OK;
}


void colvarbias_abf::write_gradients_samples(const std::string &prefix, bool close, bool local)
{
  colvarproxy *proxy = cvm::main()->proxy;

  // The following are local aliases for the class' unique pointers
  colvar_grid_count *samples_out, *z_samples_out;
  colvar_grid_gradient *gradients_out, *z_gradients_out, *czar_gradients_out;
  integrate_potential *pmf_out, *czar_pmf_out;

  // In shared ABF, write grids containing local data only if requested
  if (local) {
    samples_out = local_samples.get();
    gradients_out = local_gradients.get();
    pmf_out = local_pmf.get();
    // Update the divergence before integrating the local PMF below
    // only needs to happen here, just before output
    local_pmf->set_div();
    z_samples_out = z_samples.get();
    z_gradients_out = z_gradients.get();
    czar_gradients_out = czar_gradients.get();
    czar_pmf_out = czar_pmf.get();
  } else {
    samples_out = samples.get();
    gradients_out = gradients.get();
    pmf_out = pmf.get();
    // Note: outside of shared ABF, "global" CZAR grids are just the local ones
    z_samples_out = global_z_samples.get();
    z_gradients_out = global_z_gradients.get();
    czar_gradients_out = global_czar_gradients.get();
    czar_pmf_out = global_czar_pmf.get();
  }

  write_grid_to_file<colvar_grid_count>(samples_out, prefix + ".count", close);
  write_grid_to_file<colvar_grid_gradient>(gradients_out, prefix + ".grad", close);

  if (b_integrate) {
    // Do numerical integration (to high precision) and output a PMF
    cvm::real err;
    // Divergence has already been updated on the fly for the global PMF (member data 'pmf')
    pmf_out->integrate(integrate_iterations, integrate_tol, err);
    pmf_out->set_zero_minimum();
    write_grid_to_file<colvar_grid_scalar>(pmf_out, prefix + ".pmf", close);
  }

  if (b_CZAR_estimator) {
    // Write eABF CZAR-related quantities
    write_grid_to_file<colvar_grid_count>(z_samples_out, prefix + ".zcount", close);
    if (b_czar_window_file) {
      write_grid_to_file<colvar_grid_gradient>(z_gradients_out, prefix + ".zgrad", close);
    }

    // Update the CZAR estimator of gradients, except at step 0
    // in which case we preserve any existing data (e.g. read via inputPrefix, used to join strata in stratified eABF)
    if (cvm::step_relative() > 0) {
      for (std::vector<int> iz_bin = czar_gradients_out->new_index();
            czar_gradients_out->index_ok(iz_bin); czar_gradients_out->incr(iz_bin)) {
        for (size_t n = 0; n < czar_gradients_out->multiplicity(); n++) {
          czar_gradients_out->set_value(iz_bin, z_gradients_out->value_output(iz_bin, n)
            - proxy->target_temperature() * proxy->boltzmann() * z_samples_out->log_gradient_finite_diff(iz_bin, n), n);
        }
      }
    }
    write_grid_to_file<colvar_grid_gradient>(czar_gradients_out, prefix + ".czar.grad", close);

    if (b_integrate) {
      // Do numerical integration (to high precision) and output a PMF
      cvm::real err;
      czar_pmf_out->set_div();
      czar_pmf_out->integrate(integrate_iterations, integrate_tol, err);
      czar_pmf_out->set_zero_minimum();
      write_grid_to_file<colvar_grid_scalar>(czar_pmf_out, prefix + ".czar.pmf", close);
    }
  }
  return;
}


// For Tcl implementation of selection rules.
/// Give the total number of bins for a given bias.
int colvarbias_abf::bin_num() {
  return samples->number_of_points();
}

/// Calculate the bin index for a given bias.
int colvarbias_abf::current_bin() {
  return samples->current_bin_flat_bound();
}

/// Give the count at a given bin index.
int colvarbias_abf::bin_count(int bin_index) {
  if (bin_index < 0 || bin_index >= bin_num()) {
    cvm::error("Error: Tried to get bin count from invalid bin index "+cvm::to_str(bin_index));
    return -1;
  }
  return int(samples->get_value(bin_index));
}

// Return the average number of samples in a given "radius" around current bin
int colvarbias_abf::colvarbias_abf::local_sample_count(int radius) {
  return samples->local_sample_count(radius);
}

int colvarbias_abf::read_gradients_samples()
{
  int err = COLVARS_OK;
  // Reading the CZAR gradients is necessary for joining strata in stratified eABF
  std::unique_ptr<colvar_grid_gradient> czar_gradients_in;

  if (b_CZAR_estimator) {
    // CZAR gradients are usually computed as needed from z-gradients and z_samples
    // Therefore the czar_gradients grid is not linked to a sampling grid
    // Here we define a temporary czar_gradients grid linked to z_samples,
    // to correctly average input gradients if overlapping
    czar_gradients_in.reset(new colvar_grid_gradient(colvars, z_samples));
  }

  for ( size_t i = 0; i < input_prefix.size(); i++ ) {
    std::string prefix = input_prefix[i];

    // For user-provided files, the per-bias naming scheme may not apply
    err |= samples->read_multicol(prefix + ".count", "ABF samples file", true);
    err |= gradients->read_multicol(prefix + ".grad", "ABF gradient file", true);

    if (shared_on) {
      last_gradients->copy_grid(*gradients);
      last_samples->copy_grid(*samples);
    }
    if (b_CZAR_estimator) {
      // Read eABF z-averaged data for CZAR
      err |= z_samples->read_multicol(prefix + ".zcount", "eABF z-histogram file", true);
      err |= z_gradients->read_multicol(prefix + ".zgrad", "eABF z-gradient file", true);
      err |= czar_gradients_in->read_multicol(prefix + ".czar.grad", "eABF CZAR gradient file", true);
    }
  }

  if (b_CZAR_estimator) {
    // Now copy real CZAR gradients (divided by total count) to the final grid
    for (std::vector<int> ix = czar_gradients->new_index();
          czar_gradients->index_ok(ix); czar_gradients->incr(ix)) {
      for (size_t n = 0; n < czar_gradients->multiplicity(); n++) {
        czar_gradients->set_value(ix, czar_gradients_in->value_output(ix, n), n);
      }
    }
  }
  return err;
}


template <typename OST> OST & colvarbias_abf::write_state_data_template_(OST &os)
{
  auto flags = os.flags();

  os.setf(std::ios::fmtflags(0), std::ios::floatfield); // default floating-point format

  write_state_data_key(os, "samples");
  samples->write_raw(os, 8);

  write_state_data_key(os, "gradient");
  gradients->write_raw(os, 8);

  if (shared_on) {
    write_state_data_key(os, "local_samples");
    local_samples->write_raw(os, 8);
    write_state_data_key(os, "local_gradient");
    local_gradients->write_raw(os, 8);
  }

  if (b_CZAR_estimator) {
    os.setf(std::ios::fmtflags(0), std::ios::floatfield); // default floating-point format
    write_state_data_key(os, "z_samples");
    z_samples->write_raw(os, 8);
    write_state_data_key(os, "z_gradient");
    z_gradients->write_raw(os, 8);
  }

  os.flags(flags);
  return os;
}


std::ostream & colvarbias_abf::write_state_data(std::ostream& os)
{
  return write_state_data_template_<std::ostream>(os);
}


cvm::memory_stream & colvarbias_abf::write_state_data(cvm::memory_stream& os)
{
  return write_state_data_template_<cvm::memory_stream>(os);
}


template <typename IST> IST &colvarbias_abf::read_state_data_template_(IST &is)
{
  if ( input_prefix.size() > 0 ) {
    cvm::error("ERROR: cannot provide both inputPrefix and a colvars state file.\n", COLVARS_INPUT_ERROR);
  }

  if (! read_state_data_key(is, "samples")) {
    return is;
  }
  if (! samples->read_raw(is)) {
    return is;
  }

  if (! read_state_data_key(is, "gradient")) {
    return is;
  }
  if (! gradients->read_raw(is)) {
    return is;
  }
  if (b_integrate) {
    // Update divergence to account for restart data
    pmf->set_div();
  }

  if (shared_on) {
    if (! read_state_data_key(is, "local_samples")) {
      return is;
    }
    if (! local_samples->read_raw(is)) {
      return is;
    }
    if (! read_state_data_key(is, "local_gradient")) {
      return is;
    }
    if (! local_gradients->read_raw(is)) {
      return is;
    }
  }

  if (b_CZAR_estimator) {

    if (! read_state_data_key(is, "z_samples")) {
      return is;
    }
    if (! z_samples->read_raw(is)) {
      return is;
    }

    if (! read_state_data_key(is, "z_gradient")) {
      return is;
    }
    if (! z_gradients->read_raw(is)) {
      return is;
    }
  }

  // Last samples / gradients must be updated after restart
  // reproducing the state after the last sharing step of previous run
  if (shared_on) {
    last_gradients->copy_grid(*gradients);
    last_samples->copy_grid(*samples);
    shared_last_step = cvm::step_absolute();
  }

  return is;
}


std::istream & colvarbias_abf::read_state_data(std::istream& is)
{
  return read_state_data_template_<std::istream>(is);
}


cvm::memory_stream & colvarbias_abf::read_state_data(cvm::memory_stream& is)
{
  return read_state_data_template_<cvm::memory_stream>(is);
}


int colvarbias_abf::write_output_files()
{
  if (cvm::debug()) {
    cvm::log("ABF bias trying to write gradients and samples to disk");
  }

  // In shared eABF/CZAR, the communication routine needs to run on all ranks
  if (shared_on) {
    if (b_CZAR_estimator) replica_share_CZAR();
  }

  // In shared setting, output local data for all replicas
  if (shared_on) {
    // Write local data on all replicas
    write_gradients_samples(output_prefix, true, true);

    if (cvm::main()->proxy->replica_index() > 0) {
      // No need to report the same data as replica 0, let it do the I/O job
      return COLVARS_OK;
    }
  }
  // In shared ABF, only replica 0 reaches this
  // filename prefix for master replica
  // used in mwABF to distinguish local from complete data
  std::string master_prefix = (shared_on ? output_prefix + ".all" : output_prefix);
  write_gradients_samples(master_prefix);

  if ((history_freq > 0) &&
      (!shared_on || cvm::main()->proxy->replica_index() == 0) && // if shared, only on replica 0
      (cvm::step_absolute() % history_freq == 0) &&               // at requested frequency
      (cvm::step_absolute() != history_last_step)) {              // not twice the same timestep
    write_gradients_samples(master_prefix + ".hist", false);
    history_last_step = cvm::step_absolute();
  }

  if (b_UI_estimator) {
    eabf_UI.calc_pmf();
    eabf_UI.write_files();
  }

  return COLVARS_OK;
}


int colvarbias_abf::calc_energy(std::vector<colvarvalue> const *values)
{
  bias_energy = 0.0; // default value, overridden if a value can be calculated

  if (num_variables() > 1 || values != NULL) {
    // Use simple estimate: neglect effect of fullSamples,
    // return value at center of bin
    if (pmf) {
      std::vector<int> curr_bin = values ?
        pmf->get_colvars_index(*values) :
        pmf->get_colvars_index();
      pmf->set_zero_minimum();
      pmf->wrap_to_edge(curr_bin, curr_bin); // Closest edge if outside of grid
      bias_energy = pmf->value(curr_bin);
    }
    return COLVARS_OK;
  }

  // Get the home bin.
  int home0 = gradients->current_bin_scalar(0);
  if (home0 < 0) return COLVARS_OK;
  int gradient_len = (int)(gradients->number_of_points(0));
  int home = (home0 < gradient_len) ? home0 : (gradient_len-1);

  // Integrate the gradient up to the home bin.
  cvm::real sum = 0.0;
  for (int i = 0; i < home; i++) {
    std::vector<int> ix(1,i);
    // Include the smoothing factor if necessary.
    sum += gradients->value_output_smoothed(ix, true) * gradients->widths[0];
  }

  // Integrate the gradient up to the current position in the home interval, a fractional portion of a bin.
  std::vector<int> ix(1,home);
  cvm::real frac = gradients->current_bin_scalar_fraction(0);
  sum += gradients->value_output_smoothed(ix, true) * gradients->widths[0] * frac;

  // The applied potential is the negative integral of force samples.
  bias_energy = -sum;
  return COLVARS_OK;
}
