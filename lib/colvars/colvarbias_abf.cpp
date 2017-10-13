// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarmodule.h"
#include "colvar.h"
#include "colvarbias_abf.h"


colvarbias_abf::colvarbias_abf(char const *key)
  : colvarbias(key),
    b_UI_estimator(false),
    b_CZAR_estimator(false),
    system_force(NULL),
    gradients(NULL),
    samples(NULL),
    z_gradients(NULL),
    z_samples(NULL),
    czar_gradients(NULL),
    last_gradients(NULL),
    last_samples(NULL)
{
}


int colvarbias_abf::init(std::string const &conf)
{
  colvarbias::init(conf);

  enable(f_cvb_scalar_variables);
  enable(f_cvb_calc_pmf);

  // TODO relax this in case of VMD plugin
  if (cvm::temperature() == 0.0)
    cvm::log("WARNING: ABF should not be run without a thermostat or at 0 Kelvin!\n");

  // ************* parsing general ABF options ***********************

  get_keyval_feature((colvarparse *)this, conf, "applyBias",  f_cvb_apply_force, true);
  if (!is_enabled(f_cvb_apply_force)){
    cvm::log("WARNING: ABF biases will *not* be applied!\n");
  }

  get_keyval(conf, "updateBias",  update_bias, true);
  if (update_bias) {
    enable(f_cvb_history_dependent);
  } else {
    cvm::log("WARNING: ABF biases will *not* be updated!\n");
  }

  get_keyval(conf, "hideJacobian", hide_Jacobian, false);
  if (hide_Jacobian) {
    cvm::log("Jacobian (geometric) forces will be handled internally.\n");
  } else {
    cvm::log("Jacobian (geometric) forces will be included in reported free energy gradients.\n");
  }

  get_keyval(conf, "fullSamples", full_samples, 200);
  if ( full_samples <= 1 ) full_samples = 1;
  min_samples = full_samples / 2;
  // full_samples - min_samples >= 1 is guaranteed

  get_keyval(conf, "inputPrefix",  input_prefix, std::vector<std::string>());
  get_keyval(conf, "outputFreq", output_freq, cvm::restart_out_freq);
  get_keyval(conf, "historyFreq", history_freq, 0);
  b_history_files = (history_freq > 0);

  // shared ABF
  get_keyval(conf, "shared", shared_on, false);
  if (shared_on) {
    if (!cvm::replica_enabled() || cvm::replica_num() <= 1) {
      cvm::error("Error: shared ABF requires more than one replica.");
      return COLVARS_ERROR;
    }
    cvm::log("shared ABF will be applied among "+ cvm::to_str(cvm::replica_num()) + " replicas.\n");
    if (cvm::proxy->smp_enabled() == COLVARS_OK) {
      cvm::error("Error: shared ABF is currently not available with SMP parallelism; "
                 "please set \"SMP off\" at the top of the Colvars configuration file.\n",
                 COLVARS_NOT_IMPLEMENTED);
      return COLVARS_NOT_IMPLEMENTED;
    }

    // If shared_freq is not set, we default to output_freq
    get_keyval(conf, "sharedFreq", shared_freq, output_freq);
  }

  // ************* checking the associated colvars *******************

  if (colvars.size() == 0) {
    cvm::error("Error: no collective variables specified for the ABF bias.\n");
    return COLVARS_ERROR;
  }

  if (update_bias) {
    // Request calculation of total force
    if(enable(f_cvb_get_total_force)) return cvm::get_error();
  }

  bool b_extended = false;
  for (size_t i = 0; i < colvars.size(); i++) {

    if (colvars[i]->value().type() != colvarvalue::type_scalar) {
      cvm::error("Error: ABF bias can only use scalar-type variables.\n");
    }
    colvars[i]->enable(f_cv_grid);
    if (hide_Jacobian) {
      colvars[i]->enable(f_cv_hide_Jacobian);
    }

    // If any colvar is extended-system, we need to collect the extended
    // system gradient
    if (colvars[i]->is_enabled(f_cv_extended_Lagrangian))
      b_extended = true;

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

  if (get_keyval(conf, "maxForce", max_force)) {
    if (max_force.size() != colvars.size()) {
      cvm::error("Error: Number of parameters to maxForce does not match number of colvars.");
    }
    for (size_t i = 0; i < colvars.size(); i++) {
      if (max_force[i] < 0.0) {
        cvm::error("Error: maxForce should be non-negative.");
      }
    }
    cap_force = true;
  } else {
    cap_force = false;
  }

  bin.assign(colvars.size(), 0);
  force_bin.assign(colvars.size(), 0);
  system_force = new cvm::real [colvars.size()];

  // Construct empty grids based on the colvars
  if (cvm::debug()) {
    cvm::log("Allocating count and free energy gradient grids.\n");
  }

  samples   = new colvar_grid_count(colvars);
  gradients = new colvar_grid_gradient(colvars);
  gradients->samples = samples;
  samples->has_parent_data = true;

  // Data for eABF z-based estimator
  if (b_extended) {
    get_keyval(conf, "CZARestimator", b_CZAR_estimator, true);
    // CZAR output files for stratified eABF
    get_keyval(conf, "writeCZARwindowFile", b_czar_window_file, false,
               colvarparse::parse_silent);

    z_bin.assign(colvars.size(), 0);
    z_samples   = new colvar_grid_count(colvars);
    z_samples->request_actual_value();
    z_gradients = new colvar_grid_gradient(colvars);
    z_gradients->request_actual_value();
    z_gradients->samples = z_samples;
    z_samples->has_parent_data = true;
    czar_gradients = new colvar_grid_gradient(colvars);
  }

  // For shared ABF, we store a second set of grids.
  // This used to be only if "shared" was defined,
  // but now we allow calling share externally (e.g. from Tcl).
  last_samples   = new colvar_grid_count(colvars);
  last_gradients = new colvar_grid_gradient(colvars);
  last_gradients->samples = last_samples;
  last_samples->has_parent_data = true;
  shared_last_step = -1;

  // If custom grids are provided, read them
  if ( input_prefix.size() > 0 ) {
    read_gradients_samples();
  }

  // if extendedLangrangian is on, then call UI estimator
  if (b_extended) {
    get_keyval(conf, "UIestimator", b_UI_estimator, false);

    if (b_UI_estimator) {
    std::vector<double> UI_lowerboundary;
    std::vector<double> UI_upperboundary;
    std::vector<double> UI_width;
    std::vector<double> UI_krestr;

    bool UI_restart = (input_prefix.size() > 0);

    for (size_t i = 0; i < colvars.size(); i++)
    {
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
                                         cvm::temperature());
    }
  }

  cvm::log("Finished ABF setup.\n");
  return COLVARS_OK;
}

/// Destructor
colvarbias_abf::~colvarbias_abf()
{
  if (samples) {
    delete samples;
    samples = NULL;
  }

  if (gradients) {
    delete gradients;
    gradients = NULL;
  }

  if (z_samples) {
    delete z_samples;
    z_samples = NULL;
  }

  if (z_gradients) {
    delete z_gradients;
    z_gradients = NULL;
  }

  if (czar_gradients) {
    delete czar_gradients;
    czar_gradients = NULL;
  }

  // shared ABF
  // We used to only do this if "shared" was defined,
  // but now we can call shared externally
  if (last_samples) {
    delete last_samples;
    last_samples = NULL;
  }

  if (last_gradients) {
    delete last_gradients;
    last_gradients = NULL;
  }

  if (system_force) {
    delete [] system_force;
    system_force = NULL;
  }
}


/// Update the FE gradient, compute and apply biasing force
/// also output data to disk if needed

int colvarbias_abf::update()
{
  if (cvm::debug()) cvm::log("Updating ABF bias " + this->name);

  if (cvm::step_relative() == 0) {

    // At first timestep, do only:
    // initialization stuff (file operations relying on n_abf_biases
    // compute current value of colvars

    for (size_t i = 0; i < colvars.size(); i++) {
      bin[i] = samples->current_bin_scalar(i);
    }

  } else {

    for (size_t i = 0; i < colvars.size(); i++) {
      bin[i] = samples->current_bin_scalar(i);
    }

    if ( update_bias && samples->index_ok(force_bin) ) {
      // Only if requested and within bounds of the grid...

      for (size_t i = 0; i < colvars.size(); i++) {
        // get total forces (lagging by 1 timestep) from colvars
        // and subtract previous ABF force if necessary
        update_system_force(i);
      }
      gradients->acc_force(force_bin, system_force);
    }
    if ( z_gradients && update_bias ) {
      for (size_t i = 0; i < colvars.size(); i++) {
        z_bin[i] = z_samples->current_bin_scalar(i);
      }
      if ( z_samples->index_ok(z_bin) ) {
        for (size_t i = 0; i < colvars.size(); i++) {
          // If we are outside the range of xi, the force has not been obtained above
          // the function is just an accessor, so cheap to call again anyway
          update_system_force(i);
        }
        z_gradients->acc_force(z_bin, system_force);
      }
    }
  }

  // save bin for next timestep
  force_bin = bin;

  // Reset biasing forces from previous timestep
  for (size_t i = 0; i < colvars.size(); i++) {
    colvar_forces[i].reset();
  }

  // Compute and apply the new bias, if applicable
  if (is_enabled(f_cvb_apply_force) && samples->index_ok(bin)) {

    size_t count = samples->value(bin);
    cvm::real fact = 1.0;

    // Factor that ensures smooth introduction of the force
    if ( count < full_samples ) {
      fact = (count < min_samples) ? 0.0 :
        (cvm::real(count - min_samples)) / (cvm::real(full_samples - min_samples));
    }

    const cvm::real * grad  = &(gradients->value(bin));

    if ( fact != 0.0 ) {
      if ( (colvars.size() == 1) && colvars[0]->periodic_boundaries() ) {
        // Enforce a zero-mean bias on periodic, 1D coordinates
        // in other words: boundary condition is that the biasing potential is periodic
        colvar_forces[0].real_value = fact * (grad[0] / cvm::real(count) - gradients->average());
      } else {
        for (size_t i = 0; i < colvars.size(); i++) {
          // subtracting the mean force (opposite of the FE gradient) means adding the gradient
          colvar_forces[i].real_value = fact * grad[i] / cvm::real(count);
        }
      }
      if (cap_force) {
        for (size_t i = 0; i < colvars.size(); i++) {
          if ( colvar_forces[i].real_value * colvar_forces[i].real_value > max_force[i] * max_force[i] ) {
            colvar_forces[i].real_value = (colvar_forces[i].real_value > 0 ? max_force[i] : -1.0 * max_force[i]);
          }
        }
      }
    }
  }

  // update the output prefix; TODO: move later to setup_output() function
  if (cvm::main()->num_biases_feature(colvardeps::f_cvb_calc_pmf) == 1) {
    // This is the only bias computing PMFs
    output_prefix = cvm::output_prefix();
  } else {
    output_prefix = cvm::output_prefix() + "." + this->name;
  }

  if (output_freq && (cvm::step_absolute() % output_freq) == 0) {
    if (cvm::debug()) cvm::log("ABF bias trying to write gradients and samples to disk");
    write_gradients_samples(output_prefix);
  }

  if (b_history_files && (cvm::step_absolute() % history_freq) == 0) {
    // file already exists iff cvm::step_relative() > 0
    // otherwise, backup and replace
    write_gradients_samples(output_prefix + ".hist", (cvm::step_relative() > 0));
  }

  if (shared_on && shared_last_step >= 0 && cvm::step_absolute() % shared_freq == 0) {
    // Share gradients and samples for shared ABF.
    replica_share();
  }

  // Prepare for the first sharing.
  if (shared_last_step < 0) {
    // Copy the current gradient and count values into last.
    last_gradients->copy_grid(*gradients);
    last_samples->copy_grid(*samples);
    shared_last_step = cvm::step_absolute();
    cvm::log("Prepared sample and gradient buffers at step "+cvm::to_str(cvm::step_absolute())+".");
  }

  // update UI estimator every step
  if (b_UI_estimator)
  {
    std::vector<double> x(colvars.size(),0);
    std::vector<double> y(colvars.size(),0);
    for (size_t i = 0; i < colvars.size(); i++)
    {
      x[i] = colvars[i]->actual_value();
      y[i] = colvars[i]->value();
    }
    eabf_UI.update_output_filename(output_prefix);
    eabf_UI.update(cvm::step_absolute(), x, y);
  }

  return COLVARS_OK;
}


int colvarbias_abf::replica_share() {
  int p;

  if ( !cvm::replica_enabled() ) {
    cvm::error("Error: shared ABF: No replicas.\n");
    return COLVARS_ERROR;
  }
  // We must have stored the last_gradients and last_samples.
  if (shared_last_step < 0 ) {
    cvm::error("Error: shared ABF: Tried to apply shared ABF before any sampling had occurred.\n");
    return COLVARS_ERROR;
  }

  // Share gradients for shared ABF.
  cvm::log("shared ABF: Sharing gradient and samples among replicas at step "+cvm::to_str(cvm::step_absolute()) );

  // Count of data items.
  size_t data_n = gradients->raw_data_num();
  size_t samp_start = data_n*sizeof(cvm::real);
  size_t msg_total = data_n*sizeof(size_t) + samp_start;
  char* msg_data = new char[msg_total];

  if (cvm::replica_index() == 0) {
    // Replica 0 collects the delta gradient and count from the others.
    for (p = 1; p < cvm::replica_num(); p++) {
      // Receive the deltas.
      cvm::replica_comm_recv(msg_data, msg_total, p);

      // Map the deltas from the others into the grids.
      last_gradients->raw_data_in((cvm::real*)(&msg_data[0]));
      last_samples->raw_data_in((size_t*)(&msg_data[samp_start]));

      // Combine the delta gradient and count of the other replicas
      // with Replica 0's current state (including its delta).
      gradients->add_grid( *last_gradients );
      samples->add_grid( *last_samples );
    }

    // Now we must send the combined gradient to the other replicas.
    gradients->raw_data_out((cvm::real*)(&msg_data[0]));
    samples->raw_data_out((size_t*)(&msg_data[samp_start]));
    for (p = 1; p < cvm::replica_num(); p++) {
      cvm::replica_comm_send(msg_data, msg_total, p);
    }

  } else {
    // All other replicas send their delta gradient and count.
    // Calculate the delta gradient and count.
    last_gradients->delta_grid(*gradients);
    last_samples->delta_grid(*samples);

    // Cast the raw char data to the gradient and samples.
    last_gradients->raw_data_out((cvm::real*)(&msg_data[0]));
    last_samples->raw_data_out((size_t*)(&msg_data[samp_start]));
    cvm::replica_comm_send(msg_data, msg_total, 0);

    // We now receive the combined gradient from Replica 0.
    cvm::replica_comm_recv(msg_data, msg_total, 0);
    // We sync to the combined gradient computed by Replica 0.
    gradients->raw_data_in((cvm::real*)(&msg_data[0]));
    samples->raw_data_in((size_t*)(&msg_data[samp_start]));
  }

  // Without a barrier it's possible that one replica starts
  // share 2 when other replicas haven't finished share 1.
  cvm::replica_comm_barrier();
  // Done syncing the replicas.
  delete[] msg_data;

  // Copy the current gradient and count values into last.
  last_gradients->copy_grid(*gradients);
  last_samples->copy_grid(*samples);
  shared_last_step = cvm::step_absolute();

  return COLVARS_OK;
}

void colvarbias_abf::write_gradients_samples(const std::string &prefix, bool append)
{
  std::string  samples_out_name = prefix + ".count";
  std::string  gradients_out_name = prefix + ".grad";
  std::ios::openmode mode = (append ? std::ios::app : std::ios::out);

  std::ostream *samples_os =
    cvm::proxy->output_stream(samples_out_name, mode);
  if (!samples_os) {
    cvm::error("Error opening ABF samples file " + samples_out_name + " for writing");
  }
  samples->write_multicol(*samples_os);
  cvm::proxy->close_output_stream(samples_out_name);

  std::ostream *gradients_os =
    cvm::proxy->output_stream(gradients_out_name, mode);
  if (!gradients_os) {
    cvm::error("Error opening ABF gradient file " + gradients_out_name + " for writing");
  }
  gradients->write_multicol(*gradients_os);
  cvm::proxy->close_output_stream(gradients_out_name);

  if (colvars.size() == 1) {
    // Do numerical integration and output a PMF
    std::string  pmf_out_name = prefix + ".pmf";
    std::ostream *pmf_os = cvm::proxy->output_stream(pmf_out_name, mode);
    if (!pmf_os) {
      cvm::error("Error opening pmf file " + pmf_out_name + " for writing");
    }
    gradients->write_1D_integral(*pmf_os);
    *pmf_os << std::endl;
    cvm::proxy->close_output_stream(pmf_out_name);
  }

  if (b_CZAR_estimator) {
    // Write eABF CZAR-related quantities

    std::string  z_samples_out_name = prefix + ".zcount";

    std::ostream *z_samples_os =
      cvm::proxy->output_stream(z_samples_out_name, mode);
    if (!z_samples_os) {
      cvm::error("Error opening eABF z-histogram file " + z_samples_out_name + " for writing");
    }
    z_samples->write_multicol(*z_samples_os);
    cvm::proxy->close_output_stream(z_samples_out_name);

    if (b_czar_window_file) {
      std::string  z_gradients_out_name = prefix + ".zgrad";

      std::ostream *z_gradients_os =
        cvm::proxy->output_stream(z_gradients_out_name, mode);
      if (!z_gradients_os) {
        cvm::error("Error opening eABF z-gradient file " + z_gradients_out_name + " for writing");
      }
      z_gradients->write_multicol(*z_gradients_os);
      cvm::proxy->close_output_stream(z_gradients_out_name);
    }

    // Calculate CZAR estimator of gradients
    for (std::vector<int> ix = czar_gradients->new_index();
          czar_gradients->index_ok(ix); czar_gradients->incr(ix)) {
      for (size_t n = 0; n < czar_gradients->multiplicity(); n++) {
        czar_gradients->set_value(ix, z_gradients->value_output(ix, n)
            - cvm::temperature() * cvm::boltzmann() * z_samples->log_gradient_finite_diff(ix, n),
            n);
      }
    }

    std::string  czar_gradients_out_name = prefix + ".czar.grad";

    std::ostream *czar_gradients_os =
      cvm::proxy->output_stream(czar_gradients_out_name, mode);
    if (!czar_gradients_os) {
      cvm::error("Error opening CZAR gradient file " + czar_gradients_out_name + " for writing");
    }
    czar_gradients->write_multicol(*czar_gradients_os);
    cvm::proxy->close_output_stream(czar_gradients_out_name);

    if (colvars.size() == 1) {
      // Do numerical integration and output a PMF
      std::string  czar_pmf_out_name = prefix + ".czar.pmf";
      std::ostream *czar_pmf_os =
        cvm::proxy->output_stream(czar_pmf_out_name, mode);
      if (!czar_pmf_os)  cvm::error("Error opening CZAR pmf file " + czar_pmf_out_name + " for writing");
      czar_gradients->write_1D_integral(*czar_pmf_os);
      *czar_pmf_os << std::endl;
      cvm::proxy->close_output_stream(czar_pmf_out_name);
    }
  }
  return;
}


// For Tcl implementation of selection rules.
/// Give the total number of bins for a given bias.
int colvarbias_abf::bin_num() {
  return samples->number_of_points(0);
}
/// Calculate the bin index for a given bias.
int colvarbias_abf::current_bin() {
  return samples->current_bin_scalar(0);
}
/// Give the count at a given bin index.
int colvarbias_abf::bin_count(int bin_index) {
  if (bin_index < 0 || bin_index >= bin_num()) {
    cvm::error("Error: Tried to get bin count from invalid bin index "+cvm::to_str(bin_index));
    return -1;
  }
  std::vector<int> ix(1,(int)bin_index);
  return samples->value(ix);
}


void colvarbias_abf::read_gradients_samples()
{
  std::string samples_in_name, gradients_in_name, z_samples_in_name, z_gradients_in_name;

  for ( size_t i = 0; i < input_prefix.size(); i++ ) {
    samples_in_name = input_prefix[i] + ".count";
    gradients_in_name = input_prefix[i] + ".grad";
    z_samples_in_name = input_prefix[i] + ".zcount";
    z_gradients_in_name = input_prefix[i] + ".zgrad";
    // For user-provided files, the per-bias naming scheme may not apply

    std::ifstream is;

    cvm::log("Reading sample count from " + samples_in_name + " and gradient from " + gradients_in_name);
    is.open(samples_in_name.c_str());
    if (!is.is_open()) cvm::error("Error opening ABF samples file " + samples_in_name + " for reading");
    samples->read_multicol(is, true);
    is.close();
    is.clear();

    is.open(gradients_in_name.c_str());
    if (!is.is_open()) {
      cvm::error("Error opening ABF gradient file " +
                 gradients_in_name + " for reading", INPUT_ERROR);
    } else {
      gradients->read_multicol(is, true);
      is.close();
    }

    if (b_CZAR_estimator) {
      // Read eABF z-averaged data for CZAR
      cvm::log("Reading z-histogram from " + z_samples_in_name + " and z-gradient from " + z_gradients_in_name);

      is.clear();
      is.open(z_samples_in_name.c_str());
      if (!is.is_open())  cvm::error("Error opening eABF z-histogram file " + z_samples_in_name + " for reading");
      z_samples->read_multicol(is, true);
      is.close();
      is.clear();

      is.open(z_gradients_in_name.c_str());
      if (!is.is_open())  cvm::error("Error opening eABF z-gradient file " + z_gradients_in_name + " for reading");
      z_gradients->read_multicol(is, true);
      is.close();
    }
  }
  return;
}


std::ostream & colvarbias_abf::write_state_data(std::ostream& os)
{
  std::ios::fmtflags flags(os.flags());

  os.setf(std::ios::fmtflags(0), std::ios::floatfield); // default floating-point format
  os << "\nsamples\n";
  samples->write_raw(os, 8);
  os.flags(flags);

  os << "\ngradient\n";
  gradients->write_raw(os, 8);

  if (b_CZAR_estimator) {
    os.setf(std::ios::fmtflags(0), std::ios::floatfield); // default floating-point format
    os << "\nz_samples\n";
    z_samples->write_raw(os, 8);
    os.flags(flags);
    os << "\nz_gradient\n";
    z_gradients->write_raw(os, 8);
  }

  os.flags(flags);
  return os;
}


std::istream & colvarbias_abf::read_state_data(std::istream& is)
{
  if ( input_prefix.size() > 0 ) {
    cvm::error("ERROR: cannot provide both inputPrefix and a colvars state file.\n", INPUT_ERROR);
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

  return is;
}
