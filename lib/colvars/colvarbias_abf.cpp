/********************************************************************************
 * Implementation of the ABF and histogram biases                               *
 ********************************************************************************/

#include "colvarmodule.h"
#include "colvar.h"
#include "colvarbias_abf.h"

/// ABF bias constructor; parses the config file

colvarbias_abf::colvarbias_abf (std::string const &conf, char const *key)
  : colvarbias (conf, key),
    gradients (NULL),
    samples (NULL)
{
  if (cvm::temperature() == 0.0)
    cvm::log ("WARNING: ABF should not be run without a thermostat or at 0 Kelvin!\n");

  // ************* parsing general ABF options ***********************

  get_keyval (conf, "applyBias",  apply_bias, true);
  if (!apply_bias) cvm::log ("WARNING: ABF biases will *not* be applied!\n");

  get_keyval (conf, "updateBias",  update_bias, true);
  if (!update_bias) cvm::log ("WARNING: ABF biases will *not* be updated!\n");

  get_keyval (conf, "hideJacobian", hide_Jacobian, false);
  if (hide_Jacobian) {
    cvm::log ("Jacobian (geometric) forces will be handled internally.\n");
  } else {
    cvm::log ("Jacobian (geometric) forces will be included in reported free energy gradients.\n");
  }

  get_keyval (conf, "fullSamples", full_samples, 200);
  if ( full_samples <= 1 ) full_samples = 1;
  min_samples = full_samples / 2;
  // full_samples - min_samples >= 1 is guaranteed

  get_keyval (conf, "inputPrefix",  input_prefix, std::vector<std::string> ());
  get_keyval (conf, "outputFreq", output_freq, cvm::restart_out_freq);
  get_keyval (conf, "historyFreq", history_freq, 0);
  b_history_files = (history_freq > 0);

  // ************* checking the associated colvars *******************

  if (colvars.size() == 0) {
    cvm::fatal_error ("Error: no collective variables specified for the ABF bias.\n");
  }

  for (size_t i = 0; i < colvars.size(); i++) {

    if (colvars[i]->type() != colvarvalue::type_scalar) {
      cvm::fatal_error ("Error: ABF bias can only use scalar-type variables.\n");
    }

    colvars[i]->enable (colvar::task_gradients);

    if (update_bias) {
      // Request calculation of system force (which also checks for availability)
      colvars[i]->enable (colvar::task_system_force);

      if (!colvars[i]->tasks[colvar::task_extended_lagrangian]) {
        // request computation of Jacobian force
        colvars[i]->enable (colvar::task_Jacobian_force);

        // request Jacobian force as part as system force
        // except if the user explicitly requires the "silent" Jacobian
        // correction AND the colvar has a single component
        if (hide_Jacobian) {
          if (colvars[i]->n_components() > 1) {
            cvm::log ("WARNING: colvar \"" + colvars[i]->name
            + "\" has multiple components; reporting its Jacobian forces\n");
            colvars[i]->enable (colvar::task_report_Jacobian_force);
          }
        } else {
          colvars[i]->enable (colvar::task_report_Jacobian_force);
        }
      }
    }

    // Here we could check for orthogonality of the Cartesian coordinates
    // and make it just a warning if some parameter is set?
  }

  if (get_keyval (conf, "maxForce", max_force)) {
    if (max_force.size() != colvars.size()) {
      cvm::fatal_error ("Error: Number of parameters to maxForce does not match number of colvars.");
    }
    for (size_t i=0; i<colvars.size(); i++) {
      if (max_force[i] < 0.0) {
        cvm::fatal_error ("Error: maxForce should be non-negative.");
      }
    }
    cap_force = true;
  } else {
    cap_force = false;
  }

  bin.assign (colvars.size(), 0);
  force_bin.assign (colvars.size(), 0);
  force = new cvm::real [colvars.size()];

  // Construct empty grids based on the colvars
  samples   = new colvar_grid_count    (colvars);
  gradients = new colvar_grid_gradient (colvars);
  gradients->samples = samples;
  samples->has_parent_data = true;

  // If custom grids are provided, read them
  if ( input_prefix.size() > 0 ) {
    read_gradients_samples ();
  }

  cvm::log ("Finished ABF setup.\n");
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

  delete [] force;

  if (cvm::n_abf_biases > 0)
    cvm::n_abf_biases -= 1;
}


/// Update the FE gradient, compute and apply biasing force
/// also output data to disk if needed

cvm::real colvarbias_abf::update()
{
  if (cvm::debug()) cvm::log ("Updating ABF bias " + this->name);

  if (cvm::step_relative() == 0) {

    // At first timestep, do only:
    // initialization stuff (file operations relying on n_abf_biases
    // compute current value of colvars

    if ( cvm::n_abf_biases == 1 && cvm::n_meta_biases == 0 ) {
      // This is the only ABF bias
      output_prefix = cvm::output_prefix;
    } else {
      output_prefix = cvm::output_prefix + "." + this->name;
    }

    for (size_t i=0; i<colvars.size(); i++) {
      bin[i] = samples->current_bin_scalar(i);
    }

  } else {

    for (size_t i=0; i<colvars.size(); i++) {
      bin[i] = samples->current_bin_scalar(i);
    }

    if ( update_bias && samples->index_ok (force_bin) ) {
      // Only if requested and within bounds of the grid...

      for (size_t i=0; i<colvars.size(); i++) {	  // get forces (lagging by 1 timestep) from colvars
        force[i] = colvars[i]->system_force();
      }
      gradients->acc_force (force_bin, force);
    }
  }

  // save bin for next timestep
  force_bin = bin;

  // Reset biasing forces from previous timestep
  for (size_t i=0; i<colvars.size(); i++) {
    colvar_forces[i].reset();
  }

  // Compute and apply the new bias, if applicable
  if ( apply_bias && samples->index_ok (bin) ) {

    size_t  count = samples->value (bin);
    cvm::real	fact = 1.0;

    // Factor that ensures smooth introduction of the force
    if ( count < full_samples ) {
      fact = ( count < min_samples) ? 0.0 :
        (cvm::real (count - min_samples)) / (cvm::real (full_samples - min_samples));
    }

    const cvm::real * grad  = &(gradients->value (bin));

    if ( fact != 0.0 ) {
      if ( (colvars.size() == 1) && colvars[0]->periodic_boundaries() ) {
        // Enforce a zero-mean bias on periodic, 1D coordinates
        // in other words: boundary condition is that the biasing potential is periodic
        colvar_forces[0].real_value = fact * (grad[0] / cvm::real (count) - gradients->average ());
      } else {
        for (size_t i=0; i<colvars.size(); i++) {
          // subtracting the mean force (opposite of the FE gradient) means adding the gradient
          colvar_forces[i].real_value = fact * grad[i] / cvm::real (count);
        }
      }
      if (cap_force) {
        for (size_t i=0; i<colvars.size(); i++) {
          if ( colvar_forces[i].real_value * colvar_forces[i].real_value > max_force[i] * max_force[i] ) {
            colvar_forces[i].real_value = (colvar_forces[i].real_value > 0 ? max_force[i] : -1.0 * max_force[i]);
          }
        }
      }
    }
  }

  if (output_freq && (cvm::step_absolute() % output_freq) == 0) {
    if (cvm::debug()) cvm::log ("ABF bias trying to write gradients and samples to disk");
    write_gradients_samples (output_prefix);
  }
  if (b_history_files && (cvm::step_absolute() % history_freq) == 0) {
    // append to existing file only if cvm::step_absolute() > 0
    // otherwise, backup and replace
    write_gradients_samples (output_prefix + ".hist", (cvm::step_absolute() > 0));
  }
  return 0.0;
}


void colvarbias_abf::write_gradients_samples (const std::string &prefix, bool append)
{
  std::string  samples_out_name = prefix + ".count";
  std::string  gradients_out_name = prefix + ".grad";
  std::ios::openmode mode = (append ? std::ios::app : std::ios::out);

  std::ofstream samples_os;
  std::ofstream gradients_os;

  if (!append) cvm::backup_file (samples_out_name.c_str());
  samples_os.open (samples_out_name.c_str(), mode);
  if (!samples_os.good()) cvm::fatal_error ("Error opening ABF samples file " + samples_out_name + " for writing");
  samples->write_multicol (samples_os);
  samples_os.close ();

  if (!append) cvm::backup_file (gradients_out_name.c_str());
  gradients_os.open (gradients_out_name.c_str(), mode);
  if (!gradients_os.good())	cvm::fatal_error ("Error opening ABF gradient file " + gradients_out_name + " for writing");
  gradients->write_multicol (gradients_os);
  gradients_os.close ();

  if (colvars.size () == 1) {
    std::string  pmf_out_name = prefix + ".pmf";
    if (!append) cvm::backup_file (pmf_out_name.c_str());
    std::ofstream pmf_os;
    // Do numerical integration and output a PMF
    pmf_os.open (pmf_out_name.c_str(), mode);
    if (!pmf_os.good())	cvm::fatal_error ("Error opening pmf file " + pmf_out_name + " for writing");
    gradients->write_1D_integral (pmf_os);
    pmf_os << std::endl;
    pmf_os.close ();
  }
  return;
}

void colvarbias_abf::read_gradients_samples ()
{
  std::string samples_in_name, gradients_in_name;

  for ( size_t i = 0; i < input_prefix.size(); i++ ) {
    samples_in_name = input_prefix[i] + ".count";
    gradients_in_name = input_prefix[i] + ".grad";
    // For user-provided files, the per-bias naming scheme may not apply

    std::ifstream is;

    cvm::log ("Reading sample count from " + samples_in_name + " and gradients from " + gradients_in_name);
    is.open (samples_in_name.c_str());
    if (!is.good()) cvm::fatal_error ("Error opening ABF samples file " + samples_in_name + " for reading");
    samples->read_multicol (is, true);
    is.close ();
    is.clear();

    is.open (gradients_in_name.c_str());
    if (!is.good())	cvm::fatal_error ("Error opening ABF gradient file " + gradients_in_name + " for reading");
    gradients->read_multicol (is, true);
    is.close ();
  }
  return;
}


std::ostream & colvarbias_abf::write_restart (std::ostream& os)
{

  std::ios::fmtflags flags (os.flags ());
  os.setf(std::ios::fmtflags (0), std::ios::floatfield); // default floating-point format

  os << "abf {\n"
     << "  configuration {\n"
     << "    name " << this->name << "\n";
  os << "  }\n";

  os << "samples\n";
  samples->write_raw (os, 8);

  os << "\ngradient\n";
  gradients->write_raw (os);

  os << "}\n\n";

  os.flags (flags);
  return os;
}


std::istream & colvarbias_abf::read_restart (std::istream& is)
{
  if ( input_prefix.size() > 0 ) {
    cvm::fatal_error ("ERROR: cannot provide both inputPrefix and restart information (colvarsInput)");
  }

  size_t const start_pos = is.tellg();

  cvm::log ("Restarting ABF bias \""+
            this->name+"\".\n");
  std::string key, brace, conf;

  if ( !(is >> key)   || !(key == "abf") ||
       !(is >> brace) || !(brace == "{") ||
       !(is >> colvarparse::read_block ("configuration", conf)) ) {
    cvm::log ("Error: in reading restart configuration for ABF bias \""+
              this->name+"\" at position "+
              cvm::to_str (is.tellg())+" in stream.\n");
    is.clear();
    is.seekg (start_pos, std::ios::beg);
    is.setstate (std::ios::failbit);
    return is;
  }

  std::string name = "";
  if ( (colvarparse::get_keyval (conf, "name", name, std::string (""), colvarparse::parse_silent)) &&
         (name != this->name) )
    cvm::fatal_error ("Error: in the restart file, the "
                      "\"abf\" block has wrong name (" + name + ")\n");
  if ( name == "" ) {
    cvm::fatal_error ("Error: \"abf\" block in the restart file has no name.\n");
  }

  if ( !(is >> key)   || !(key == "samples")) {
    cvm::log ("Error: in reading restart configuration for ABF bias \""+
              this->name+"\" at position "+
              cvm::to_str (is.tellg())+" in stream.\n");
    is.clear();
    is.seekg (start_pos, std::ios::beg);
    is.setstate (std::ios::failbit);
    return is;
  }
  if (! samples->read_raw (is)) {
    samples->read_raw_error();
  }

  if ( !(is >> key)   || !(key == "gradient")) {
    cvm::log ("Error: in reading restart configuration for ABF bias \""+
              this->name+"\" at position "+
              cvm::to_str (is.tellg())+" in stream.\n");
    is.clear();
    is.seekg (start_pos, std::ios::beg);
    is.setstate (std::ios::failbit);
    return is;
  }
  if (! gradients->read_raw (is)) {
    gradients->read_raw_error();
  }

  is >> brace;
  if (brace != "}") {
    cvm::fatal_error ("Error: corrupt restart information for ABF bias \""+
                      this->name+"\": no matching brace at position "+
                      cvm::to_str (is.tellg())+" in the restart file.\n");
    is.setstate (std::ios::failbit);
  }
  return is;
}




/// Histogram "bias" constructor

colvarbias_histogram::colvarbias_histogram (std::string const &conf, char const *key)
  : colvarbias (conf, key),
    grid (NULL)
{
  get_keyval (conf, "outputfreq", output_freq, cvm::restart_out_freq);

  if ( output_freq == 0 ) {
    cvm::fatal_error ("User required histogram with zero output frequency");
  }

  grid   = new colvar_grid_count    (colvars);
  bin.assign (colvars.size(), 0);

  out_name = cvm::output_prefix + "." + this->name + ".dat";
  cvm::log ("Histogram will be written to file " + out_name);

  cvm::log ("Finished histogram setup.\n");
}

/// Destructor
colvarbias_histogram::~colvarbias_histogram()
{
  if (grid_os.good())	grid_os.close();

  if (grid) {
    delete grid;
    grid = NULL;
  }

  if (cvm::n_histo_biases > 0)
    cvm::n_histo_biases -= 1;
}

/// Update the grid
cvm::real colvarbias_histogram::update()
{
  if (cvm::debug()) cvm::log ("Updating Grid bias " + this->name);

  for (size_t i=0; i<colvars.size(); i++) {
    bin[i] = grid->current_bin_scalar(i);
  }

  if ( grid->index_ok (bin) ) {	  // Only within bounds of the grid...
    grid->incr_count (bin);
  }

  if (output_freq && (cvm::step_absolute() % output_freq) == 0) {
    if (cvm::debug()) cvm::log ("Histogram bias trying to write grid to disk");

    grid_os.open (out_name.c_str());
    if (!grid_os.good()) cvm::fatal_error ("Error opening histogram file " + out_name + " for writing");
    grid->write_multicol (grid_os);
    grid_os.close ();
  }
  return 0.0; // no bias energy for histogram
}


std::istream & colvarbias_histogram::read_restart (std::istream& is)
{
  size_t const start_pos = is.tellg();

  cvm::log ("Restarting collective variable histogram \""+
            this->name+"\".\n");
  std::string key, brace, conf;

  if ( !(is >> key)   || !(key == "histogram") ||
       !(is >> brace) || !(brace == "{") ||
       !(is >> colvarparse::read_block ("configuration", conf)) ) {
    cvm::log ("Error: in reading restart configuration for histogram \""+
              this->name+"\" at position "+
              cvm::to_str (is.tellg())+" in stream.\n");
    is.clear();
    is.seekg (start_pos, std::ios::beg);
    is.setstate (std::ios::failbit);
    return is;
  }

  int id = -1;
  std::string name = "";
  if ( (colvarparse::get_keyval (conf, "name", name, std::string (""), colvarparse::parse_silent)) &&
         (name != this->name) )
    cvm::fatal_error ("Error: in the restart file, the "
                      "\"histogram\" block has a wrong name: different system?\n");
  if ( (id == -1) && (name == "") ) {
    cvm::fatal_error ("Error: \"histogram\" block in the restart file "
                      "has no name.\n");
  }

  if ( !(is >> key)   || !(key == "grid")) {
    cvm::log ("Error: in reading restart configuration for histogram \""+
              this->name+"\" at position "+
              cvm::to_str (is.tellg())+" in stream.\n");
    is.clear();
    is.seekg (start_pos, std::ios::beg);
    is.setstate (std::ios::failbit);
    return is;
  }
  if (! grid->read_raw (is)) {
    grid->read_raw_error();
  }

  is >> brace;
  if (brace != "}") {
    cvm::fatal_error ("Error: corrupt restart information for ABF bias \""+
                      this->name+"\": no matching brace at position "+
                      cvm::to_str (is.tellg())+" in the restart file.\n");
    is.setstate (std::ios::failbit);
  }
  return is;
}

std::ostream & colvarbias_histogram::write_restart (std::ostream& os)
{
  os << "histogram {\n"
     << "  configuration {\n"
     << "    name " << this->name << "\n";
  os << "  }\n";

  os << "grid\n";
  grid->write_raw (os, 8);

  os << "}\n\n";

  return os;
}
