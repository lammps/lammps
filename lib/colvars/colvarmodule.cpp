/// -*- c++ -*-

#include <sstream>
#include <string.h>

#include "colvarmodule.h"
#include "colvarparse.h"
#include "colvarproxy.h"
#include "colvar.h"
#include "colvarbias.h"
#include "colvarbias_alb.h"
#include "colvarbias_meta.h"
#include "colvarbias_abf.h"
#include "colvarbias_restraint.h"
#include "colvarscript.h"

colvarmodule::colvarmodule(colvarproxy *proxy_in)
{
  // pointer to the proxy object
  if (proxy == NULL) {
    proxy = proxy_in;
    parse = new colvarparse();
  } else {
    // TODO relax this error to handle multiple molecules in VMD
    // once the module is not static anymore
    cvm::error("Error: trying to allocate the collective "
                      "variable module twice.\n");
    return;
  }
  cvm::log(cvm::line_marker);
  cvm::log("Initializing the collective variables module, version "+
            cvm::to_str(COLVARS_VERSION)+".\n");

  // set initial default values

  // "it_restart" will be set by the input state file, if any;
  // "it" should be updated by the proxy
  colvarmodule::it = colvarmodule::it_restart = 0;
  colvarmodule::it_restart_from_state_file = true;

  colvarmodule::use_scripted_forces = false;

  colvarmodule::b_analysis = false;
  colvarmodule::debug_gradients_step_size = 1.0e-07;
  colvarmodule::rotation::crossing_threshold = 1.0e-02;

  colvarmodule::cv_traj_freq = 100;
  colvarmodule::restart_out_freq = proxy->restart_frequency();

  // by default overwrite the existing trajectory file
  colvarmodule::cv_traj_append = false;
}


int colvarmodule::config_file(char const  *config_filename)
{
  cvm::log(cvm::line_marker);
  cvm::log("Reading new configuration from file \""+
            std::string(config_filename)+"\":\n");

  // open the configfile
  config_s.open(config_filename);
  if (!config_s) {
    cvm::error("Error: in opening configuration file \""+
                      std::string(config_filename)+"\".\n",
                FILE_ERROR);
    return COLVARS_ERROR;
  }

  // read the config file into a string
  std::string conf = "";
  std::string line;
  while (colvarparse::getline_nocomments(config_s, line)) {
    conf.append(line+"\n");
  }
  config_s.close();

  return config(conf);
}


int colvarmodule::config_string(std::string const &config_str)
{
  cvm::log(cvm::line_marker);
  cvm::log("Reading new configuration:\n");
  std::istringstream config_s(config_str);

  // strip the comments away
  std::string conf = "";
  std::string line;
  while (colvarparse::getline_nocomments(config_s, line)) {
    conf.append(line+"\n");
  }
  return config(conf);
}

int colvarmodule::config(std::string &conf)
{
  int error_code = 0;

  // parse global options
  error_code |= parse_global_params(conf);

  if (error_code != COLVARS_OK || cvm::get_error()) {
    set_error_bits(INPUT_ERROR);
    return COLVARS_ERROR;
  }

  // parse the options for collective variables
  error_code |= parse_colvars(conf);

  if (error_code != COLVARS_OK || cvm::get_error()) {
    set_error_bits(INPUT_ERROR);
    return COLVARS_ERROR;
  }

  // parse the options for biases
  error_code |= parse_biases(conf);

  if (error_code != COLVARS_OK || cvm::get_error()) {
    set_error_bits(INPUT_ERROR);
    return COLVARS_ERROR;
  }

  // done parsing known keywords, check that all keywords found were valid ones
  error_code |= parse->check_keywords(conf, "colvarmodule");

  if (error_code != COLVARS_OK || cvm::get_error()) {
    set_error_bits(INPUT_ERROR);
    return COLVARS_ERROR;
  }

  cvm::log(cvm::line_marker);
  cvm::log("Collective variables module (re)initialized.\n");
  cvm::log(cvm::line_marker);

  // configuration might have changed, better redo the labels
  write_traj_label(cv_traj_os);

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


int colvarmodule::parse_global_params(std::string const &conf)
{
  std::string index_file_name;
  if (parse->get_keyval(conf, "indexFile", index_file_name)) {
    read_index_file(index_file_name.c_str());
  }

  parse->get_keyval(conf, "analysis", b_analysis, b_analysis);

  parse->get_keyval(conf, "debugGradientsStepSize", debug_gradients_step_size,
                     debug_gradients_step_size,
                     colvarparse::parse_silent);

  parse->get_keyval(conf, "eigenvalueCrossingThreshold",
                     colvarmodule::rotation::crossing_threshold,
                     colvarmodule::rotation::crossing_threshold,
                     colvarparse::parse_silent);

  parse->get_keyval(conf, "colvarsTrajFrequency", cv_traj_freq, cv_traj_freq);
  parse->get_keyval(conf, "colvarsRestartFrequency",
                     restart_out_freq, restart_out_freq);

  // if this is true when initializing, it means
  // we are continuing after a reset(): default to true
  parse->get_keyval(conf, "colvarsTrajAppend", cv_traj_append, cv_traj_append);

  parse->get_keyval(conf, "scriptedColvarForces", use_scripted_forces, false,
                     colvarparse::parse_silent);

  if (use_scripted_forces && !proxy->force_script_defined) {
    cvm::fatal_error("User script for scripted colvar forces not found.");
  }

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


int colvarmodule::parse_colvars(std::string const &conf)
{
  if (cvm::debug())
    cvm::log("Initializing the collective variables.\n");

  std::string colvar_conf = "";
  size_t pos = 0;
  while (parse->key_lookup(conf, "colvar", colvar_conf, pos)) {

    if (colvar_conf.size()) {
      cvm::log(cvm::line_marker);
      cvm::increase_depth();
      colvars.push_back(new colvar(colvar_conf));
      if (cvm::get_error() ||
          ((colvars.back())->check_keywords(colvar_conf, "colvar") != COLVARS_OK)) {
            cvm::log("Error while constructing colvar number " +
        cvm::to_str(colvars.size()) + " : deleting.");
        delete colvars.back();  // the colvar destructor updates the colvars array
        return COLVARS_ERROR;
      }
      cvm::decrease_depth();
    } else {
      cvm::error("Error: \"colvar\" keyword found without any configuration.\n", INPUT_ERROR);
      return COLVARS_ERROR;
    }
    cvm::decrease_depth();
    colvar_conf = "";
  }

  if (!colvars.size()) {
    cvm::log("Warning: no collective variables defined.\n");
  }

  if (colvars.size())
    cvm::log(cvm::line_marker);
  cvm::log("Collective variables initialized, "+
            cvm::to_str(colvars.size())+
            " in total.\n");

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}

bool colvarmodule::check_new_bias(std::string &conf, char const *key)
{
  if (cvm::get_error() ||
     (biases.back()->check_keywords(conf, key) != COLVARS_OK)) {
    cvm::log("Error while constructing bias number " +
        cvm::to_str(biases.size()) + " : deleting.\n");
    delete biases.back(); // the bias destructor updates the biases array
    return true;
  }
  return false;
}


template <class bias_type>
int colvarmodule::parse_biases_type(std::string const &conf,
                                    char const *keyword,
                                    size_t &bias_count)
{
  std::string bias_conf = "";
  size_t conf_saved_pos = 0;
  while (parse->key_lookup(conf, keyword, bias_conf, conf_saved_pos)) {
    if (bias_conf.size()) {
      cvm::log(cvm::line_marker);
      cvm::increase_depth();
      biases.push_back(new bias_type(bias_conf, keyword));
      if (cvm::check_new_bias(bias_conf, keyword)) {
        return COLVARS_ERROR;
      }
      cvm::decrease_depth();
      bias_count++;
    } else {
      cvm::error("Error: keyword \""+std::string(keyword)+"\" found without configuration.\n",
                 INPUT_ERROR);
      return COLVARS_ERROR;
    }
    bias_conf = "";
  }
  return COLVARS_OK;
}


int colvarmodule::parse_biases(std::string const &conf)
{
  if (cvm::debug())
    cvm::log("Initializing the collective variables biases.\n");

  /// initialize ABF instances
  parse_biases_type<colvarbias_abf>(conf, "abf", n_abf_biases);

  /// initialize adaptive linear biases
  parse_biases_type<colvarbias_alb>(conf, "ALB", n_rest_biases);

  /// initialize harmonic restraints
  parse_biases_type<colvarbias_restraint_harmonic>(conf, "harmonic", n_rest_biases);

  /// initialize histograms
  parse_biases_type<colvarbias_histogram>(conf, "histogram", n_histo_biases);

  /// initialize linear restraints
  parse_biases_type<colvarbias_restraint_linear>(conf, "linear", n_rest_biases);

  /// initialize metadynamics instances
  parse_biases_type<colvarbias_meta>(conf, "metadynamics", n_meta_biases);

  if (use_scripted_forces) {
    cvm::log(cvm::line_marker);
    cvm::increase_depth();
    cvm::log("User forces script will be run at each bias update.");
    cvm::decrease_depth();
  }

  if (biases.size() || use_scripted_forces) {
    cvm::log(cvm::line_marker);
    cvm::log("Collective variables biases initialized, "+
             cvm::to_str(biases.size())+" in total.\n");
  } else {
    if (!use_scripted_forces) {
      cvm::log("No collective variables biases were defined.\n");
    }
  }

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


colvarbias * colvarmodule::bias_by_name(std::string const &name) {
  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    if ((*bi)->name == name) {
      return (*bi);
    }
  }
  return NULL;
}


colvar *colvarmodule::colvar_by_name(std::string const &name) {
  for (std::vector<colvar *>::iterator cvi = colvars.begin();
       cvi != colvars.end();
       cvi++) {
    if ((*cvi)->name == name) {
      return (*cvi);
    }
  }
  return NULL;
}


int colvarmodule::change_configuration(std::string const &bias_name,
                                         std::string const &conf)
{
  // This is deprecated; supported strategy is to delete the bias
  // and parse the new config
  cvm::increase_depth();
  colvarbias *b;
  b = bias_by_name(bias_name);
  if (b == NULL) { cvm::error("Error: bias not found: " + bias_name); }
  b->change_configuration(conf);
  cvm::decrease_depth();
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


std::string colvarmodule::read_colvar(std::string const &name)
{
  cvm::increase_depth();
  colvar *c;
  std::stringstream ss;
  c = colvar_by_name(name);
  if (c == NULL) { cvm::fatal_error("Error: colvar not found: " + name); }
  ss << c->value();
  cvm::decrease_depth();
  return ss.str();
}

cvm::real colvarmodule::energy_difference(std::string const &bias_name,
                                           std::string const &conf)
{
  cvm::increase_depth();
  colvarbias *b;
  cvm::real energy_diff = 0.;
  b = bias_by_name(bias_name);
  if (b == NULL) { cvm::fatal_error("Error: bias not found: " + bias_name); }
  energy_diff = b->energy_difference(conf);
  cvm::decrease_depth();
  return energy_diff;
}

int colvarmodule::bias_current_bin(std::string const &bias_name)
{
  cvm::increase_depth();
  int ret;
  colvarbias *b = bias_by_name(bias_name);

  if (b != NULL) {
    ret = b->current_bin();
  } else {
    cvm::error("Error: bias not found.\n");
    ret = COLVARS_ERROR;
  }

  cvm::decrease_depth();
  return ret;
}

int colvarmodule::bias_bin_num(std::string const &bias_name)
{
  cvm::increase_depth();
  int ret;
  colvarbias *b = bias_by_name(bias_name);

  if (b != NULL) {
    ret = b->bin_num();
  } else {
    cvm::error("Error: bias not found.\n");
    ret = COLVARS_ERROR;
  }

  cvm::decrease_depth();
  return ret;
}

int colvarmodule::bias_bin_count(std::string const &bias_name, size_t bin_index)
{
  cvm::increase_depth();
  int ret;
  colvarbias *b = bias_by_name(bias_name);

  if (b != NULL) {
    ret = b->bin_count(bin_index);
  } else {
    cvm::error("Error: bias not found.\n");
    ret = COLVARS_ERROR;
  }

  cvm::decrease_depth();
  return ret;
}

int colvarmodule::bias_share(std::string const &bias_name)
{
  cvm::increase_depth();
  int ret;
  colvarbias *b = bias_by_name(bias_name);

  if (b != NULL) {
    b->replica_share();
    ret = COLVARS_OK;
  } else {
    cvm::error("Error: bias not found.\n");
    ret = COLVARS_ERROR;
  }

  cvm::decrease_depth();
  return ret;
}


int colvarmodule::calc() {
  cvm::real total_bias_energy = 0.0;
  cvm::real total_colvar_energy = 0.0;

  std::vector<colvar *>::iterator cvi;
  std::vector<colvarbias *>::iterator bi;

  if (cvm::debug()) {
    cvm::log(cvm::line_marker);
    cvm::log("Collective variables module, step no. "+
              cvm::to_str(cvm::step_absolute())+"\n");
  }

  // calculate collective variables and their gradients
  if (cvm::debug())
    cvm::log("Calculating collective variables.\n");
  cvm::increase_depth();
  for (cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
    (*cvi)->calc();
    if (cvm::get_error()) {
      return COLVARS_ERROR;
    }
  }
  cvm::decrease_depth();

  // update the biases and communicate their forces to the collective
  // variables
  if (cvm::debug() && biases.size())
    cvm::log("Updating collective variable biases.\n");
  cvm::increase_depth();
  for (bi = biases.begin(); bi != biases.end(); bi++) {
    total_bias_energy += (*bi)->update();
    if (cvm::get_error()) {
      return COLVARS_ERROR;
    }
  }
  cvm::decrease_depth();

  // sum the forces from all biases for each collective variable
  if (cvm::debug() && biases.size())
    cvm::log("Collecting forces from all biases.\n");
  cvm::increase_depth();
  for (cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
    (*cvi)->reset_bias_force();
  }
  for (bi = biases.begin(); bi != biases.end(); bi++) {
    (*bi)->communicate_forces();
    if (cvm::get_error()) {
      return COLVARS_ERROR;
    }
  }

  // Run user force script, if provided,
  // potentially adding scripted forces to the colvars
  if (use_scripted_forces) {
    int res;
    res = proxy->run_force_callback();
    if (res == COLVARS_NOT_IMPLEMENTED) {
      cvm::error("Colvar forces scripts are not implemented.");
      return COLVARS_ERROR;
    }
    if (res != COLVARS_OK) {
      cvm::error("Error running user colvar forces script");
      return COLVARS_ERROR;
    }
  }

  cvm::decrease_depth();

  if (cvm::b_analysis) {
    // perform runtime analysis of colvars and biases
    if (cvm::debug() && biases.size())
      cvm::log("Perform runtime analyses.\n");
    cvm::increase_depth();
    for (cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
      (*cvi)->analyse();
      if (cvm::get_error()) {
        return COLVARS_ERROR;
      }
    }
    for (bi = biases.begin(); bi != biases.end(); bi++) {
      (*bi)->analyse();
      if (cvm::get_error()) {
        return COLVARS_ERROR;
      }
    }
    cvm::decrease_depth();
  }

  // sum up the forces for each colvar, including wall forces
  // and integrate any internal
  // equation of motion (extended system)
  if (cvm::debug())
    cvm::log("Updating the internal degrees of freedom "
              "of colvars (if they have any).\n");
  cvm::increase_depth();
  for (cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
    total_colvar_energy += (*cvi)->update();
    if (cvm::get_error()) {
      return COLVARS_ERROR;
    }
  }
  cvm::decrease_depth();
  proxy->add_energy(total_bias_energy + total_colvar_energy);

  // make collective variables communicate their forces to their
  // coupled degrees of freedom (i.e. atoms)
  if (cvm::debug())
    cvm::log("Communicating forces from the colvars to the atoms.\n");
  cvm::increase_depth();
  for (cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
    if ((*cvi)->tasks[colvar::task_gradients]) {
      (*cvi)->communicate_forces();
      if (cvm::get_error()) {
        return COLVARS_ERROR;
      }
    }
  }
  cvm::decrease_depth();

  // write restart file, if needed
  if (restart_out_freq && restart_out_name.size()) {
    if ( (cvm::step_relative() > 0) &&
         ((cvm::step_absolute() % restart_out_freq) == 0) ) {
      cvm::log("Writing the state file \""+
                restart_out_name+"\".\n");
      proxy->backup_file(restart_out_name.c_str());
      restart_out_os.open(restart_out_name.c_str());
      if (!write_restart(restart_out_os))
        cvm::error("Error: in writing restart file.\n");
      restart_out_os.close();
    }
  }

  // write trajectory file, if needed
  if (cv_traj_freq && cv_traj_name.size()) {

    if (!cv_traj_os.is_open()) {
      open_traj_file(cv_traj_name);
    }

    // write labels in the traj file every 1000 lines and at first timestep
    if ((cvm::step_absolute() % (cv_traj_freq * 1000)) == 0 || cvm::step_relative() == 0) {
      write_traj_label(cv_traj_os);
    }

    if ((cvm::step_absolute() % cv_traj_freq) == 0) {
      write_traj(cv_traj_os);
    }

    if (restart_out_freq) {
      // flush the trajectory file if we are at the restart frequency
      if ( (cvm::step_relative() > 0) &&
           ((cvm::step_absolute() % restart_out_freq) == 0) ) {
        cvm::log("Synchronizing (emptying the buffer of) trajectory file \""+
                  cv_traj_name+"\".\n");
        cv_traj_os.flush();
      }
    }
  } // end if (cv_traj_freq)

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


int colvarmodule::analyze()
{
  if (cvm::debug()) {
    cvm::log("colvarmodule::analyze(), step = "+cvm::to_str(it)+".\n");
  }

  if (cvm::step_relative() == 0)
    cvm::log("Performing analysis.\n");

  // perform colvar-specific analysis
  for (std::vector<colvar *>::iterator cvi = colvars.begin();
       cvi != colvars.end();
       cvi++) {
    cvm::increase_depth();
    (*cvi)->analyse();
    cvm::decrease_depth();
  }

  // perform bias-specific analysis
  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    cvm::increase_depth();
    (*bi)->analyse();
    cvm::decrease_depth();
  }

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}

int colvarmodule::setup()
{
  // loop over all components of all colvars to reset masses of all groups
  for (std::vector<colvar *>::iterator cvi = colvars.begin();
       cvi != colvars.end();  cvi++) {
    (*cvi)->setup();
  }
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);

}

colvarmodule::~colvarmodule()
{
  reset();
  delete parse;
  proxy = NULL;
}

int colvarmodule::reset()
{
  cvm::log("Resetting the Collective Variables Module.\n");
  // Iterate backwards because we are deleting the elements as we go
  for (std::vector<colvar *>::reverse_iterator cvi = colvars.rbegin();
       cvi != colvars.rend();
       cvi++) {
    delete *cvi; // the colvar destructor updates the colvars array
  }
  colvars.clear();

  // Iterate backwards because we are deleting the elements as we go
  for (std::vector<colvarbias *>::reverse_iterator bi = biases.rbegin();
       bi != biases.rend();
       bi++) {
    delete *bi; // the bias destructor updates the biases array
  }
  biases.clear();

  index_groups.clear();
  index_group_names.clear();

  // Do not close file here, as we might not be done with it yet.
  cv_traj_os.flush();

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


int colvarmodule::setup_input()
{
  // name of input state file
  restart_in_name = proxy->input_prefix().size() ?
    std::string(proxy->input_prefix()+".colvars.state") :
    std::string("") ;

  // read the restart configuration, if available
  if (restart_in_name.size()) {
    // read the restart file
    std::ifstream input_is(restart_in_name.c_str());
    if (!input_is.good()) {
      cvm::error("Error: in opening restart file \""+
                        std::string(restart_in_name)+"\".\n",
                FILE_ERROR);
      return COLVARS_ERROR;
    } else {
      cvm::log("Restarting from file \""+restart_in_name+"\".\n");
      read_restart(input_is);
      if (cvm::get_error() != COLVARS_OK) {
        return COLVARS_ERROR;
      }
      cvm::log(cvm::line_marker);
    }
  }
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);

}


int colvarmodule::setup_output()
{
  // output state file (restart)
  restart_out_name = proxy->restart_output_prefix().size() ?
    std::string(proxy->restart_output_prefix()+".colvars.state") :
    std::string("");

  if (restart_out_name.size()) {
    cvm::log("The restart output state file will be \""+restart_out_name+"\".\n");
  }

  output_prefix = proxy->output_prefix();
  if (output_prefix.size()) {
    cvm::log("The final output state file will be \""+
              (output_prefix.size() ?
               std::string(output_prefix+".colvars.state") :
               std::string("colvars.state"))+"\".\n");
    // cvm::log (cvm::line_marker);
  }

  cv_traj_name =
    (output_prefix.size() ?
     std::string(output_prefix+".colvars.traj") :
     std::string(""));

  if (cv_traj_freq && cv_traj_name.size()) {
    open_traj_file(cv_traj_name);
  }

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


std::istream & colvarmodule::read_restart(std::istream &is)
{
  {
    // read global restart information
    std::string restart_conf;
    if (is >> colvarparse::read_block("configuration", restart_conf)) {
      if (it_restart_from_state_file) {
        parse->get_keyval(restart_conf, "step",
                           it_restart, (size_t) 0,
                           colvarparse::parse_silent);
        it = it_restart;
      }
    }
    is.clear();
  }

  // colvars restart
  cvm::increase_depth();
  for (std::vector<colvar *>::iterator cvi = colvars.begin();
       cvi != colvars.end();
       cvi++) {
    if ( !((*cvi)->read_restart(is)) ) {
      cvm::error("Error: in reading restart configuration for collective variable \""+
                        (*cvi)->name+"\".\n",
                INPUT_ERROR);
    }
  }

  // biases restart
  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    if (!((*bi)->read_restart(is)))
      cvm::error("Error: in reading restart configuration for bias \""+
                   (*bi)->name+"\".\n",
                INPUT_ERROR);
  }
  cvm::decrease_depth();

  return is;
}


int colvarmodule::backup_file(char const *filename)
{
  return proxy->backup_file(filename);
}


int colvarmodule::write_output_files()
{
  // if this is a simulation run (i.e. not a postprocessing), output data
  // must be written to be able to restart the simulation
  std::string const out_name =
    (output_prefix.size() ?
     std::string(output_prefix+".colvars.state") :
     std::string("colvars.state"));
  cvm::log("Saving collective variables state to \""+out_name+"\".\n");
  proxy->backup_file(out_name.c_str());
  std::ofstream out(out_name.c_str());
  out.setf(std::ios::scientific, std::ios::floatfield);
  this->write_restart(out);
  out.close();

  cvm::increase_depth();
  for (std::vector<colvar *>::iterator cvi = colvars.begin();
       cvi != colvars.end();
       cvi++) {
    (*cvi)->write_output_files();
  }
  cvm::decrease_depth();

  // do not close to avoid problems with multiple NAMD runs
  cv_traj_os.flush();
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}



int colvarmodule::read_traj(char const *traj_filename,
                              size_t      traj_read_begin,
                              size_t      traj_read_end)
{
  cvm::log("Opening trajectory file \""+
            std::string(traj_filename)+"\".\n");
  std::ifstream traj_is(traj_filename);

  while (true) {
    while (true) {

      std::string line("");

      do {
        if (!colvarparse::getline_nocomments(traj_is, line)) {
          cvm::log("End of file \""+std::string(traj_filename)+
                    "\" reached, or corrupted file.\n");
          traj_is.close();
          return false;
        }
      } while (line.find_first_not_of(colvarparse::white_space) == std::string::npos);

      std::istringstream is(line);

      if (!(is >> it)) return false;

      if ( (it < traj_read_begin) ) {

        if ((it % 1000) == 0)
          std::cerr << "Skipping trajectory step " << it
                    << "                    \r";

        continue;

      } else {

        if ((it % 1000) == 0)
          std::cerr << "Reading from trajectory, step = " << it
                    << "                    \r";

        if ( (traj_read_end > traj_read_begin) &&
             (it > traj_read_end) ) {
          std::cerr << "\n";
          cvm::error("Reached the end of the trajectory, "
                    "read_end = "+cvm::to_str(traj_read_end)+"\n",
                    FILE_ERROR);
          return COLVARS_ERROR;
        }

        for (std::vector<colvar *>::iterator cvi = colvars.begin();
             cvi != colvars.end();
             cvi++) {
          if (!(*cvi)->read_traj(is)) {
            cvm::error("Error: in reading colvar \""+(*cvi)->name+
                      "\" from trajectory file \""+
                      std::string(traj_filename)+"\".\n",
                      FILE_ERROR);
            return COLVARS_ERROR;
          }
        }

        break;
      }
    }
  }
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


std::ostream & colvarmodule::write_restart(std::ostream &os)
{
  os.setf(std::ios::scientific, std::ios::floatfield);
  os << "configuration {\n"
     << "  step " << std::setw(it_width)
     << it << "\n"
     << "  dt " << dt() << "\n"
     << "}\n\n";

  cvm::increase_depth();
  for (std::vector<colvar *>::iterator cvi = colvars.begin();
       cvi != colvars.end();
       cvi++) {
    (*cvi)->write_restart(os);
  }

  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    (*bi)->write_restart(os);
  }
  cvm::decrease_depth();

  return os;
}

int colvarmodule::open_traj_file(std::string const &file_name)
{
  if (cv_traj_os.is_open()) {
    return COLVARS_OK;
  }

  // (re)open trajectory file
  if (cv_traj_append) {
    cvm::log("Appending to colvar trajectory file \""+file_name+
              "\".\n");
    cv_traj_os.open(file_name.c_str(), std::ios::app);
  } else {
    cvm::log("Writing to colvar trajectory file \""+file_name+
              "\".\n");
    proxy->backup_file(file_name.c_str());
    cv_traj_os.open(file_name.c_str(), std::ios::out);
  }

  if (!cv_traj_os.is_open()) {
    cvm::error("Error: cannot write to file \""+file_name+"\".\n",
                FILE_ERROR);
  }

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}

int colvarmodule::close_traj_file()
{
  if (cv_traj_os.good()) {
    cv_traj_os.close();
  }
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}

std::ostream & colvarmodule::write_traj_label(std::ostream &os)
{
  if (!os.good()) {
    cvm::error("Cannot write to trajectory file.");
    return os;
  }
  os.setf(std::ios::scientific, std::ios::floatfield);

  os << "# " << cvm::wrap_string("step", cvm::it_width-2)
     << " ";

  cvm::increase_depth();
  for (std::vector<colvar *>::iterator cvi = colvars.begin();
       cvi != colvars.end();
       cvi++) {
    (*cvi)->write_traj_label(os);
  }
  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    (*bi)->write_traj_label(os);
  }
  os << "\n";
  if (cvm::debug())
    os.flush();
  cvm::decrease_depth();
  return os;
}

std::ostream & colvarmodule::write_traj(std::ostream &os)
{
  os.setf(std::ios::scientific, std::ios::floatfield);

  os << std::setw(cvm::it_width) << it
     << " ";

  cvm::increase_depth();
  for (std::vector<colvar *>::iterator cvi = colvars.begin();
       cvi != colvars.end();
       cvi++) {
    (*cvi)->write_traj(os);
  }
  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    (*bi)->write_traj(os);
  }
  os << "\n";
  if (cvm::debug())
    os.flush();
  cvm::decrease_depth();
  return os;
}


void cvm::log(std::string const &message)
{
  if (depth > 0)
    proxy->log((std::string(2*depth, ' '))+message);
  else
    proxy->log(message);
}

void cvm::increase_depth()
{
  depth++;
}

void cvm::decrease_depth()
{
  if (depth) depth--;
}

void cvm::error(std::string const &message, int code)
{
  set_error_bits(code);
  proxy->error(message);
}

void cvm::fatal_error(std::string const &message)
{
  set_error_bits(FATAL_ERROR);
  proxy->fatal_error(message);
}

void cvm::exit(std::string const &message)
{
  proxy->exit(message);
}


int cvm::read_index_file(char const *filename)
{
  std::ifstream is(filename, std::ios::binary);
  if (!is.good())
    cvm::error("Error: in opening index file \""+
                      std::string(filename)+"\".\n",
                      FILE_ERROR);

  while (is.good()) {
    char open, close;
    std::string group_name;
    if ( (is >> open) && (open == '[') &&
         (is >> group_name) &&
         (is >> close) && (close == ']') ) {
      for (std::list<std::string>::iterator names_i = index_group_names.begin();
           names_i != index_group_names.end();
           names_i++) {
        if (*names_i == group_name) {
          cvm::error("Error: the group name \""+group_name+
                      "\" appears in multiple index files.\n",
                      FILE_ERROR);
        }
      }
      cvm::index_group_names.push_back(group_name);
      cvm::index_groups.push_back(std::vector<int> ());
    } else {
      cvm::error("Error: in parsing index file \""+
                 std::string(filename)+"\".\n",
                 INPUT_ERROR);
    }

    int atom_number = 1;
    size_t pos = is.tellg();
    while ( (is >> atom_number) && (atom_number > 0) ) {
      (cvm::index_groups.back()).push_back(atom_number);
      pos = is.tellg();
    }
    is.clear();
    is.seekg(pos, std::ios::beg);
    std::string delim;
    if ( (is >> delim) && (delim == "[") ) {
      // new group
      is.clear();
      is.seekg(pos, std::ios::beg);
    } else {
      break;
    }
  }

  cvm::log("The following index groups were read from the index file \""+
            std::string(filename)+"\":\n");
  std::list<std::string>::iterator names_i = index_group_names.begin();
  std::list<std::vector<int> >::iterator lists_i = index_groups.begin();
  for ( ; names_i != index_group_names.end() ; names_i++, lists_i++) {
    cvm::log("  "+(*names_i)+" ("+cvm::to_str(lists_i->size())+" atoms).\n");
  }
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}

int cvm::load_atoms(char const *file_name,
                             std::vector<cvm::atom> &atoms,
                             std::string const &pdb_field,
                             double const pdb_field_value)
{
  return proxy->load_atoms(file_name, atoms, pdb_field, pdb_field_value);
}

int cvm::load_coords(char const *file_name,
                              std::vector<cvm::atom_pos> &pos,
                              const std::vector<int> &indices,
                              std::string const &pdb_field,
                              double const pdb_field_value)
{
  // Differentiate between PDB and XYZ files
  // for XYZ files, use CVM internal parser
  // otherwise call proxy function for PDB

  std::string const ext(strlen(file_name) > 4 ? (file_name + (strlen(file_name) - 4)) : file_name);
  if (colvarparse::to_lower_cppstr(ext) == std::string(".xyz")) {
    if ( pdb_field.size() > 0 ) {
      cvm::error("Error: PDB column may not be specified for XYZ coordinate file.\n", INPUT_ERROR);
      return COLVARS_ERROR;
    }
    return cvm::load_coords_xyz(file_name, pos, indices);
  } else {
    return proxy->load_coords(file_name, pos, indices, pdb_field, pdb_field_value);
  }
}


int cvm::load_coords_xyz(char const *filename,
                           std::vector<atom_pos> &pos,
                           const std::vector<int> &indices)
{
  std::ifstream xyz_is(filename);
  unsigned int natoms;
  char symbol[256];
  std::string line;

  if ( ! (xyz_is >> natoms) ) {
    cvm::error("Error: cannot parse XYZ file "
                 + std::string(filename) + ".\n", INPUT_ERROR);
  }
  // skip comment line
  std::getline(xyz_is, line);
  std::getline(xyz_is, line);
  xyz_is.width(255);
  std::vector<atom_pos>::iterator pos_i = pos.begin();

  if (pos.size() != natoms) { // Use specified indices
    int next = 0; // indices are zero-based
    std::vector<int>::const_iterator index = indices.begin();
    for ( ; pos_i != pos.end() ; pos_i++, index++) {

      while ( next < *index ) {
        std::getline(xyz_is, line);
        next++;
      }
      xyz_is >> symbol;
      xyz_is >> (*pos_i)[0] >> (*pos_i)[1] >> (*pos_i)[2];
    }
  } else {          // Use all positions
    for ( ; pos_i != pos.end() ; pos_i++) {
      xyz_is >> symbol;
      xyz_is >> (*pos_i)[0] >> (*pos_i)[1] >> (*pos_i)[2];
    }
  }
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}

// static pointers
std::vector<colvar *>     colvarmodule::colvars;
std::vector<colvarbias *> colvarmodule::biases;
size_t                    colvarmodule::n_abf_biases = 0;
size_t                    colvarmodule::n_rest_biases = 0;
size_t                    colvarmodule::n_histo_biases = 0;
size_t                    colvarmodule::n_meta_biases = 0;
colvarproxy              *colvarmodule::proxy = NULL;


// static runtime data
cvm::real colvarmodule::debug_gradients_step_size = 1.0e-03;
int       colvarmodule::errorCode = 0;
size_t    colvarmodule::it = 0;
size_t    colvarmodule::it_restart = 0;
size_t    colvarmodule::restart_out_freq = 0;
size_t    colvarmodule::cv_traj_freq = 0;
size_t    colvarmodule::depth = 0;
bool      colvarmodule::b_analysis = false;
cvm::real colvarmodule::rotation::crossing_threshold = 1.0E-04;
std::list<std::string> colvarmodule::index_group_names;
std::list<std::vector<int> > colvarmodule::index_groups;
bool     colvarmodule::use_scripted_forces = false;

// file name prefixes
std::string colvarmodule::output_prefix = "";
std::string colvarmodule::restart_in_name = "";


// i/o constants
size_t const colvarmodule::it_width = 12;
size_t const colvarmodule::cv_prec  = 14;
size_t const colvarmodule::cv_width = 21;
size_t const colvarmodule::en_prec  = 14;
size_t const colvarmodule::en_width = 21;
std::string const colvarmodule::line_marker =
  "----------------------------------------------------------------------\n";
