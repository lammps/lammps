// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <sstream>
#include <cstring>

#include "colvarmodule.h"
#include "colvarparse.h"
#include "colvarproxy.h"
#include "colvar.h"
#include "colvarbias.h"
#include "colvarbias_abf.h"
#include "colvarbias_alb.h"
#include "colvarbias_histogram.h"
#include "colvarbias_meta.h"
#include "colvarbias_restraint.h"
#include "colvarscript.h"
#include "colvaratoms.h"
#include "colvarcomp.h"


colvarmodule::colvarmodule(colvarproxy *proxy_in)
{
  depth_s = 0;
  log_level_ = 10;

  cv_traj_os = NULL;

  xyz_reader_use_count = 0;

  if (proxy == NULL) {
    proxy = proxy_in; // Pointer to the proxy object
    parse = new colvarparse(); // Parsing object for global options
    version_int = proxy->get_version_from_string(COLVARS_VERSION);
  } else {
    // TODO relax this error to handle multiple molecules in VMD
    // once the module is not static anymore
    cvm::error("Error: trying to allocate the collective "
               "variable module twice.\n", BUG_ERROR);
    return;
  }

  cvm::log(cvm::line_marker);
  cvm::log("Initializing the collective variables module, version "+
           cvm::to_str(COLVARS_VERSION)+".\n");
  cvm::log("Please cite Fiorin et al, Mol Phys 2013:\n "
           "https://doi.org/10.1080/00268976.2013.813594\n"
           "in any publication based on this calculation.\n");

  if (proxy->smp_enabled() == COLVARS_OK) {
    cvm::log("SMP parallelism is enabled; if needed, use \"smp off\" to override this.\n");
  }

#if (__cplusplus >= 201103L)
  cvm::log("This version was built with the C++11 standard or higher.");
#else
  cvm::log("This version was built without the C++11 standard: some features are disabled.\n"
    "Please see the following link for details:\n"
    "https://colvars.github.io/README-c++11.html");
#endif

  // set initial default values

  // "it_restart" will be set by the input state file, if any;
  // "it" should be updated by the proxy
  colvarmodule::it = colvarmodule::it_restart = 0;

  use_scripted_forces = false;
  scripting_after_biases = false;

  colvarmodule::debug_gradients_step_size = 1.0e-07;

  colvarmodule::rotation::monitor_crossings = false;
  colvarmodule::rotation::crossing_threshold = 1.0e-02;

  cv_traj_freq = 100;
  restart_out_freq = proxy->restart_frequency();

  // by default overwrite the existing trajectory file
  cv_traj_append = false;

  cv_traj_write_labels = true;
}


colvarmodule * colvarmodule::main()
{
  return proxy->colvars;
}


std::vector<colvar *> *colvarmodule::variables()
{
  return &colvars;
}


std::vector<colvar *> *colvarmodule::variables_active()
{
  return &colvars_active;
}


std::vector<colvar *> *colvarmodule::variables_active_smp()
{
  return &colvars_smp;
}


std::vector<int> *colvarmodule::variables_active_smp_items()
{
  return &colvars_smp_items;
}


std::vector<colvarbias *> *colvarmodule::biases_active()
{
  return &(biases_active_);
}


size_t colvarmodule::size() const
{
  return colvars.size() + biases.size();
}


int colvarmodule::read_config_file(char const  *config_filename)
{
  cvm::log(cvm::line_marker);
  cvm::log("Reading new configuration from file \""+
           std::string(config_filename)+"\":\n");

  // open the configfile
  config_s.open(config_filename);
  if (!config_s.is_open()) {
    cvm::error("Error: in opening configuration file \""+
               std::string(config_filename)+"\".\n",
               FILE_ERROR);
    return COLVARS_ERROR;
  }

  // read the config file into a string
  std::string conf = "";
  std::string line;
  while (parse->read_config_line(config_s, line)) {
    // Delete lines that contain only white space after removing comments
    if (line.find_first_not_of(colvarparse::white_space) != std::string::npos)
      conf.append(line+"\n");
  }
  config_s.close();

  return parse_config(conf);
}


int colvarmodule::read_config_string(std::string const &config_str)
{
  cvm::log(cvm::line_marker);
  cvm::log("Reading new configuration:\n");
  std::istringstream new_config_s(config_str);

  // strip the comments away
  std::string conf = "";
  std::string line;
  while (parse->read_config_line(new_config_s, line)) {
    // Delete lines that contain only white space after removing comments
    if (line.find_first_not_of(colvarparse::white_space) != std::string::npos)
      conf.append(line+"\n");
  }

  return parse_config(conf);
}


std::istream & colvarmodule::getline(std::istream &is, std::string &line)
{
  std::string l;
  if (std::getline(is, l)) {
    size_t const sz = l.size();
    if (sz > 0) {
      if (l[sz-1] == '\r' ) {
        line = l.substr(0, sz-1);
      } else {
        line = l;
      }
    } else {
      line.clear();
    }
  }
  return is;
}


int colvarmodule::parse_config(std::string &conf)
{
  extra_conf.clear();

  // Check that the input has matching braces
  if (colvarparse::check_braces(conf, 0) != COLVARS_OK) {
    return cvm::error("Error: unmatched curly braces in configuration.\n",
                      INPUT_ERROR);
  }

  // Parse global options
  if (catch_input_errors(parse_global_params(conf))) {
    return get_error();
  }

  // Parse the options for collective variables
  if (catch_input_errors(parse_colvars(conf))) {
    return get_error();
  }

  // Parse the options for biases
  if (catch_input_errors(parse_biases(conf))) {
    return get_error();
  }

  // Done parsing known keywords, check that all keywords found were valid ones
  if (catch_input_errors(parse->check_keywords(conf, "colvarmodule"))) {
    return get_error();
  }

  // Parse auto-generated configuration (e.g. for back-compatibility)
  if (extra_conf.size()) {
    catch_input_errors(parse_global_params(extra_conf));
    catch_input_errors(parse_colvars(extra_conf));
    catch_input_errors(parse_biases(extra_conf));
    parse->check_keywords(extra_conf, "colvarmodule");
    extra_conf.clear();
    if (get_error() != COLVARS_OK) return get_error();
  }

  cvm::log(cvm::line_marker);
  cvm::log("Collective variables module (re)initialized.\n");
  cvm::log(cvm::line_marker);

  // Update any necessary proxy data
  proxy->setup();

  // configuration might have changed, better redo the labels
  cv_traj_write_labels = true;

  return get_error();
}


std::string const & colvarmodule::get_config() const
{
  return parse->get_config();
}


int colvarmodule::append_new_config(std::string const &new_conf)
{
  extra_conf += new_conf;
  return COLVARS_OK;
}


int colvarmodule::parse_global_params(std::string const &conf)
{
  // TODO document and then echo this keyword
  parse->get_keyval(conf, "logLevel", log_level_, log_level_,
                    colvarparse::parse_silent);
  {
    std::string units;
    if (parse->get_keyval(conf, "units", units)) {
      units = colvarparse::to_lower_cppstr(units);
      int error_code = proxy->set_unit_system(units, (colvars.size() != 0));
      if (error_code != COLVARS_OK) {
        return error_code;
      }
    }
  }

  {
    std::string index_file_name;
    size_t pos = 0;
    while (parse->key_lookup(conf, "indexFile", &index_file_name, &pos)) {
      cvm::log("# indexFile = \""+index_file_name+"\"\n");
      read_index_file(index_file_name.c_str());
      index_file_name.clear();
    }
  }

  if (parse->get_keyval(conf, "smp", proxy->b_smp_active, proxy->b_smp_active)) {
    if (proxy->b_smp_active == false) {
      cvm::log("SMP parallelism has been disabled.\n");
    }
  }

  bool b_analysis = true;
  if (parse->get_keyval(conf, "analysis", b_analysis, true,
                        colvarparse::parse_silent)) {
    cvm::log("Warning: keyword \"analysis\" is deprecated: it is now set "
             "to true; individual analyses are performed only if requested.");
  }

  parse->get_keyval(conf, "debugGradientsStepSize", debug_gradients_step_size,
                    debug_gradients_step_size,
                    colvarparse::parse_silent);

  parse->get_keyval(conf, "monitorEigenvalueCrossing",
                    colvarmodule::rotation::monitor_crossings,
                    colvarmodule::rotation::monitor_crossings,
                    colvarparse::parse_silent);
  parse->get_keyval(conf, "eigenvalueCrossingThreshold",
                    colvarmodule::rotation::crossing_threshold,
                    colvarmodule::rotation::crossing_threshold,
                    colvarparse::parse_silent);

  parse->get_keyval(conf, "colvarsTrajFrequency", cv_traj_freq, cv_traj_freq);
  parse->get_keyval(conf, "colvarsRestartFrequency",
                    restart_out_freq, restart_out_freq);

  // Deprecate append flag
  parse->get_keyval(conf, "colvarsTrajAppend",
                    cv_traj_append, cv_traj_append, colvarparse::parse_silent);

  parse->get_keyval(conf, "scriptedColvarForces",
                    use_scripted_forces, use_scripted_forces);

  parse->get_keyval(conf, "scriptingAfterBiases",
                    scripting_after_biases, scripting_after_biases);

  if (use_scripted_forces && !proxy->force_script_defined) {
    return cvm::error("User script for scripted colvar forces not found.",
                      INPUT_ERROR);
  }

  return cvm::get_error();
}


int colvarmodule::parse_colvars(std::string const &conf)
{
  if (cvm::debug())
    cvm::log("Initializing the collective variables.\n");

  std::string colvar_conf = "";
  size_t pos = 0;
  while (parse->key_lookup(conf, "colvar", &colvar_conf, &pos)) {

    if (colvar_conf.size()) {
      cvm::log(cvm::line_marker);
      cvm::increase_depth();
      colvars.push_back(new colvar());
      if (((colvars.back())->init(colvar_conf) != COLVARS_OK) ||
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
                                    char const *keyword)
{
  std::string bias_conf = "";
  size_t conf_saved_pos = 0;
  while (parse->key_lookup(conf, keyword, &bias_conf, &conf_saved_pos)) {
    if (bias_conf.size()) {
      cvm::log(cvm::line_marker);
      cvm::increase_depth();
      biases.push_back(new bias_type(keyword));
      biases.back()->init(bias_conf);
      if (cvm::check_new_bias(bias_conf, keyword) != COLVARS_OK) {
        return COLVARS_ERROR;
      }
      cvm::decrease_depth();
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
  parse_biases_type<colvarbias_abf>(conf, "abf");

  /// initialize adaptive linear biases
  parse_biases_type<colvarbias_alb>(conf, "ALB");

  /// initialize harmonic restraints
  parse_biases_type<colvarbias_restraint_harmonic>(conf, "harmonic");

  /// initialize harmonic walls restraints
  parse_biases_type<colvarbias_restraint_harmonic_walls>(conf, "harmonicWalls");

  /// initialize histograms
  parse_biases_type<colvarbias_histogram>(conf, "histogram");

  /// initialize histogram restraints
  parse_biases_type<colvarbias_restraint_histogram>(conf, "histogramRestraint");

  /// initialize linear restraints
  parse_biases_type<colvarbias_restraint_linear>(conf, "linear");

  /// initialize metadynamics instances
  parse_biases_type<colvarbias_meta>(conf, "metadynamics");

  if (use_scripted_forces) {
    cvm::log(cvm::line_marker);
    cvm::increase_depth();
    cvm::log("User forces script will be run at each bias update.");
    cvm::decrease_depth();
  }

  std::vector<std::string> const time_biases = time_dependent_biases();
  if (time_biases.size() > 1) {
    cvm::log("WARNING: there are "+cvm::to_str(time_biases.size())+
             " time-dependent biases with non-zero force parameters:\n"+
             cvm::to_str(time_biases)+"\n"+
             "Please ensure that their forces do not counteract each other.\n");
  }

  if (num_biases() || use_scripted_forces) {
    cvm::log(cvm::line_marker);
    cvm::log("Collective variables biases initialized, "+
             cvm::to_str(num_biases())+" in total.\n");
  } else {
    if (!use_scripted_forces) {
      cvm::log("No collective variables biases were defined.\n");
    }
  }

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


size_t colvarmodule::num_variables() const
{
  return colvars.size();
}


size_t colvarmodule::num_variables_feature(int feature_id) const
{
  size_t n = 0;
  for (std::vector<colvar *>::const_iterator cvi = colvars.begin();
       cvi != colvars.end();
       cvi++) {
    if ((*cvi)->is_enabled(feature_id)) {
      n++;
    }
  }
  return n;
}


size_t colvarmodule::num_biases() const
{
  return biases.size();
}


size_t colvarmodule::num_biases_feature(int feature_id) const
{
  size_t n = 0;
  for (std::vector<colvarbias *>::const_iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    if ((*bi)->is_enabled(feature_id)) {
      n++;
    }
  }
  return n;
}


size_t colvarmodule::num_biases_type(std::string const &type) const
{
  size_t n = 0;
  for (std::vector<colvarbias *>::const_iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    if ((*bi)->bias_type == type) {
      n++;
    }
  }
  return n;
}


std::vector<std::string> const colvarmodule::time_dependent_biases() const
{
  size_t i;
  std::vector<std::string> biases_names;
  for (i = 0; i < num_biases(); i++) {
    if (biases[i]->is_enabled(colvardeps::f_cvb_apply_force) &&
        biases[i]->is_enabled(colvardeps::f_cvb_active) &&
        (biases[i]->is_enabled(colvardeps::f_cvb_history_dependent) ||
         biases[i]->is_enabled(colvardeps::f_cvb_time_dependent))) {
      biases_names.push_back(biases[i]->name);
    }
  }
  return biases_names;
}


int colvarmodule::catch_input_errors(int result)
{
  if (result != COLVARS_OK || get_error()) {
    set_error_bits(result);
    set_error_bits(INPUT_ERROR);
    parse->init();
    return get_error();
  }
  return COLVARS_OK;
}


colvarbias * colvarmodule::bias_by_name(std::string const &name)
{
  colvarmodule *cv = cvm::main();
  for (std::vector<colvarbias *>::iterator bi = cv->biases.begin();
       bi != cv->biases.end();
       bi++) {
    if ((*bi)->name == name) {
      return (*bi);
    }
  }
  return NULL;
}


colvar *colvarmodule::colvar_by_name(std::string const &name)
{
  colvarmodule *cv = cvm::main();
  for (std::vector<colvar *>::iterator cvi = cv->colvars.begin();
       cvi != cv->colvars.end();
       cvi++) {
    if ((*cvi)->name == name) {
      return (*cvi);
    }
  }
  return NULL;
}


cvm::atom_group *colvarmodule::atom_group_by_name(std::string const &name)
{
  colvarmodule *cv = cvm::main();
  for (std::vector<cvm::atom_group *>::iterator agi = cv->named_atom_groups.begin();
       agi != cv->named_atom_groups.end();
       agi++) {
    if ((*agi)->name == name) {
      return (*agi);
    }
  }
  return NULL;
}


void colvarmodule::register_named_atom_group(atom_group *ag) {
  named_atom_groups.push_back(ag);
}


void colvarmodule::unregister_named_atom_group(cvm::atom_group *ag)
{
  for (std::vector<cvm::atom_group *>::iterator agi = named_atom_groups.begin();
       agi != named_atom_groups.end();
       agi++) {
    if (*agi == ag) {
      named_atom_groups.erase(agi);
      break;
    }
  }
}


int colvarmodule::change_configuration(std::string const &bias_name,
                                       std::string const &conf)
{
  // This is deprecated; supported strategy is to delete the bias
  // and parse the new config
  cvm::increase_depth();
  colvarbias *b;
  b = bias_by_name(bias_name);
  if (b == NULL) {
    cvm::error("Error: bias not found: " + bias_name);
    return COLVARS_ERROR;
  }
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
  if (c == NULL) {
    cvm::error("Error: colvar not found: " + name);
    return std::string();
  }
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
  if (b == NULL) {
    cvm::error("Error: bias not found: " + bias_name);
    return 0.;
  }
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


int colvarmodule::calc()
{
  int error_code = COLVARS_OK;

  if (cvm::debug()) {
    cvm::log(cvm::line_marker);
    cvm::log("Collective variables module, step no. "+
             cvm::to_str(cvm::step_absolute())+"\n");
  }

  error_code |= calc_colvars();
  // set biasing forces to zero before biases are calculated and summed over
  for (std::vector<colvar *>::iterator cvi = colvars.begin();
       cvi != colvars.end(); cvi++) {
    (*cvi)->reset_bias_force();
  }
  error_code |= calc_biases();
  error_code |= update_colvar_forces();

  error_code |= analyze();

  // write trajectory files, if needed
  if (cv_traj_freq && cv_traj_name.size()) {
    error_code |= write_traj_files();
  }

  // write restart files, if needed
  if (restart_out_freq && (cvm::step_relative() > 0) &&
      ((cvm::step_absolute() % restart_out_freq) == 0) ) {
    if (restart_out_name.size()) {
      // Write restart file, if different from main output
      error_code |= write_restart_file(restart_out_name);
    } else {
      error_code |= write_restart_file(output_prefix()+".colvars.state");
    }
    write_output_files();
  }

  error_code |= end_of_step();

  return error_code;
}


int colvarmodule::calc_colvars()
{
  if (cvm::debug())
    cvm::log("Calculating collective variables.\n");
  // calculate collective variables and their gradients

  // First, we need to decide which biases are awake
  // so they can activate colvars as needed
  std::vector<colvarbias *>::iterator bi;
  for (bi = biases.begin(); bi != biases.end(); bi++) {
    int tsf = (*bi)->get_time_step_factor();
    if (tsf > 0 && (step_absolute() % tsf == 0)) {
      (*bi)->enable(colvardeps::f_cvb_awake);
    } else {
      (*bi)->disable(colvardeps::f_cvb_awake);
    }
  }

  int error_code = COLVARS_OK;
  std::vector<colvar *>::iterator cvi;

  // Determine which colvars are active at this iteration
  variables_active()->clear();
  variables_active()->reserve(variables()->size());
  for (cvi = variables()->begin(); cvi != variables()->end(); cvi++) {
    // Wake up or put to sleep variables
    int tsf = (*cvi)->get_time_step_factor();
    if (tsf > 0 && (step_absolute() % tsf == 0)) {
      (*cvi)->enable(colvardeps::f_cv_awake);
    } else {
      (*cvi)->disable(colvardeps::f_cv_awake);
    }

    if ((*cvi)->is_enabled()) {
      variables_active()->push_back(*cvi);
    }
  }

  // if SMP support is available, split up the work
  if (proxy->smp_enabled() == COLVARS_OK) {

    // first, calculate how much work (currently, how many active CVCs) each colvar has

    variables_active_smp()->clear();
    variables_active_smp_items()->clear();

    variables_active_smp()->reserve(variables_active()->size());
    variables_active_smp_items()->reserve(variables_active()->size());

    // set up a vector containing all components
    cvm::increase_depth();
    for (cvi = variables_active()->begin(); cvi != variables_active()->end(); cvi++) {

      error_code |= (*cvi)->update_cvc_flags();

      size_t num_items = (*cvi)->num_active_cvcs();
      variables_active_smp()->reserve(variables_active_smp()->size() + num_items);
      variables_active_smp_items()->reserve(variables_active_smp_items()->size() + num_items);
      for (size_t icvc = 0; icvc < num_items; icvc++) {
        variables_active_smp()->push_back(*cvi);
        variables_active_smp_items()->push_back(icvc);
      }
    }
    cvm::decrease_depth();

    // calculate colvar components in parallel
    error_code |= proxy->smp_colvars_loop();

    cvm::increase_depth();
    for (cvi = variables_active()->begin(); cvi != variables_active()->end(); cvi++) {
      error_code |= (*cvi)->collect_cvc_data();
    }
    cvm::decrease_depth();

  } else {

    // calculate colvars one at a time
    cvm::increase_depth();
    for (cvi = variables_active()->begin(); cvi != variables_active()->end(); cvi++) {
      error_code |= (*cvi)->calc();
      if (cvm::get_error()) {
        return COLVARS_ERROR;
      }
    }
    cvm::decrease_depth();
  }

  error_code |= cvm::get_error();
  return error_code;
}


int colvarmodule::calc_biases()
{
  // update the biases and communicate their forces to the collective
  // variables
  if (cvm::debug() && num_biases())
    cvm::log("Updating collective variable biases.\n");

  std::vector<colvarbias *>::iterator bi;
  int error_code = COLVARS_OK;

  // Total bias energy is reset before calling scripted biases
  total_bias_energy = 0.0;

  // update the list of active biases
  // which may have changed based on f_cvb_awake in calc_colvars()
  biases_active()->clear();
  biases_active()->reserve(biases.size());
  for (bi = biases.begin(); bi != biases.end(); bi++) {
    if ((*bi)->is_enabled()) {
      biases_active()->push_back(*bi);
    }
  }

  // if SMP support is available, split up the work
  if (proxy->smp_enabled() == COLVARS_OK) {

    if (use_scripted_forces && !scripting_after_biases) {
      // calculate biases and scripted forces in parallel
      error_code |= proxy->smp_biases_script_loop();
    } else {
      // calculate biases in parallel
      error_code |= proxy->smp_biases_loop();
    }

  } else {

    if (use_scripted_forces && !scripting_after_biases) {
      error_code |= calc_scripted_forces();
    }

    cvm::increase_depth();
    for (bi = biases_active()->begin(); bi != biases_active()->end(); bi++) {
      error_code |= (*bi)->update();
      if (cvm::get_error()) {
        return error_code;
      }
    }
    cvm::decrease_depth();
  }

  for (bi = biases_active()->begin(); bi != biases_active()->end(); bi++) {
    total_bias_energy += (*bi)->get_energy();
  }

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


int colvarmodule::update_colvar_forces()
{
  int error_code = COLVARS_OK;

  std::vector<colvar *>::iterator cvi;
  std::vector<colvarbias *>::iterator bi;

  // sum the forces from all biases for each collective variable
  if (cvm::debug() && num_biases())
    cvm::log("Collecting forces from all biases.\n");
  cvm::increase_depth();
  for (bi = biases_active()->begin(); bi != biases_active()->end(); bi++) {
    (*bi)->communicate_forces();
    if (cvm::get_error()) {
      return COLVARS_ERROR;
    }
  }
  cvm::decrease_depth();

  if (use_scripted_forces && scripting_after_biases) {
    error_code |= calc_scripted_forces();
  }

  // Now we have collected energies from both built-in and scripted biases
  if (cvm::debug())
    cvm::log("Adding total bias energy: " + cvm::to_str(total_bias_energy) + "\n");
  proxy->add_energy(total_bias_energy);

  cvm::real total_colvar_energy = 0.0;
  // sum up the forces for each colvar, including wall forces
  // and integrate any internal
  // equation of motion (extended system)
  if (cvm::debug())
    cvm::log("Updating the internal degrees of freedom "
             "of colvars (if they have any).\n");
  cvm::increase_depth();
  for (cvi = variables()->begin(); cvi != variables()->end(); cvi++) {
    // Inactive colvars will only reset their forces and return 0 energy
    total_colvar_energy += (*cvi)->update_forces_energy();
    if (cvm::get_error()) {
      return COLVARS_ERROR;
    }
  }
  cvm::decrease_depth();
  if (cvm::debug())
    cvm::log("Adding total colvar energy: " + cvm::to_str(total_colvar_energy) + "\n");
  proxy->add_energy(total_colvar_energy);

  // make collective variables communicate their forces to their
  // coupled degrees of freedom (i.e. atoms)
  if (cvm::debug())
    cvm::log("Communicating forces from the colvars to the atoms.\n");
  cvm::increase_depth();
  for (cvi = variables_active()->begin(); cvi != variables_active()->end(); cvi++) {
    if ((*cvi)->is_enabled(colvardeps::f_cv_gradient)) {
      (*cvi)->communicate_forces();
      if (cvm::get_error()) {
        return COLVARS_ERROR;
      }
    }
  }
  cvm::decrease_depth();

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


int colvarmodule::calc_scripted_forces()
{
  // Run user force script, if provided,
  // potentially adding scripted forces to the colvars
  int res;
  res = proxy->run_force_callback();
  if (res == COLVARS_NOT_IMPLEMENTED) {
    cvm::error("Colvar forces scripts are not implemented.");
    return COLVARS_NOT_IMPLEMENTED;
  }
  if (res != COLVARS_OK) {
    cvm::error("Error running user colvar forces script");
    return COLVARS_ERROR;
  }
  return COLVARS_OK;
}


int colvarmodule::write_restart_file(std::string const &out_name)
{
  cvm::log("Saving collective variables state to \""+out_name+"\".\n");
  proxy->backup_file(out_name);
  std::ostream *restart_out_os = proxy->output_stream(out_name);
  if (!restart_out_os) return cvm::get_error();
  if (!write_restart(*restart_out_os)) {
    return cvm::error("Error: in writing restart file.\n", FILE_ERROR);
  }
  proxy->close_output_stream(out_name);
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


int colvarmodule::write_traj_files()
{
  if (cv_traj_os == NULL) {
    if (open_traj_file(cv_traj_name) != COLVARS_OK) {
      return cvm::get_error();
    } else {
      cv_traj_write_labels = true;
    }
  }

  // write labels in the traj file every 1000 lines and at first timestep
  if ((cvm::step_absolute() % (cv_traj_freq * 1000)) == 0 ||
      cvm::step_relative() == 0 ||
      cv_traj_write_labels) {
    write_traj_label(*cv_traj_os);
  }
  cv_traj_write_labels = false;

  if ((cvm::step_absolute() % cv_traj_freq) == 0) {
    write_traj(*cv_traj_os);
  }

  if (restart_out_freq && (cv_traj_os != NULL)) {
    // flush the trajectory file if we are at the restart frequency
    if ( (cvm::step_relative() > 0) &&
         ((cvm::step_absolute() % restart_out_freq) == 0) ) {
      cvm::log("Synchronizing (emptying the buffer of) trajectory file \""+
               cv_traj_name+"\".\n");
      proxy->flush_output_stream(cv_traj_os);
    }
  }

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


int colvarmodule::analyze()
{
  if (cvm::debug()) {
    cvm::log("colvarmodule::analyze(), step = "+cvm::to_str(it)+".\n");
  }

  // perform colvar-specific analysis
  for (std::vector<colvar *>::iterator cvi = variables_active()->begin();
       cvi != variables_active()->end();
       cvi++) {
    cvm::increase_depth();
    (*cvi)->analyze();
    cvm::decrease_depth();
  }

  // perform bias-specific analysis
  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    cvm::increase_depth();
    (*bi)->analyze();
    cvm::decrease_depth();
  }

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


int colvarmodule::end_of_step()
{
  if (cvm::debug()) {
    cvm::log("colvarmodule::end_of_step(), step = "+cvm::to_str(it)+".\n");
  }

  for (std::vector<colvar *>::iterator cvi = variables_active()->begin();
       cvi != variables_active()->end();
       cvi++) {
    cvm::increase_depth();
    (*cvi)->end_of_step();
    cvm::decrease_depth();
  }

  // perform bias-specific analysis
  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    cvm::increase_depth();
    (*bi)->end_of_step();
    cvm::decrease_depth();
  }

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


int colvarmodule::setup()
{
  if (this->size() == 0) return cvm::get_error();
  for (std::vector<colvar *>::iterator cvi = variables()->begin();
       cvi != variables()->end();  cvi++) {
    (*cvi)->setup();
  }
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


colvarmodule::~colvarmodule()
{
  if ((proxy->smp_thread_id() == COLVARS_NOT_IMPLEMENTED) ||
      (proxy->smp_thread_id() == 0)) {

    reset();

    // Delete contents of static arrays
    colvarbias::delete_features();
    colvar::delete_features();
    colvar::cvc::delete_features();
    atom_group::delete_features();

    delete parse;
    parse = NULL;
    proxy = NULL;
  }
}


int colvarmodule::reset()
{
  cvm::log("Resetting the Collective Variables module.\n");

  parse->init();

  // Iterate backwards because we are deleting the elements as we go
  for (std::vector<colvarbias *>::reverse_iterator bi = biases.rbegin();
       bi != biases.rend();
       bi++) {
    delete *bi; // the bias destructor updates the biases array
  }
  biases.clear();
  biases_active_.clear();

  // Iterate backwards because we are deleting the elements as we go
  for (std::vector<colvar *>::reverse_iterator cvi = colvars.rbegin();
       cvi != colvars.rend();
       cvi++) {
    delete *cvi; // the colvar destructor updates the colvars array
  }
  colvars.clear();

  reset_index_groups();

  proxy->reset();

  if (cv_traj_os != NULL) {
    // Do not close traj file here, as we might not be done with it yet.
    proxy->flush_output_stream(cv_traj_os);
  }

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}


int colvarmodule::setup_input()
{
  std::string restart_in_name("");

  // read the restart configuration, if available
  if (proxy->input_prefix().size()) {
    // read the restart file
    restart_in_name = proxy->input_prefix();
    std::ifstream input_is(restart_in_name.c_str());
    if (!input_is.good()) {
      // try by adding the suffix
      input_is.clear();
      restart_in_name = restart_in_name+std::string(".colvars.state");
      input_is.open(restart_in_name.c_str());
    }

    if (!input_is.good()) {
      cvm::error("Error: in opening input file \""+
                 std::string(restart_in_name)+"\".\n",
                 FILE_ERROR);
      return COLVARS_ERROR;
    } else {
      cvm::log(cvm::line_marker);
      cvm::log("Restarting from file \""+restart_in_name+"\".\n");
      read_restart(input_is);
      if (cvm::get_error() != COLVARS_OK) {
        return COLVARS_ERROR;
      } else {
        proxy->input_prefix().clear();
      }
      cvm::log(cvm::line_marker);
    }
  }

  return cvm::get_error();
}


int colvarmodule::setup_output()
{
  int error_code = COLVARS_OK;

  // output state file (restart)
  restart_out_name = proxy->restart_output_prefix().size() ?
    std::string(proxy->restart_output_prefix()+".colvars.state") :
    std::string("");

  if (restart_out_name.size()) {
    cvm::log("The restart output state file will be \""+
             restart_out_name+"\".\n");
  }

  output_prefix() = proxy->output_prefix();
  if (output_prefix().size()) {
    cvm::log("The final output state file will be \""+
             (output_prefix().size() ?
              std::string(output_prefix()+".colvars.state") :
              std::string("colvars.state"))+"\".\n");
    // cvm::log (cvm::line_marker);
  }

  cv_traj_name =
    (output_prefix().size() ?
     std::string(output_prefix()+".colvars.traj") :
     std::string(""));

  if (cv_traj_freq && cv_traj_name.size()) {
    error_code |= open_traj_file(cv_traj_name);
  }

  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    error_code |= (*bi)->setup_output();
  }

  if (error_code != COLVARS_OK || cvm::get_error()) {
    set_error_bits(FILE_ERROR);
  }

  return cvm::get_error();
}


std::istream & colvarmodule::read_restart(std::istream &is)
{
  bool warn_total_forces = false;

  {
    // read global restart information
    std::string restart_conf;
    if (is >> colvarparse::read_block("configuration", &restart_conf)) {

      parse->get_keyval(restart_conf, "step",
                        it_restart, static_cast<step_number>(0),
                        colvarparse::parse_restart);
      it = it_restart;

      std::string restart_version;
      int restart_version_int = 0;
      parse->get_keyval(restart_conf, "version",
                        restart_version, std::string(""),
                        colvarparse::parse_restart);
      if (restart_version.size()) {
        if (restart_version != std::string(COLVARS_VERSION)) {
          cvm::log("This state file was generated with version "+
                   restart_version+"\n");
        }
        restart_version_int =
          proxy->get_version_from_string(restart_version.c_str());
      }

      if (restart_version_int < 20160810) {
        // check for total force change
        if (proxy->total_forces_enabled()) {
          warn_total_forces = true;
        }
      }

      std::string units_restart;
      if (parse->get_keyval(restart_conf, "units",
                            units_restart, std::string(""),
                            colvarparse::parse_restart)) {
        units_restart = colvarparse::to_lower_cppstr(units_restart);
        if ((proxy->units.size() > 0) && (units_restart != proxy->units)) {
          cvm::error("Error: the state file has units \""+units_restart+
                     "\", but the current unit system is \""+proxy->units+
                     "\".\n", INPUT_ERROR);
        }
      }

    }
    is.clear();
    parse->clear_keyword_registry();
  }

  print_total_forces_errning(warn_total_forces);

  read_objects_state(is);

  return is;
}



std::istream & colvarmodule::read_objects_state(std::istream &is)
{
  size_t pos = 0;
  std::string word;

  while (is.good()) {
    pos = is.tellg();
    word.clear();
    is >> word;

    if (word.size()) {

      is.seekg(pos, std::ios::beg);

      if (word == "colvar") {

        cvm::increase_depth();
        for (std::vector<colvar *>::iterator cvi = colvars.begin();
             cvi != colvars.end();
             cvi++) {
          if ( !((*cvi)->read_state(is)) ) {
            // Here an error signals that the variable is a match, but the
            // state is corrupt; otherwise, the variable rewinds is silently
            cvm::error("Error: in reading restart configuration for "
                       "collective variable \""+(*cvi)->name+"\".\n",
                       INPUT_ERROR);
          }
          if (static_cast<size_t>(is.tellg()) > pos) break; // found it
        }
        cvm::decrease_depth();

      } else {

        cvm::increase_depth();
        for (std::vector<colvarbias *>::iterator bi = biases.begin();
             bi != biases.end();
             bi++) {
          if (((*bi)->state_keyword != word) && (*bi)->bias_type != word) {
            // Skip biases with different type; state_keyword is used to
            // support different versions of the state file format
            continue;
          }
          if (!((*bi)->read_state(is))) {
            // Same as above, an error means a match but the state is incorrect
            cvm::error("Error: in reading restart configuration for bias \""+
                       (*bi)->name+"\".\n",
                       INPUT_ERROR);
          }
          if (static_cast<size_t>(is.tellg()) > pos) break; // found it
        }
        cvm::decrease_depth();
      }
    }

    if (static_cast<size_t>(is.tellg()) == pos) {
      // This block has not been read by any object: discard it and move on
      // to the next one
      is >> colvarparse::read_block(word, NULL);
    }

    if (!is) break;
  }

  return is;
}


int colvarmodule::print_total_forces_errning(bool warn_total_forces)
{
  if (warn_total_forces) {
    cvm::log(cvm::line_marker);
    cvm::log("WARNING: The definition of system forces has changed.  Please see:\n");
    cvm::log("  https://colvars.github.io/README-totalforce.html\n");
    // update this ahead of time in this special case
    output_prefix() = proxy->input_prefix();
    cvm::log("All output files will now be saved with the prefix \""+output_prefix()+".tmp.*\".\n");
    cvm::log("Please review the important warning above. After that, you may rename:\n\
\""+output_prefix()+".tmp.colvars.state\"\n\
to:\n\
\""+proxy->input_prefix()+".colvars.state\"\n\
and load it to continue this simulation.\n");
    output_prefix() = output_prefix()+".tmp";
    write_restart_file(output_prefix()+".colvars.state");
    return cvm::error("Exiting with error until issue is addressed.\n",
                      INPUT_ERROR);
  }

  return COLVARS_OK;
}


int colvarmodule::backup_file(char const *filename)
{
  return proxy->backup_file(filename);
}


int colvarmodule::write_output_files()
{
  int error_code = COLVARS_OK;

  cvm::increase_depth();
  for (std::vector<colvar *>::iterator cvi = colvars.begin();
       cvi != colvars.end();
       cvi++) {
    error_code |= (*cvi)->write_output_files();
  }
  cvm::decrease_depth();

  cvm::increase_depth();
  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    error_code |= (*bi)->write_output_files();
    error_code |= (*bi)->write_state_to_replicas();
  }
  cvm::decrease_depth();

  if (cv_traj_os != NULL) {
    // do not close, there may be another run command
    proxy->flush_output_stream(cv_traj_os);
  }

  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}



int colvarmodule::read_traj(char const *traj_filename,
                            long        traj_read_begin,
                            long        traj_read_end)
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
     << "  version " << std::string(COLVARS_VERSION) << "\n";
  if (proxy->units.size() > 0) {
    os << "  units " << proxy->units << "\n";
  }
  os << "}\n\n";

  int error_code = COLVARS_OK;

  cvm::increase_depth();
  for (std::vector<colvar *>::iterator cvi = colvars.begin();
       cvi != colvars.end();
       cvi++) {
    (*cvi)->write_state(os);
  }

  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    (*bi)->write_state(os);
  }
  cvm::decrease_depth();

  if (error_code != COLVARS_OK) {
    // TODO make this function return an int instead
    os.setstate(std::ios::failbit);
  }

  return os;
}


int colvarmodule::open_traj_file(std::string const &file_name)
{
  if (cv_traj_os != NULL) {
    return COLVARS_OK;
  }

  // (re)open trajectory file
  if (cv_traj_append) {
    cvm::log("Appending to colvar trajectory file \""+file_name+
             "\".\n");
    cv_traj_os = (cvm::proxy)->output_stream(file_name, std::ios::app);
  } else {
    cvm::log("Writing to colvar trajectory file \""+file_name+
             "\".\n");
    proxy->backup_file(file_name.c_str());
    cv_traj_os = (cvm::proxy)->output_stream(file_name);
  }

  if (cv_traj_os == NULL) {
    cvm::error("Error: cannot write to file \""+file_name+"\".\n",
               FILE_ERROR);
  }

  return cvm::get_error();
}


int colvarmodule::close_traj_file()
{
  if (cv_traj_os != NULL) {
    proxy->close_output_stream(cv_traj_name);
    cv_traj_os = NULL;
  }
  return cvm::get_error();
}


std::ostream & colvarmodule::write_traj_label(std::ostream &os)
{
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

  if (cvm::debug()) {
    proxy->flush_output_stream(&os);
  }

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

  if (cvm::debug()) {
    proxy->flush_output_stream(&os);
  }

  cvm::decrease_depth();
  return os;
}


void cvm::log(std::string const &message, int min_log_level)
{
  if (cvm::log_level() < min_log_level) return;
  // allow logging when the module is not fully initialized
  size_t const d = (cvm::main() != NULL) ? depth() : 0;
  if (d > 0) {
    proxy->log((std::string(2*d, ' '))+message);
  } else {
    proxy->log(message);
  }
}


void cvm::increase_depth()
{
  (depth())++;
}


void cvm::decrease_depth()
{
  if (depth() > 0) {
    (depth())--;
  }
}


size_t & cvm::depth()
{
  // NOTE: do not call log() or error() here, to avoid recursion
  colvarmodule *cv = cvm::main();
  if (proxy->smp_enabled() == COLVARS_OK) {
    int const nt = proxy->smp_num_threads();
    if (int(cv->depth_v.size()) != nt) {
      proxy->smp_lock();
      // update array of depths
      if (cv->depth_v.size() > 0) { cv->depth_s = cv->depth_v[0]; }
      cv->depth_v.clear();
      cv->depth_v.assign(nt, cv->depth_s);
      proxy->smp_unlock();
    }
    return cv->depth_v[proxy->smp_thread_id()];
  }
  return cv->depth_s;
}


void colvarmodule::set_error_bits(int code)
{
  if (code < 0) {
    cvm::fatal_error("Error: set_error_bits() received negative error code.\n");
    return;
  }
  proxy->smp_lock();
  errorCode |= code | COLVARS_ERROR;
  proxy->smp_unlock();
}

bool colvarmodule::get_error_bit(int code)
{
  return bool(errorCode & code);
}

void colvarmodule::clear_error()
{
  proxy->smp_lock();
  errorCode = COLVARS_OK;
  proxy->smp_unlock();
}


int colvarmodule::error(std::string const &message, int code)
{
  set_error_bits(code);
  proxy->error(message);
  return get_error();
}


int colvarmodule::fatal_error(std::string const &message)
{
  set_error_bits(FATAL_ERROR);
  proxy->fatal_error(message);
  return get_error();
}


int cvm::read_index_file(char const *filename)
{
  std::ifstream is(filename, std::ios::binary);
  if (!is.good()) {
    cvm::error("Error: in opening index file \""+
               std::string(filename)+"\".\n",
               FILE_ERROR);
  }

  while (is.good()) {
    char open, close;
    std::string group_name;
    int index_of_group = -1;
    if ( (is >> open) && (open == '[') &&
         (is >> group_name) &&
         (is >> close) && (close == ']') ) {
      size_t i = 0;
      for ( ; i < index_group_names.size(); i++) {
        if (index_group_names[i] == group_name) {
          // Found a group with the same name
          index_of_group = i;
        }
      }
      if (index_of_group < 0) {
        index_group_names.push_back(group_name);
        index_groups.push_back(NULL);
        index_of_group = index_groups.size()-1;
      }
    } else {
      return cvm::error("Error: in parsing index file \""+
                        std::string(filename)+"\".\n",
                        INPUT_ERROR);
    }

    std::vector<int> *old_index_group = index_groups[index_of_group];
    std::vector<int> *new_index_group = new std::vector<int>();

    int atom_number = 1;
    size_t pos = is.tellg();
    while ( (is >> atom_number) && (atom_number > 0) ) {
      new_index_group->push_back(atom_number);
      pos = is.tellg();
    }

    if (old_index_group != NULL) {
      bool equal = false;
      if (new_index_group->size() == old_index_group->size()) {
        if (std::equal(new_index_group->begin(), new_index_group->end(),
                       old_index_group->begin())) {
          equal = true;
        }
      }
      if (! equal) {
        new_index_group->clear();
        delete new_index_group;
        new_index_group = NULL;
        return cvm::error("Error: the index group \""+group_name+
                          "\" was redefined.\n", INPUT_ERROR);
      } else {
        old_index_group->clear();
        delete old_index_group;
        old_index_group = NULL;
      }
    }

    index_groups[index_of_group] = new_index_group;

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

  cvm::log("The following index groups are currently defined:\n");
  size_t i = 0;
  for ( ; i < index_group_names.size(); i++) {
    cvm::log("  "+(index_group_names[i])+" ("+
             cvm::to_str((index_groups[i])->size())+" atoms)\n");
  }

  return COLVARS_OK;
}


int colvarmodule::reset_index_groups()
{
  size_t i = 0;
  for ( ; i < index_groups.size(); i++) {
    delete index_groups[i];
    index_groups[i] = NULL;
  }
  index_group_names.clear();
  index_groups.clear();
  return COLVARS_OK;
}


int cvm::load_atoms(char const *file_name,
                    cvm::atom_group &atoms,
                    std::string const &pdb_field,
                    double pdb_field_value)
{
  return proxy->load_atoms(file_name, atoms, pdb_field, pdb_field_value);
}


int cvm::load_coords(char const *file_name,
                     std::vector<cvm::rvector> *pos,
                     cvm::atom_group *atoms,
                     std::string const &pdb_field,
                     double pdb_field_value)
{
  int error_code = COLVARS_OK;

  std::string const ext(strlen(file_name) > 4 ?
                        (file_name + (strlen(file_name) - 4)) :
                        file_name);

  atoms->create_sorted_ids();

  std::vector<cvm::rvector> sorted_pos(atoms->size(), cvm::rvector(0.0));

  // Differentiate between PDB and XYZ files
  if (colvarparse::to_lower_cppstr(ext) == std::string(".xyz")) {
    if (pdb_field.size() > 0) {
      return cvm::error("Error: PDB column may not be specified "
                        "for XYZ coordinate files.\n", INPUT_ERROR);
    }
    // For XYZ files, use internal parser
    error_code |= cvm::main()->load_coords_xyz(file_name, &sorted_pos, atoms);
  } else {
    // Otherwise, call proxy function for PDB
    error_code |= proxy->load_coords(file_name,
                                     sorted_pos, atoms->sorted_ids(),
                                     pdb_field, pdb_field_value);
  }

  std::vector<int> const &map = atoms->sorted_ids_map();
  for (size_t i = 0; i < atoms->size(); i++) {
    (*pos)[map[i]] = sorted_pos[i];
  }

  return error_code;
}


int cvm::load_coords_xyz(char const *filename,
                         std::vector<rvector> *pos,
                         cvm::atom_group *atoms)
{
  std::ifstream xyz_is(filename);
  unsigned int natoms;
  char symbol[256];
  std::string line;
  cvm::real x = 0.0, y = 0.0, z = 0.0;

  if ( ! (xyz_is >> natoms) ) {
    cvm::error("Error: cannot parse XYZ file "
               + std::string(filename) + ".\n", INPUT_ERROR);
  }

  ++xyz_reader_use_count;
  if (xyz_reader_use_count < 2) {
    cvm::log("Warning: beginning from 2019-11-26 the XYZ file reader assumes Angstrom units.");
  }

  // skip comment line
  cvm::getline(xyz_is, line);
  cvm::getline(xyz_is, line);
  xyz_is.width(255);
  std::vector<atom_pos>::iterator pos_i = pos->begin();

  if (pos->size() != natoms) { // Use specified indices
    int next = 0; // indices are zero-based
    std::vector<int>::const_iterator index = atoms->sorted_ids().begin();
    for ( ; pos_i != pos->end() ; pos_i++, index++) {

      while ( next < *index ) {
        cvm::getline(xyz_is, line);
        next++;
      }
      xyz_is >> symbol;
      xyz_is >> x >> y >> z;
      // XYZ files are assumed to be in Angstrom (as eg. VMD will)
      (*pos_i)[0] = proxy->angstrom_to_internal(x);
      (*pos_i)[1] = proxy->angstrom_to_internal(y);
      (*pos_i)[2] = proxy->angstrom_to_internal(z);
    }
  } else {          // Use all positions
    for ( ; pos_i != pos->end() ; pos_i++) {
      xyz_is >> symbol;
      xyz_is >> x >> y >> z;
      (*pos_i)[0] = proxy->angstrom_to_internal(x);
      (*pos_i)[1] = proxy->angstrom_to_internal(y);
      (*pos_i)[2] = proxy->angstrom_to_internal(z);
    }
  }
  return (cvm::get_error() ? COLVARS_ERROR : COLVARS_OK);
}



// Wrappers to proxy functions: these may go in the future


cvm::real cvm::boltzmann()
{
  return proxy->boltzmann();
}


cvm::real cvm::temperature()
{
  return proxy->temperature();
}


cvm::real cvm::dt()
{
  return proxy->dt();
}


void cvm::request_total_force()
{
  proxy->request_total_force(true);
}


cvm::rvector cvm::position_distance(cvm::atom_pos const &pos1,
                                    cvm::atom_pos const &pos2)
{
  return proxy->position_distance(pos1, pos2);
}


cvm::real cvm::rand_gaussian(void)
{
  return proxy->rand_gaussian();
}


template<typename T> std::string _to_str(T const &x,
                                         size_t width, size_t prec)
{
  std::ostringstream os;
  if (width) os.width(width);
  if (prec) {
    os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(prec);
  }
  os << x;
  return os.str();
}


template<typename T> std::string _to_str_vector(std::vector<T> const &x,
                                                size_t width, size_t prec)
{
  if (!x.size()) return std::string("");
  std::ostringstream os;
  if (prec) {
    os.setf(std::ios::scientific, std::ios::floatfield);
  }
  os << "{ ";
  if (width) os.width(width);
  if (prec) os.precision(prec);
  os << x[0];
  for (size_t i = 1; i < x.size(); i++) {
    os << ", ";
    if (width) os.width(width);
    if (prec) os.precision(prec);
    os << x[i];
  }
  os << " }";
  return os.str();
}



std::string colvarmodule::to_str(std::string const &x)
{
  return std::string("\"")+x+std::string("\"");
}

std::string colvarmodule::to_str(char const *x)
{
  return std::string("\"")+std::string(x)+std::string("\"");
}

std::string colvarmodule::to_str(bool x)
{
  return (x ? "on" : "off");
}

std::string colvarmodule::to_str(int const &x,
                                 size_t width, size_t prec)
{
  return _to_str<int>(x, width, prec);
}

std::string colvarmodule::to_str(size_t const &x,
                                 size_t width, size_t prec)
{
  return _to_str<size_t>(x, width, prec);
}

std::string colvarmodule::to_str(long int const &x,
                                 size_t width, size_t prec)
{
  return _to_str<long int>(x, width, prec);
}

std::string colvarmodule::to_str(step_number const &x,
                                 size_t width, size_t prec)
{
  return _to_str<step_number>(x, width, prec);
}

std::string colvarmodule::to_str(cvm::real const &x,
                                 size_t width, size_t prec)
{
  return _to_str<cvm::real>(x, width, prec);
}

std::string colvarmodule::to_str(cvm::rvector const &x,
                                 size_t width, size_t prec)
{
  return _to_str<cvm::rvector>(x, width, prec);
}

std::string colvarmodule::to_str(cvm::quaternion const &x,
                                 size_t width, size_t prec)
{
  return _to_str<cvm::quaternion>(x, width, prec);
}

std::string colvarmodule::to_str(colvarvalue const &x,
                                 size_t width, size_t prec)
{
  return _to_str<colvarvalue>(x, width, prec);
}

std::string colvarmodule::to_str(cvm::vector1d<cvm::real> const &x,
                                 size_t width, size_t prec)
{
  return _to_str< cvm::vector1d<cvm::real> >(x, width, prec);
}

std::string colvarmodule::to_str(cvm::matrix2d<cvm::real> const &x,
                                 size_t width, size_t prec)
{
  return _to_str< cvm::matrix2d<cvm::real> >(x, width, prec);
}


std::string colvarmodule::to_str(std::vector<int> const &x,
                                 size_t width, size_t prec)
{
  return _to_str_vector<int>(x, width, prec);
}

std::string colvarmodule::to_str(std::vector<size_t> const &x,
                                 size_t width, size_t prec)
{
  return _to_str_vector<size_t>(x, width, prec);
}

std::string colvarmodule::to_str(std::vector<long int> const &x,
                                 size_t width, size_t prec)
{
  return _to_str_vector<long int>(x, width, prec);
}

std::string colvarmodule::to_str(std::vector<cvm::real> const &x,
                                 size_t width, size_t prec)
{
  return _to_str_vector<cvm::real>(x, width, prec);
}

std::string colvarmodule::to_str(std::vector<cvm::rvector> const &x,
                                 size_t width, size_t prec)
{
  return _to_str_vector<cvm::rvector>(x, width, prec);
}

std::string colvarmodule::to_str(std::vector<cvm::quaternion> const &x,
                                 size_t width, size_t prec)
{
  return _to_str_vector<cvm::quaternion>(x, width, prec);
}

std::string colvarmodule::to_str(std::vector<colvarvalue> const &x,
                                 size_t width, size_t prec)
{
  return _to_str_vector<colvarvalue>(x, width, prec);
}

std::string colvarmodule::to_str(std::vector<std::string> const &x,
                                 size_t width, size_t prec)
{
  return _to_str_vector<std::string>(x, width, prec);
}


std::string cvm::wrap_string(std::string const &s, size_t nchars)
{
  if (!s.size()) {
    return std::string(nchars, ' ');
  } else {
    return ( (s.size() <= nchars) ?
             (s+std::string(nchars-s.size(), ' ')) :
             (std::string(s, 0, nchars)) );
  }
}


// shared pointer to the proxy object
colvarproxy *colvarmodule::proxy = NULL;

// static runtime data
cvm::real colvarmodule::debug_gradients_step_size = 1.0e-07;
int       colvarmodule::errorCode = 0;
int       colvarmodule::log_level_ = 10;
cvm::step_number colvarmodule::it = 0;
cvm::step_number colvarmodule::it_restart = 0;
size_t    colvarmodule::restart_out_freq = 0;
size_t    colvarmodule::cv_traj_freq = 0;
bool      colvarmodule::use_scripted_forces = false;
bool      colvarmodule::scripting_after_biases = true;

// i/o constants
size_t const colvarmodule::it_width = 12;
size_t const colvarmodule::cv_prec  = 14;
size_t const colvarmodule::cv_width = 21;
size_t const colvarmodule::en_prec  = 14;
size_t const colvarmodule::en_width = 21;
const char * const colvarmodule::line_marker = (const char *)
  "----------------------------------------------------------------------\n";
