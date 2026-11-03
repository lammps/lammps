// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>
#include <cstdio>

#include "colvarmodule.h"
#include "colvar_gpu_support.h"
#include "colvarparse.h"
#include "colvarproxy.h"
#include "colvar.h"
#include "colvarbias.h"
#include "colvarbias_abf.h"
#include "colvarbias_abmd.h"
#include "colvarbias_alb.h"
#include "colvarbias_histogram.h"
#include "colvarbias_histogram_reweight_amd.h"
#include "colvarbias_meta.h"
#include "colvarbias_restraint.h"
#include "colvarbias_opes.h"
#include "colvarscript.h"
#include "colvaratoms.h"
#include "colvarcomp.h"
#include "colvars_memstream.h"
#include "colvars_version.h"



std::string colvarmodule::version() const
{
  return std::string(COLVARS_VERSION);
}


int colvarmodule::version_number() const
{
  return version_int;
}


int colvarmodule::patch_version_number() const
{
  return COLVARS_PATCH_VERSION;
}


/// Track usage of Colvars features
class colvarmodule::usage {

public:

  /// Constructor
  usage(colvarmodule *cvmodule_in);

  /// Increment usage count for the given feature; return error if not found
  int cite_feature(std::string const &feature);

  /// Increment usage count for the given paper; return error if not found
  int cite_paper(std::string const &paper);

  /// Generate a report for used features (0 = URL, 1 = BibTeX)
  std::string report(int flag);

  colvarmodule *cvmodule;

protected:

  /// Usage count for each feature
  std::map<std::string, int> feature_count_;

  /// Usage count for each cited paper
  std::map<std::string, int> paper_count_;

  /// URL for each paper
  std::map<std::string, std::string> paper_url_;

  /// BibTeX entry for each paper
  std::map<std::string, std::string> paper_bibtex_;

  /// Map code features to the relevant papers
  std::map<std::string, std::string> feature_paper_map_;

};


namespace {
  constexpr uint32_t colvars_magic_number = 2013813594;
}


colvarmodule::colvarmodule(colvarproxy *proxy_in)
{
  depth_s = 0;
  log_level_ = 10;

  xyz_reader_use_count = 0;

  num_biases_types_used_ =
    reinterpret_cast<void *>(new std::map<std::string, int>());

  restart_version_str.clear();
  restart_version_int = 0;

  usage_ = new usage(this);
  usage_->cite_feature("Colvars module");

  if (proxy != NULL) {
    // TODO relax this error to handle multiple molecules in VMD
    // once the module is not static anymore
    this->error("Error: trying to allocate the collective "
               "variable module twice.\n", COLVARS_BUG_ERROR);
    return;
  }

  proxy = proxy_in; // Pointer to the proxy object
  proxy_static = proxy_in; // Temporary - assume single proxy & module objects

  parse = new colvarparse(this); // Parsing object for global options
  version_int = proxy->get_version_from_string(COLVARS_VERSION);

  this->log(line_marker);
  this->log(
      "Initializing the collective variables module, version " + version() +
      (patch_version_number() ? (" (patch " + this->to_str(patch_version_number()) + ")") : "") +
      ".\n");
  this->log("Please cite Fiorin et al, Mol Phys 2013:\n"
           "  https://doi.org/10.1080/00268976.2013.813594\n"
           "as well as all other papers listed below for individual features used.\n");

  this->log("Summary of compile-time features available in this build:\n");

  std::string cxx_lang_msg("  - C++ language version: " + this->to_str(__cplusplus));
#if defined(_WIN32) && !defined(__CYGWIN__)
  cxx_lang_msg += std::string(" (warning: may not be accurate for this build)");
#endif
  cxx_lang_msg += std::string("\n");
  this->log(cxx_lang_msg);

  if (proxy->check_replicas_enabled() == COLVARS_NOT_IMPLEMENTED) {
    this->log("  - Multiple replicas: not available\n");
  } else {
    if (proxy->check_replicas_enabled() == COLVARS_OK) {
      this->log("  - Multiple replicas: enabled (replica number " +
               to_str(proxy->replica_index() + 1) + " of " + to_str(proxy->num_replicas()) + ")\n");
    } else {
      this->log("  - Multiple replicas: available, but not (yet) enabled\n");
    }
  }

#if defined(LEPTON)
  this->log("  - Lepton custom functions: available\n");
#else
  this->log("  - Lepton custom functions: not available\n");
#endif

#if defined(COLVARS_TCL)
  this->log("  - Tcl interpreter: available\n");
#else
  this->log("  - Tcl interpreter: not available\n");
#endif

  // set initial default values

  binary_restart = false;
  char const *env_var = getenv("COLVARS_BINARY_RESTART");
  if (env_var && atoi(env_var)) {
    binary_restart = true;
  }

  // "it_restart" will be set by the input state file, if any;
  // "it" should be updated by the proxy
  it = it_restart = 0;

  use_scripted_forces = false;
  scripting_after_biases = false;

  colvarmodule::debug_gradients_step_size = 1.0e-05;

  colvarmodule::rotation::monitor_crossings = false;
  colvarmodule::rotation::crossing_threshold = 1.0e-02;

  cv_traj_freq = 100;
  restart_out_freq = proxy->default_restart_frequency();

  cv_traj_write_labels = true;

  // Removes the need for proxy specializations to create this
  proxy->script = new colvarscript(proxy, this);

#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  gpu_calc = nullptr;
#endif
}


colvarmodule * colvarmodule::main()
{
  return proxy_static ? proxy_static->cvmodule : nullptr;
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


int colvarmodule::calc_component_smp(int i)
{
  colvar *x = (*(variables_active_smp()))[i];
  int x_item = (*(variables_active_smp_items()))[i];
  if (this->debug()) {
    this->log("Thread "+this->to_str(proxy->smp_thread_id())+"/"+
             this->to_str(proxy->smp_num_threads())+
             ": calc_component_smp(), i = "+this->to_str(i)+", cv = "+
             x->name+", cvc = "+this->to_str(x_item)+"\n");
  }
  return x->calc_cvcs(x_item, 1);
}


std::vector<colvarbias *> *colvarmodule::biases_active()
{
  return &(biases_active_);
}


size_t colvarmodule::size() const
{
  return colvars.size() + biases.size();
}


void colvarmodule::set_initial_step(step_number it_in)
{
  this->log("Setting initial step number from MD engine: " + this->to_str(it_in) + "\n");
  it = it_restart = it_in;
}


int colvarmodule::read_config_file(char const  *config_filename)
{
  this->log(line_marker);
  this->log("Reading new configuration from file \""+
           std::string(config_filename)+"\":\n");

  // open the configfile
  std::istream &config_s = proxy->input_stream(config_filename,
                                               "configuration file/string");
  if (!config_s) {
    return this->error("Error: in opening configuration file \""+
                      std::string(config_filename)+"\".\n",
                      COLVARS_FILE_ERROR);
  }

  // read the config file into a string
  std::string conf = "";
  std::string line;
  while (parse->read_config_line(config_s, line)) {
    // Delete lines that contain only white space after removing comments
    if (line.find_first_not_of(colvarparse::white_space) != std::string::npos)
      conf.append(line+"\n");
  }
  proxy->close_input_stream(config_filename);

  return parse_config(conf);
}


int colvarmodule::read_config_string(std::string const &config_str)
{
  this->log(line_marker);
  this->log("Reading new configuration:\n");
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
        // Replace Windows newlines with Unix newlines
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
  // Auto-generated additional configuration
  extra_conf.clear();

  // Check that the input has matching braces
  if (colvarparse::check_braces(conf, 0) != COLVARS_OK) {
    return this->error("Error: unmatched curly braces in configuration.\n",
                      COLVARS_INPUT_ERROR);
  }

  // Check that the input has only ASCII characters, and warn otherwise
  parse->check_ascii(conf);

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

  this->log(line_marker);
  this->log("Collective variables module (re)initialized.\n");
  this->log(line_marker);

  if (source_Tcl_script.size() > 0) {
    run_tcl_script(source_Tcl_script);
  }

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


void colvarmodule::config_changed()
{
  cv_traj_write_labels = true;
}


int colvarmodule::parse_global_params(std::string const &conf)
{
  int error_code = COLVARS_OK;
  // TODO document and then echo this keyword
  parse->get_keyval(conf, "logLevel", log_level_, log_level_,
                    colvarparse::parse_silent);
  {
    std::string units;
    if (parse->get_keyval(conf, "units", units)) {
      units = colvarparse::to_lower_cppstr(units);
      error_code |= proxy->set_unit_system(units, (colvars.size() != 0));
    }
  }

  {
    std::string index_file_name;
    size_t pos = 0;
    while (parse->key_lookup(conf, "indexFile", &index_file_name, &pos)) {
      this->log("# indexFile = \""+index_file_name+"\"\n");
      error_code |= read_index_file(index_file_name.c_str());
      index_file_name.clear();
    }
  }

  std::string smp;
  if (parse->get_keyval(conf, "smp", smp)) {
    if (smp == "cvcs" || smp == "on" || smp == "yes") {
      if (proxy->set_smp_mode(colvarproxy_smp::smp_mode_t::cvcs) != COLVARS_OK) {
        this->error("Colvars component-based parallelism is not implemented.\n");
        return COLVARS_INPUT_ERROR;
      }
    } else if (smp == "inner_loop") {
      if (proxy->set_smp_mode(colvarproxy_smp::smp_mode_t::inner_loop) != COLVARS_OK) {
        this->error("SMP parallelism inside the calculation of Colvars components is not implemented.\n");
        return COLVARS_INPUT_ERROR;
      }
    } else if (smp == "gpu") {
      if (proxy->set_smp_mode(colvarproxy_smp::smp_mode_t::gpu) != COLVARS_OK) {
        this->error("GPU parallelism is not implemented.\n");
        return COLVARS_INPUT_ERROR;
      } else {
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
        gpu_calc = std::unique_ptr<colvars_gpu::colvarmodule_gpu_calc>(
          new colvars_gpu::colvarmodule_gpu_calc);
        gpu_calc->init();
        this->log("EXPERIMENTAL GPU parallelism will be applied inside:\n");
        this->log("   - atom groups\n");
#endif
      }
    } else {
      proxy->set_smp_mode(colvarproxy_smp::smp_mode_t::none);
      this->log("SMP parallelism has been disabled.\n");
    }
  } else {
    this->log("SMP parallelism is not set.\n");
    this->log("Available SMP parallelism modes of this proxy are:\n");
    const auto available_smp_modes = proxy->get_available_smp_modes();
    for (size_t i = 0; i < available_smp_modes.size(); ++i) {
      switch (available_smp_modes[i]) {
        case colvarproxy_smp::smp_mode_t::cvcs: {
          this->log("   - cvcs\n");
          break;
        }
        case colvarproxy_smp::smp_mode_t::inner_loop: {
          this->log("   - inner_loop\n");
          break;
        }
        case colvarproxy_smp::smp_mode_t::gpu: {
          this->log("   - gpu\n");
          break;
        }
        case colvarproxy_smp::smp_mode_t::none: {
          this->log("   - none\n");
          break;
        }
      }
    }
    this->log("Set SMP parallelism to the preferred (default) mode to the proxy.\n");
    // Find the proxy's preferred SMP mode if SMP is not defined
    colvarproxy_smp::smp_mode_t preferred_smp_mode = proxy->get_preferred_smp_mode();
    proxy->set_smp_mode(preferred_smp_mode);
    switch (preferred_smp_mode) {
      case colvarproxy_smp::smp_mode_t::cvcs: {
        this->log("SMP parallelism will be applied to Colvars components.\n");
        this->log("  - SMP parallelism: enabled (num. threads = " + to_str(proxy->smp_num_threads()) + ")\n");
        break;
      }
      case colvarproxy_smp::smp_mode_t::inner_loop: {
        this->log("SMP parallelism will be applied inside the Colvars components.\n");
        this->log("  - SMP parallelism: enabled (num. threads = " + to_str(proxy->smp_num_threads()) + ")\n");
        break;
      }
      case colvarproxy_smp::smp_mode_t::gpu: {
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
        gpu_calc = std::unique_ptr<colvars_gpu::colvarmodule_gpu_calc>(
          new colvars_gpu::colvarmodule_gpu_calc);
        gpu_calc->init();
#endif
        this->log("EXPERIMENTAL GPU parallelism will be applied inside:\n");
        this->log("   - atom groups\n");
        break;
      }
      case colvarproxy_smp::smp_mode_t::none: {
        this->log("SMP parallelism is disabled by default.\n");
        break;
      }
    }
  }

  bool b_analysis = true;
  if (parse->get_keyval(conf, "analysis", b_analysis, true, colvarparse::parse_silent)) {
    this->log("Warning: keyword \"analysis\" is deprecated: it is now always set "
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
  if (cv_traj_freq % this->proxy->time_step_factor() != 0) {
    this->error("colvarsTrajFrequency (currently " + this->to_str(cv_traj_freq)
      + ") must be a multiple of the global Colvars timestep multiplier ("
      +  this->to_str(this->proxy->time_step_factor()) + ").\n", COLVARS_INPUT_ERROR);
  }

  parse->get_keyval(conf, "colvarsRestartFrequency",
                    restart_out_freq, restart_out_freq);
  if (restart_out_freq % this->proxy->time_step_factor() != 0) {
    this->error("colvarsRestartFrequency (currently " + this->to_str(restart_out_freq)
      + ") must be a multiple of the global Colvars timestep multiplier ("
      +  this->to_str(this->proxy->time_step_factor()) + ").\n", COLVARS_INPUT_ERROR);
  }

  parse->get_keyval(conf, "scriptedColvarForces",
                    use_scripted_forces, use_scripted_forces);

  parse->get_keyval(conf, "scriptingAfterBiases",
                    scripting_after_biases, scripting_after_biases);

#if defined(COLVARS_TCL)
  parse->get_keyval(conf, "sourceTclFile", source_Tcl_script);
#endif

  if (proxy->engine_name() == "GROMACS" && proxy->version_number() >= 20231003) {
    parse->get_keyval(conf, "defaultInputStateFile", default_input_state_file_,
                      default_input_state_file_);
  }

  return error_code;
}


int colvarmodule::run_tcl_script(std::string const &filename) {

  int result = COLVARS_OK;

#if defined(COLVARS_TCL)
  result = proxy->tcl_run_file(filename);
#endif

  return result;
}


int colvarmodule::parse_colvars(std::string const &conf)
{
  if (this->debug())
    this->log("Initializing the collective variables.\n");

  std::string colvar_conf = "";
  size_t pos = 0;
  while (parse->key_lookup(conf, "colvar", &colvar_conf, &pos)) {

    if (colvar_conf.size()) {
      this->log(line_marker);
      this->increase_depth();
      colvars.push_back(new colvar(this));
      if (((colvars.back())->init(colvar_conf) != COLVARS_OK) ||
          ((colvars.back())->check_keywords(colvar_conf, "colvar") != COLVARS_OK)) {
        this->log("Error while constructing colvar number " +
                 this->to_str(colvars.size()) + " : deleting.");
        delete colvars.back();  // the colvar destructor updates the colvars array
        this->decrease_depth();
        return COLVARS_ERROR;
      }
      this->decrease_depth();
    } else {
      this->error("Error: \"colvar\" keyword found without any configuration.\n", COLVARS_INPUT_ERROR);
      return COLVARS_ERROR;
    }
    this->decrease_depth();
    colvar_conf = "";
  }

  if (pos > 0) {
    // One or more new variables were added
    config_changed();
  }

  if (!colvars.size()) {
    this->log("Warning: no collective variables defined.\n");
  }

  if (colvars.size())
    this->log(line_marker);
  this->log("Collective variables initialized, "+
           this->to_str(colvars.size())+
           " in total.\n");

  return (this->get_error() ? COLVARS_ERROR : COLVARS_OK);
}


bool colvarmodule::check_new_bias(std::string &conf, char const *key)
{
  if (this->get_error() ||
      (biases.back()->check_keywords(conf, key) != COLVARS_OK)) {
    this->log("Error while constructing bias number " +
             this->to_str(biases.size()) + " : deleting.\n");
    delete biases.back(); // the bias destructor updates the biases array
    return true;
  }
  return false;
}


template <class bias_type>
int colvarmodule::parse_biases_type(std::string const &conf,
                                    char const *keyword)
{
  // Allow camel case when calling, but use only lower case for parsing
  std::string const &type_keyword = colvarparse::to_lower_cppstr(keyword);

  // Check how many times this bias keyword was used, set default name
  // accordingly
  std::map<std::string, int> *num_biases_types_used =
    reinterpret_cast<std::map<std::string, int> *>(num_biases_types_used_);
  if (num_biases_types_used->count(type_keyword) == 0) {
    (*num_biases_types_used)[type_keyword] = 0;
  }

  std::string bias_conf = "";
  size_t conf_saved_pos = 0;
  while (parse->key_lookup(conf, keyword, &bias_conf, &conf_saved_pos)) {
    if (bias_conf.size()) {
      this->log(line_marker);
      this->increase_depth();
      int &bias_count = (*num_biases_types_used)[type_keyword];
      biases.push_back(new bias_type(this, type_keyword.c_str()));
      bias_count += 1;
      biases.back()->rank = bias_count;
      biases.back()->init(bias_conf);
      if (this->check_new_bias(bias_conf, keyword) != COLVARS_OK) {
        return COLVARS_ERROR;
      }
      this->decrease_depth();
    } else {
      this->error("Error: keyword \""+std::string(keyword)+"\" found without configuration.\n",
                 COLVARS_INPUT_ERROR);
      return COLVARS_ERROR;
    }
    bias_conf = "";
  }
  if (conf_saved_pos > 0) {
    // One or more new biases were added
    config_changed();
  }
  return COLVARS_OK;
}


int colvarmodule::parse_biases(std::string const &conf)
{
  if (this->debug())
    this->log("Initializing the collective variables biases.\n");

  /// initialize ABF instances
  parse_biases_type<colvarbias_abf>(conf, "abf");

  /// initialize ABMD instances
  parse_biases_type<colvarbias_abmd>(conf, "abmd");

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

  /// initialize reweightaMD instances
  parse_biases_type<colvarbias_reweightaMD>(conf, "reweightaMD");

  /// initialize OPES instances
  parse_biases_type<colvarbias_opes>(conf, "opes_metad");

  if (use_scripted_forces) {
    this->log(line_marker);
    this->increase_depth();
    this->log("User forces script will be run at each bias update.\n");
    this->decrease_depth();
  }

  std::vector<std::string> const time_biases = time_dependent_biases();
  if (time_biases.size() > 1) {
    this->log("WARNING: there are "+this->to_str(time_biases.size())+
             " time-dependent biases with non-zero force parameters:\n"+
             this->to_str(time_biases)+"\n"+
             "Please ensure that their forces do not counteract each other.\n");
  }

  if (num_biases() || use_scripted_forces) {
    this->log(line_marker);
    this->log("Collective variables biases initialized, "+
             this->to_str(num_biases())+" in total.\n");
  } else {
    if (!use_scripted_forces) {
      this->log("No collective variables biases were defined.\n");
    }
  }

  return (this->get_error() ? COLVARS_ERROR : COLVARS_OK);
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
    set_error_bits(COLVARS_INPUT_ERROR);
    parse->clear();
    return get_error();
  }
  return COLVARS_OK;
}


colvarbias * colvarmodule::bias_by_name(std::string const &name)
{
  for (std::vector<colvarbias *>::iterator bi = this->biases.begin();
       bi != this->biases.end();
       bi++) {
    if ((*bi)->name == name) {
      return (*bi);
    }
  }
  return NULL;
}


colvar *colvarmodule::colvar_by_name(std::string const &name)
{
  for (std::vector<colvar *>::iterator cvi = this->colvars.begin();
       cvi != this->colvars.end();
       cvi++) {
    if ((*cvi)->name == name) {
      return (*cvi);
    }
  }
  return NULL;
}

cvm::atom_group *colvarmodule::atom_group_soa_by_name(std::string const& name) {
  for (std::vector<cvm::atom_group *>::iterator agi = this->named_atom_groups_soa.begin();
       agi != this->named_atom_groups_soa.end();
       agi++) {
    if ((*agi)->name == name) {
      return (*agi);
    }
  }
  return nullptr;
}

void colvarmodule::register_named_atom_group_soa(atom_group *ag) {
  named_atom_groups_soa.push_back(ag);
}

void colvarmodule::unregister_named_atom_group_soa(atom_group *ag) {
  for (std::vector<cvm::atom_group *>::iterator agi = named_atom_groups_soa.begin();
       agi != named_atom_groups_soa.end();
       agi++) {
    if (*agi == ag) {
      named_atom_groups_soa.erase(agi);
      break;
    }
  }
}

int colvarmodule::change_configuration(std::string const &bias_name,
                                       std::string const &conf)
{
  // This is deprecated; supported strategy is to delete the bias
  // and parse the new config
  this->increase_depth();
  colvarbias *b;
  b = bias_by_name(bias_name);
  if (b == NULL) {
    this->error("Error: bias not found: " + bias_name);
    return COLVARS_ERROR;
  }
  b->change_configuration(conf);
  this->decrease_depth();
  return (this->get_error() ? COLVARS_ERROR : COLVARS_OK);
}


std::string colvarmodule::read_colvar(std::string const &name)
{
  this->increase_depth();
  colvar *c;
  std::stringstream ss;
  c = colvar_by_name(name);
  if (c == NULL) {
    this->error("Error: colvar not found: " + name);
    return std::string();
  }
  ss << c->value();
  this->decrease_depth();
  return ss.str();
}


cvm::real colvarmodule::energy_difference(std::string const &bias_name,
                                          std::string const &conf)
{
  this->increase_depth();
  colvarbias *b;
  cvm::real energy_diff = 0.;
  b = bias_by_name(bias_name);
  if (b == NULL) {
    this->error("Error: bias not found: " + bias_name);
    return 0.;
  }
  energy_diff = b->energy_difference(conf);
  this->decrease_depth();
  return energy_diff;
}


int colvarmodule::calc()
{
  int error_code = COLVARS_OK;

  if (this->debug()) {
    this->log(line_marker);
    this->log("Collective variables module, step no. "+
             this->to_str(this->step_absolute())+"\n");
  }

  error_code |= calc_colvars();
  error_code |= calc_biases();
  error_code |= update_colvar_forces();

  error_code |= analyze();

  // write trajectory files, if needed
  if (cv_traj_freq && cv_traj_name.size()) {
    error_code |= write_traj_files();
  }

  // write restart files and similar data
  if (restart_out_freq && (this->step_relative() > 0) &&
      ((this->step_absolute() % restart_out_freq) == 0)) {

    if (restart_out_name.size()) {
      // Write restart file, if different from main output
      error_code |= write_restart_file(restart_out_name);
    } else if (output_prefix().size()) {
      error_code |= write_restart_file(output_prefix() + ".colvars.state");
    }

    if (output_prefix().size()) {
      this->increase_depth();
      for (std::vector<colvar *>::iterator cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
        // TODO remove this when corrFunc becomes a bias
        error_code |= (*cvi)->write_output_files();
      }
      for (std::vector<colvarbias *>::iterator bi = biases.begin(); bi != biases.end(); bi++) {
        error_code |= (*bi)->write_state_to_replicas();
      }
      this->decrease_depth();
    }
  }

  // Write output files for biases, at the specified frequency for each
  this->increase_depth();
  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    if ((*bi)->output_freq > 0) {
      if ((this->step_relative() > 0) &&
          ((this->step_absolute() % (*bi)->output_freq) == 0) ) {
        error_code |= (*bi)->write_output_files();
      }
    }
  }
  this->decrease_depth();

  error_code |= end_of_step();

  // TODO move this to a base-class proxy method that calls this function
  error_code |= proxy->end_of_step();

  return error_code;
}


int colvarmodule::calc_colvars()
{
  if (this->debug())
    this->log("Calculating collective variables.\n");
  // calculate collective variables and their gradients

  // First, we need to decide which biases are awake
  // so they can activate colvars as needed
  std::vector<colvarbias *>::iterator bi;
  for (bi = biases.begin(); bi != biases.end(); bi++) {
    int const tsf = (*bi)->get_time_step_factor();
    if (tsf > 1) {
      if (step_absolute() % tsf == 0) {
        (*bi)->enable(colvardeps::f_cvb_awake);
      } else {
        (*bi)->disable(colvardeps::f_cvb_awake);
      }
    }
  }

  int error_code = COLVARS_OK;
  std::vector<colvar *>::iterator cvi;

  // Determine which colvars are active at this iteration
  variables_active()->clear();
  variables_active()->reserve(variables()->size());
  for (cvi = variables()->begin(); cvi != variables()->end(); cvi++) {
    // Wake up or put to sleep variables with MTS
    int tsf = (*cvi)->get_time_step_factor();
    if (tsf > 1) {
      if (step_absolute() % tsf == 0) {
        (*cvi)->enable(colvardeps::f_cv_awake);
      } else {
        (*cvi)->disable(colvardeps::f_cv_awake);
      }
    }

    if ((*cvi)->is_enabled()) {
      variables_active()->push_back(*cvi);
    }
  }

  // if SMP support is available, split up the work
  if (proxy->get_smp_mode() == colvarproxy_smp::smp_mode_t::cvcs) {

    // first, calculate how much work (currently, how many active CVCs) each colvar has

    variables_active_smp()->clear();
    variables_active_smp_items()->clear();

    variables_active_smp()->reserve(variables_active()->size());
    variables_active_smp_items()->reserve(variables_active()->size());

    // set up a vector containing all components
    this->increase_depth();
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
    this->decrease_depth();

    // calculate active colvar components in parallel
    error_code |= proxy->smp_loop(variables_active_smp()->size(), [this](int i) {
        return this->calc_component_smp(i);
      });

    this->increase_depth();
    for (cvi = variables_active()->begin(); cvi != variables_active()->end(); cvi++) {
      error_code |= (*cvi)->collect_cvc_data();
    }
    this->decrease_depth();

  } else if (proxy->get_smp_mode() == colvarproxy_smp::smp_mode_t::gpu) {
    this->increase_depth();
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
    error_code |= gpu_calc->calc_cvs(*variables_active(), this);
#else
    return this->error("GPU calculation is not implemented.\n");
#endif
    this->decrease_depth();
  } else {
    this->increase_depth();
    // calculate colvars one at a time
    for (cvi = variables_active()->begin(); cvi != variables_active()->end(); cvi++) {
      error_code |= (*cvi)->calc();
      if (this->get_error()) {
        return COLVARS_ERROR;
      }
    }
    this->decrease_depth();
  }

  error_code |= this->get_error();
  return error_code;
}


int colvarmodule::calc_biases()
{
  // update the biases and communicate their forces to the collective
  // variables
  if (this->debug() && num_biases())
    this->log("Updating collective variable biases.\n");

  // set biasing forces to zero before biases are calculated and summed over
  for (std::vector<colvar *>::iterator cvi = colvars.begin();
       cvi != colvars.end(); cvi++) {
    (*cvi)->reset_bias_force();
  }

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

  bool biases_need_main_thread = false;
  for (bi = biases_active()->begin(); bi != biases_active()->end(); bi++) {
    if ((*bi)->replica_share_freq() > 0) {
      // Biases that share data with replicas need read/write access to I/O or MPI
      biases_need_main_thread = true;
    }
  }

  // If SMP support is available, split up the work (unless biases need to use main thread's memory)
  if (proxy->get_smp_mode() == colvarproxy::smp_mode_t::cvcs && !biases_need_main_thread) {

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

    // Straight loop over biases on a single thread
    this->increase_depth();
    for (bi = biases_active()->begin(); bi != biases_active()->end(); bi++) {
      error_code |= (*bi)->update();
      if (this->get_error()) {
        this->decrease_depth();
        return error_code;
      }
    }
    this->decrease_depth();
  }

  for (bi = biases_active()->begin(); bi != biases_active()->end(); bi++) {
    total_bias_energy += (*bi)->get_energy();
  }

  return error_code;
}


int colvarmodule::update_colvar_forces()
{
  int error_code = COLVARS_OK;

  std::vector<colvar *>::iterator cvi;
  std::vector<colvarbias *>::iterator bi;

  // sum the forces from all biases for each collective variable
  if (this->debug() && num_biases())
    this->log("Collecting forces from all biases.\n");
  this->increase_depth();
  for (bi = biases_active()->begin(); bi != biases_active()->end(); bi++) {
    error_code |= (*bi)->communicate_forces();
  }
  this->decrease_depth();

  if (use_scripted_forces && scripting_after_biases) {
    error_code |= calc_scripted_forces();
  }

  // Now we have collected energies from both built-in and scripted biases
  if (this->debug())
    this->log("Adding total bias energy: " + this->to_str(total_bias_energy) + "\n");
  proxy->add_energy(total_bias_energy);

  cvm::real total_colvar_energy = 0.0;
  // sum up the forces for each colvar, including wall forces
  // and integrate any internal
  // equation of motion (extended system)
  if (this->debug())
    this->log("Updating the internal degrees of freedom "
             "of colvars (if they have any).\n");
  this->increase_depth();
  for (cvi = variables()->begin(); cvi != variables()->end(); cvi++) {
    // Inactive colvars will only reset their forces and return 0 energy
    total_colvar_energy += (*cvi)->update_forces_energy();
  }
  this->decrease_depth();
  if (this->debug())
    this->log("Adding total colvar energy: " + this->to_str(total_colvar_energy) + "\n");
  proxy->add_energy(total_colvar_energy);

  // make collective variables communicate their forces to their
  // coupled degrees of freedom (i.e. atoms)
  if (this->debug())
    this->log("Communicating forces from the colvars to the atoms.\n");
  this->increase_depth();
  if (proxy->get_smp_mode() == colvarproxy_smp::smp_mode_t::gpu) {
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
    error_code |= gpu_calc->apply_forces(*variables_active(), this);
#endif
  } else {
    for (cvi = variables_active()->begin(); cvi != variables_active()->end(); cvi++) {
      if ((*cvi)->is_enabled(colvardeps::f_cv_apply_force)) {
        (*cvi)->communicate_forces();
        if (this->get_error()) {
          return COLVARS_ERROR;
        }
      }
    }
  }
  this->decrease_depth();

  return error_code;
}


int colvarmodule::calc_scripted_forces()
{
  // Run user force script, if provided,
  // potentially adding scripted forces to the colvars
  int res;
  res = proxy->run_force_callback();
  if (res == COLVARS_NOT_IMPLEMENTED) {
    this->error("Colvar forces scripts are not implemented.");
    return COLVARS_NOT_IMPLEMENTED;
  }
  if (res != COLVARS_OK) {
    this->error("Error running user colvar forces script");
    return COLVARS_ERROR;
  }
  return COLVARS_OK;
}


int colvarmodule::write_restart_file(std::string const &out_name)
{
  this->log("Saving collective variables state to \""+out_name+"\".\n");
  std::ostream &restart_out_os = proxy->output_stream(out_name, "state file");
  if (!restart_out_os) return COLVARS_FILE_ERROR;

  if (binary_restart) {
    cvm::memory_stream mem_os;
    if (!write_state(mem_os)) {
      return this->error("Error: in writing binary state information to file.\n", COLVARS_ERROR);
    }
    if (!restart_out_os.write(reinterpret_cast<char *>(mem_os.output_buffer()),
                              mem_os.length())) {
      return this->error("Error: in writing restart file.\n", COLVARS_FILE_ERROR);
    }
  } else {
    if (!write_state(restart_out_os)) {
      return this->error("Error: in writing restart file.\n", COLVARS_FILE_ERROR);
    }
  }

  proxy->close_output_stream(out_name);

  // Take the opportunity to flush colvars.traj

  return (this->get_error() ? COLVARS_ERROR : COLVARS_OK);
}


int colvarmodule::write_restart_string(std::string &output)
{
  this->log("Saving state to output buffer.\n");
  std::ostringstream os;
  if (!write_state(os)) {
    return this->error("Error: in writing restart to buffer.\n", COLVARS_FILE_ERROR);
  }
  output = os.str();
  return COLVARS_OK;
}


int colvarmodule::write_traj_files()
{
  int error_code = COLVARS_OK;

  if (this->debug()) {
    this->log("colvarmodule::write_traj_files()\n");
  }

  std::ostream &cv_traj_os = proxy->output_stream(cv_traj_name,
                                                  "colvars trajectory");

  if (!cv_traj_os) {
    return COLVARS_FILE_ERROR;
  }

  // Write labels in the traj file at beginning and then every 1000 lines
  if ( (this->step_relative() == 0) || cv_traj_write_labels ||
       ((this->step_absolute() % (cv_traj_freq * 1000)) == 0) ) {
    error_code |=
      write_traj_label(cv_traj_os) ? COLVARS_OK : COLVARS_FILE_ERROR;
    cv_traj_write_labels = false;
  }

  if (this->debug()) {
    proxy->flush_output_stream(cv_traj_name);
  }

  if ((this->step_absolute() % cv_traj_freq) == 0) {
    error_code |= write_traj(cv_traj_os) ? COLVARS_OK : COLVARS_FILE_ERROR;
  }

  if (this->debug()) {
    proxy->flush_output_stream(cv_traj_name);
  }

  if (restart_out_freq && ((this->step_absolute() % restart_out_freq) == 0)) {
    this->log("Synchronizing (emptying the buffer of) trajectory file \""+
             cv_traj_name+"\".\n");
    error_code |= proxy->flush_output_stream(cv_traj_name);
  }

  return error_code;
}


int colvarmodule::analyze()
{
  if (this->debug()) {
    this->log("colvarmodule::analyze(), step = "+this->to_str(it)+".\n");
  }

  // perform colvar-specific analysis
  for (std::vector<colvar *>::iterator cvi = variables_active()->begin();
       cvi != variables_active()->end();
       cvi++) {
    this->increase_depth();
    (*cvi)->analyze();
    this->decrease_depth();
  }

  // perform bias-specific analysis
  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    this->increase_depth();
    (*bi)->analyze();
    this->decrease_depth();
  }

  return (this->get_error() ? COLVARS_ERROR : COLVARS_OK);
}


int colvarmodule::end_of_step()
{
  if (this->debug()) {
    this->log("colvarmodule::end_of_step(), step = "+this->to_str(it)+".\n");
  }

  for (std::vector<colvar *>::iterator cvi = variables_active()->begin();
       cvi != variables_active()->end();
       cvi++) {
    this->increase_depth();
    (*cvi)->end_of_step();
    this->decrease_depth();
  }

  // perform bias-specific analysis
  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    this->increase_depth();
    (*bi)->end_of_step();
    this->decrease_depth();
  }

  return (this->get_error() ? COLVARS_ERROR : COLVARS_OK);
}


int colvarmodule::update_engine_parameters()
{
  if (size() == 0) {
    // No-op if no variables or biases are defined
    return this->get_error();
  }
  if (proxy->simulation_running()) {
    this->log("Current simulation parameters: initial step = " + this->to_str(it) +
             ", integration timestep = " + this->to_str(dt()) + "\n");
  }
  this->log("Updating atomic parameters (masses, charges, etc).\n");
  for (std::vector<colvar *>::iterator cvi = variables()->begin(); cvi != variables()->end();
       cvi++) {
    (*cvi)->setup();
  }
  return (this->get_error() ? COLVARS_ERROR : COLVARS_OK);
}


colvarmodule::~colvarmodule()
{
  if ((proxy->smp_thread_id() < 0) ||  // not using threads
      (proxy->smp_thread_id() == 0)) { // or this is thread 0

    reset();

    // Delete contents of static arrays
    colvarbias::delete_features();
    colvar::delete_features();
    colvar::cvc::delete_features();
    atom_group::delete_features();

    delete
      reinterpret_cast<std::map<std::string, int> *>(num_biases_types_used_);
    num_biases_types_used_ = NULL;

    delete parse;
    parse = NULL;

    delete usage_;
    usage_ = NULL;

    // The proxy object will be deallocated last (if at all)
    proxy = NULL;
  }
}


int colvarmodule::reset()
{
  parse->clear();

  // Iterate backwards because we are deleting the elements as we go
  while (!biases.empty()) {
    colvarbias* tail = biases.back();
    delete tail; // the bias destructor updates the biases array
  }
  biases.clear();
  biases_active_.clear();

  // Reset counters tracking usage of each bias type
  reinterpret_cast<std::map<std::string, int> *>(num_biases_types_used_)->clear();

  // Iterate backwards because we are deleting the elements as we go
  while (!colvars.empty()) {
    colvar* cvi = colvars.back();
    delete cvi; // the colvar destructor updates the colvars array
  };
  colvars.clear();

  reset_index_groups();

  proxy->flush_output_streams();
  proxy->reset();

  clear_error();

  return COLVARS_OK;
}


int colvarmodule::setup_input()
{
  if (proxy->input_prefix().empty() && (!proxy->input_stream_exists("input state string")) &&
      input_state_buffer_.empty()) {
    // If no input sources have been defined up to this point, use defaultInputStateFile
    proxy->set_input_prefix(default_input_state_file_);
  }

  if (!proxy->input_prefix().empty()) {

    // Read state from a file

    std::string restart_in_name(proxy->input_prefix() + std::string(".colvars.state"));
    std::istream *input_is = &(proxy->input_stream(restart_in_name, "restart file/channel", false));
    if (!*input_is) {
      // Try without the suffix ".colvars.state"
      restart_in_name = proxy->input_prefix();
      input_is = &(proxy->input_stream(restart_in_name, "restart file/channel"));
      if (!*input_is) {
        // Error message has already been printed, return now
        return COLVARS_FILE_ERROR;
      }
    }

    // Now that the file has been opened, clear this field so that this block
    // will not be executed twice
    proxy->set_input_prefix("");

    this->log(line_marker);

    input_is->seekg(0, std::ios::end);
    size_t const file_size = input_is->tellg();
    input_is->seekg(0, std::ios::beg);

    bool binary_state_file = false;

    uint32_t file_magic_number = 0;
    if (file_size > sizeof(uint32_t)) {
      if (input_is->read(reinterpret_cast<char *>(&file_magic_number), sizeof(uint32_t))) {
        if (file_magic_number == colvars_magic_number) {
          binary_state_file = true;
        }
        input_is->seekg(0, std::ios::beg);
      }
    }

    if (binary_state_file) {
      this->log("Loading state from binary file \"" + restart_in_name + "\".\n");
      // TODO integrate istream.read() into memory_stream to avoid copying
      auto *buf = new unsigned char[file_size];
      if (input_is->read(reinterpret_cast<char *>(buf), file_size)) {
        cvm::memory_stream mem_is(file_size, buf);
        if (!read_state(mem_is)) {
          input_is->setstate(std::ios::failbit);
          this->error("Error: cannot interpret contents of binary file \"" + restart_in_name +
                         "\".\n",
                     COLVARS_INPUT_ERROR);
        }
      } else {
        this->error("Error: cannot read from binary file \"" + restart_in_name + "\".\n",
                   COLVARS_INPUT_ERROR);
      }
      delete[] buf;
    } else {
      this->log("Loading state from text file \"" + restart_in_name + "\".\n");
      read_state(*input_is);
    }
    this->log(line_marker);

    // Now that an explicit state file was read, we shall ignore any other restart info
    if (proxy->input_stream_exists("input state string")) {
      proxy->delete_input_stream("input state string");
    }
    input_state_buffer_.clear();

    proxy->delete_input_stream(restart_in_name);
  }

  if (proxy->input_stream_exists("input state string")) {

    if (!input_state_buffer_.empty()) {
      return this->error("Error: formatted/text and unformatted/binary input state buffers are "
                        "defined at the same time.\n",
                        COLVARS_BUG_ERROR);
    }

    this->log(line_marker);
    this->log("Loading state from formatted string.\n");
    read_state(proxy->input_stream("input state string"));
    this->log(line_marker);

    proxy->delete_input_stream("input state string");
  }

  if (!input_state_buffer_.empty()) {
    this->log(line_marker);
    this->log("Loading state from unformatted memory.\n");
    cvm::memory_stream ms(input_state_buffer_.size(), input_state_buffer_.data());
    read_state(ms);
    this->log(line_marker);

    input_state_buffer_.clear();
  }

  default_input_state_file_.clear();

  return get_error();
}


int colvarmodule::setup_output()
{
  int error_code = COLVARS_OK;

  // output state file (restart)
  restart_out_name = proxy->restart_output_prefix().size() ?
    std::string(proxy->restart_output_prefix()+".colvars.state") :
    std::string("");

  std::string const state_file_format(binary_restart ? " (binary format)" : "");

  if (restart_out_name.size()) {
    this->log("The restart output state file" + state_file_format + " will be \""+
             restart_out_name+"\".\n");
  }

  if (output_prefix() != proxy->output_prefix()) {
    output_prefix() = proxy->output_prefix();
    if (output_prefix().size()) {
      this->log("The final output state file will be \"" +
               (output_prefix().size() ? std::string(output_prefix() + ".colvars.state")
                                       : std::string("colvars.state")) +
               "\".\n");
    }

    if (proxy->output_stream_exists(cv_traj_name)) {
      // Close old file
      proxy->close_output_stream(cv_traj_name);
      cv_traj_write_labels = true;
    }

    cv_traj_name =
        (output_prefix().size() ? std::string(output_prefix() + ".colvars.traj") : std::string(""));

    for (std::vector<colvarbias *>::iterator bi = biases.begin();
         bi != biases.end();
         bi++) {
      error_code |= (*bi)->setup_output();
    }
  }

  return error_code;
}


std::string colvarmodule::state_file_prefix(char const *filename)
{
  std::string const filename_str(filename);
  std::string const prefix =
    filename_str.substr(0, filename_str.find(".colvars.state"));
  if (prefix.size() == 0) {
    this->error("Error: invalid filename/prefix value \""+filename_str+"\".",
               COLVARS_INPUT_ERROR);
  }
  return prefix;
}


template <typename IST> IST & colvarmodule::read_state_template_(IST &is)
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

      restart_version_str.clear();
      restart_version_int = 0;
      parse->get_keyval(restart_conf, "version",
                        restart_version_str, std::string(""),
                        colvarparse::parse_restart);
      if (restart_version_str.size()) {
        // Initialize integer version number of this restart file
        restart_version_int =
          proxy->get_version_from_string(restart_version_str.c_str());
      }

      if (restart_version() != version()) {
        this->log("This state file was generated with version " + restart_version() + "\n");
        if (std::is_same<IST, cvm::memory_stream>::value) {
          this->log("Warning: compatibility between differetn Colvars versions is not "
                   "guaranteed for unformatted (binary) state files.\n");
        }
      }

      if (restart_version_number() < 20160810) {
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
          this->error("Error: the state file has units \""+units_restart+
                     "\", but the current unit system is \""+proxy->units+
                     "\".\n", COLVARS_INPUT_ERROR);
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


std::istream & colvarmodule::read_state(std::istream &is)
{
  return read_state_template_<std::istream>(is);
}


cvm::memory_stream &colvarmodule::read_state(cvm::memory_stream &is)
{
  uint32_t file_magic_number = 0;
  if (!(is >> file_magic_number)) {
    return is;
  }
  if (file_magic_number == colvars_magic_number) {
    return read_state_template_<cvm::memory_stream>(is);
  } else {
    is.setstate(std::ios::failbit);
    this->error("Error: magic number of binary file (" +
                   this->to_str(static_cast<size_t>(file_magic_number)) +
                   ") does not match the expected magic number for a Colvars state file (" +
                   this->to_str(static_cast<size_t>(colvars_magic_number)) + ").\n",
               COLVARS_INPUT_ERROR);
  }
  return is;
}


int colvarmodule::set_input_state_buffer(size_t n, unsigned char *buf)
{
  input_state_buffer_.clear();
  std::copy(buf, buf + n, std::back_inserter(input_state_buffer_));
  return COLVARS_OK;
}


int colvarmodule::set_input_state_buffer(std::vector<unsigned char> &buf)
{
  input_state_buffer_ = std::move(buf);
  return COLVARS_OK;
}


std::istream & colvarmodule::read_objects_state(std::istream &is)
{
  auto pos = is.tellg();
  std::string word;

  while (is) {
    pos = is.tellg();

    if (is >> word) {

      is.seekg(pos);

      if (word == "colvar") {

        this->increase_depth();
        for (std::vector<colvar *>::iterator cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
          if (!((*cvi)->read_state(is))) {
            // Here an error signals that the variable is a match, but the
            // state is corrupt; otherwise, the variable rewinds is silently
            this->error("Error: in reading state for collective variable \"" +
                           (*cvi)->name + "\" at position " + this->to_str(is.tellg()) +
                           " in stream.\n",
                       COLVARS_INPUT_ERROR);
          }
          if (is.tellg() > pos)
            break; // found it
        }
        this->decrease_depth();

      } else {

        this->increase_depth();
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
            this->error("Error: in reading state for bias \"" + (*bi)->name + "\" at position " +
                           this->to_str(is.tellg()) + " in stream.\n",
                       COLVARS_INPUT_ERROR);
          }
          if (is.tellg() > pos)
            break; // found it
        }
        this->decrease_depth();
      }
    }

    if (is.tellg() == pos) {
      // This block has not been read by any object: discard it, print warning and move on
      std::string conf;
      is >> colvarparse::read_block(word, &conf);
      std::string name;
      colvarparse parse_tmp(this);
      if (parse_tmp.key_lookup(conf, "name", &name)) {
        cvm::log("Warning: input state file contains unknown " + word + " object \"" + name + "\", discarding.\n");
      } else {
        cvm::log("Warning: input state file contains unknown " + word + " object, discarding.\n");
      }
    }

    if (!is) break;
  }

  return is;
}


cvm::memory_stream &colvarmodule::read_objects_state(cvm::memory_stream &is)
{
  // An unformatted stream must match the objects' exact configuration
  this->increase_depth();
  for (std::vector<colvar *>::iterator cvi = colvars.begin(); cvi != colvars.end(); cvi++) {
    if (!(*cvi)->read_state(is)) {
      return is;
    }
  }
  for (std::vector<colvarbias *>::iterator bi = biases.begin(); bi != biases.end(); bi++) {
    if (!(*bi)->read_state(is)) {
      return is;
    }
  }
  this->decrease_depth();
  return is;
}


int colvarmodule::print_total_forces_errning(bool warn_total_forces)
{
  if (warn_total_forces) {
    this->log(line_marker);
    this->log("WARNING: The definition of system forces has changed.  Please see:\n");
    this->log("  https://colvars.github.io/README-totalforce.html\n");
    // update this ahead of time in this special case
    output_prefix() = proxy->input_prefix();
    this->log("All output files will now be saved with the prefix \""+output_prefix()+".tmp.*\".\n");
    this->log("Please review the important warning above. After that, you may rename:\n\
\""+output_prefix()+".tmp.colvars.state\"\n\
to:\n\
\""+proxy->input_prefix()+".colvars.state\"\n\
and load it to continue this simulation.\n");
    output_prefix() = output_prefix()+".tmp";
    write_restart_file(output_prefix()+".colvars.state");
    return this->error("Exiting with error until issue is addressed.\n",
                      COLVARS_INPUT_ERROR);
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
  this->increase_depth();
  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    // Only write output files if they have not already been written this time step
    if ((*bi)->output_freq == 0    ||
        this->step_relative() == 0  ||
        (this->step_absolute() % (*bi)->output_freq) != 0) {
      error_code |= (*bi)->write_output_files();
    }
    error_code |= (*bi)->write_state_to_replicas();
  }
  this->decrease_depth();
  return error_code;
}


int colvarmodule::read_traj(char const *traj_filename,
                            long        traj_read_begin,
                            long        traj_read_end)
{
  this->log("Opening trajectory file \""+
           std::string(traj_filename)+"\".\n");
  // NB: this function is not currently used, but when it will it should
  // retain the ability for direct file-based access (in case traj files
  // exceed memory)
  std::ifstream traj_is(traj_filename);

  while (true) {
    while (true) {

      std::string line("");

      do {
        if (!colvarparse::getline_nocomments(traj_is, line)) {
          this->log("End of file \""+std::string(traj_filename)+
                   "\" reached, or corrupted file.\n");
          traj_is.close();
          return false;
        }
      } while (line.find_first_not_of(colvarparse::white_space) == std::string::npos);

      std::istringstream is(line);

      if (!(is >> it)) return false;

      if ( (it < traj_read_begin) ) {

        if ((it % 1000) == 0)
          std::fprintf(stderr, "Skipping trajectory step %lld                    \r", (long long)it);

        continue;

      } else {

        if ((it % 1000) == 0)
          std::fprintf(stderr, "Reading from trajectory, step = %lld                    \r", (long long)it);

        if ( (traj_read_end > traj_read_begin) &&
             (it > traj_read_end) ) {
          std::fprintf(stderr, "\n");
          this->error("Reached the end of the trajectory, "
                     "read_end = "+this->to_str(traj_read_end)+"\n",
                     COLVARS_FILE_ERROR);
          return COLVARS_ERROR;
        }

        for (std::vector<colvar *>::iterator cvi = colvars.begin();
             cvi != colvars.end();
             cvi++) {
          if (!(*cvi)->read_traj(is)) {
            this->error("Error: in reading colvar \""+(*cvi)->name+
                       "\" from trajectory file \""+
                       std::string(traj_filename)+"\".\n",
                       COLVARS_FILE_ERROR);
            return COLVARS_ERROR;
          }
        }

        break;
      }
    }
  }
  return (this->get_error() ? COLVARS_ERROR : COLVARS_OK);
}


template <typename OST> OST &colvarmodule::write_state_template_(OST &os)
{
  bool const formatted = !std::is_same<OST, cvm::memory_stream>::value;

  std::ostringstream oss;
  oss.setf(std::ios::scientific, std::ios::floatfield);
  oss << "  step " << std::setw(it_width)
      << it << "\n"
      << "  dt " << dt() << "\n"
      << "  version " << std::string(COLVARS_VERSION) << "\n";
  if (proxy->units.size() > 0) {
    oss << "  units " << proxy->units << "\n";
  }

  os << std::string("configuration");
  if (formatted) os << " {\n";
  os << oss.str();
  if (formatted) os << "}\n\n";

  int error_code = COLVARS_OK;

  this->increase_depth();
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
  this->decrease_depth();

  if (error_code != COLVARS_OK) {
    // TODO make this function return an int instead
    os.setstate(std::ios::failbit);
  }

  return os;
}


std::ostream &colvarmodule::write_state(std::ostream &os)
{
  return write_state_template_<std::ostream>(os);
}


cvm::memory_stream &colvarmodule::write_state(cvm::memory_stream &os)
{
  if (os << colvars_magic_number) {
    write_state_template_<cvm::memory_stream>(os);
  }
  return os;
}


int colvarmodule::write_state_buffer(std::vector<unsigned char> &buffer)
{
  cvm::memory_stream os(buffer);
  if (os << colvars_magic_number) {
    write_state_template_<cvm::memory_stream>(os);
  }
  return os ? COLVARS_OK : COLVARS_ERROR;
}


std::ostream &colvarmodule::write_traj_label(std::ostream &os)
{
  os.setf(std::ios::scientific, std::ios::floatfield);

  os << "# " << this->wrap_string("step", this->it_width-2)
     << " ";

  this->increase_depth();
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

  this->decrease_depth();
  return os;
}


std::ostream & colvarmodule::write_traj(std::ostream &os)
{
  os.setf(std::ios::scientific, std::ios::floatfield);

  os << std::setw(this->it_width) << it
     << " ";

  this->increase_depth();
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

  this->decrease_depth();
  return os;
}


void colvarmodule::log(std::string const &message, int min_log_level)
{
  if (this->log_level() < min_log_level) return;

  std::string const trailing_newline = (message.size() > 0) ?
    (message[message.size()-1] == '\n' ? "" : "\n") : "";

  if (proxy) { // This may run before cvm construction
    // allow logging when the module is not fully initialized
    size_t const d = depth();
    if (d > 0) {
      proxy->log((std::string(2*d, ' ')) + message + trailing_newline);
    } else {
      proxy->log(message + trailing_newline);
    }
  } else { // Safe without a proxy pointer
    // std::cout << message + trailing_newline;
    std::printf("%s%s", message.c_str(), trailing_newline.c_str());
  }
}


void colvarmodule::increase_depth()
{
  (depth())++;
}


void colvarmodule::decrease_depth()
{
  if (depth() > 0) {
    (depth())--;
  }
}


size_t & colvarmodule::depth()
{
  // NOTE: do not call log() or error() here, to avoid recursion
  if (proxy->get_smp_mode() == colvarproxy::smp_mode_t::cvcs) {
    int const nt = proxy->smp_num_threads();
    if (int(this->depth_v.size()) != nt) {
      proxy->smp_lock();
      // update array of depths
      if (this->depth_v.size() > 0) { this->depth_s = this->depth_v[0]; }
      this->depth_v.clear();
      this->depth_v.assign(nt, this->depth_s);
      proxy->smp_unlock();
    }
    return this->depth_v[proxy->smp_thread_id()];
  }
  return this->depth_s;
}


void colvarmodule::set_error_bits(int code)
{
  if (code < 0) {
    this->log("Error: set_error_bits() received negative error code.\n");
    return;
  }
  if (proxy) { // This may run before cvm construction
    proxy->smp_lock();
    errorCode |= code | COLVARS_ERROR;
    proxy->smp_unlock();
  }
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
  proxy->clear_error_msgs();
}


int colvarmodule::error(std::string const &message, int code)
{
  set_error_bits(code >= 0 ? code : COLVARS_ERROR);

  std::string const trailing_newline = (message.size() > 0) ?
    (message[message.size()-1] == '\n' ? "" : "\n") : "";

  std::string prefix = "Error: ";
  if (message.size() >= 7) {
    if (message.substr(0, 7) == prefix) {
      prefix.clear();
    }
  }

  if (proxy) { // This may run before cvm construction
    size_t const d = depth();
    if (d > 0) {
      proxy->error((std::string(2*d, ' ')) + prefix + message + trailing_newline);
    } else {
      proxy->error(message + trailing_newline);
    }
    return get_error();

  } else { // Safe without a proxy pointer
    std::fprintf(stderr, "%s%s%s", prefix.c_str(), message.c_str(), trailing_newline.c_str());
    return code;
  }
}


int colvarmodule::read_index_file(char const *filename)
{
  std::istream &is = proxy->input_stream(filename, "index file");

  if (!is) {
    return COLVARS_FILE_ERROR;
  } else {
    index_file_names.push_back(std::string(filename));
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
      return this->error("Error: in parsing index file \""+
                        std::string(filename)+"\".\n",
                        COLVARS_INPUT_ERROR);
    }

    std::vector<int> *old_index_group = index_groups[index_of_group];
    std::vector<int> *new_index_group = new std::vector<int>();

    int atom_number = 1;
    std::streampos pos = is.tellg();
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
        return this->error("Error: the index group \""+group_name+
                          "\" was redefined.\n", COLVARS_INPUT_ERROR);
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

  this->log("The following index groups are currently defined:\n");
  size_t i = 0;
  for ( ; i < index_group_names.size(); i++) {
    this->log("  "+(index_group_names[i])+" ("+
             this->to_str((index_groups[i])->size())+" atoms)\n");
  }

  return proxy->close_input_stream(filename);
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
  index_file_names.clear();
  return COLVARS_OK;
}

int colvarmodule::load_coords(char const *file_name,
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

  std::vector<cvm::atom_pos> sorted_pos(atoms->size(), cvm::rvector(0.0));

  // Differentiate between PDB and XYZ files
  if (colvarparse::to_lower_cppstr(ext) == std::string(".xyz")) {
    if (pdb_field.size() > 0) {
      return this->error("Error: PDB column may not be specified "
                        "for XYZ coordinate files.\n", COLVARS_INPUT_ERROR);
    }
    // For XYZ files, use internal parser
    error_code |= this->load_coords_xyz(file_name, &sorted_pos, atoms);
  } else {
    // Otherwise, call proxy function for PDB
    error_code |= proxy->load_coords_pdb(file_name, sorted_pos, atoms->sorted_ids(), pdb_field,
                                         pdb_field_value);
  }

  if (error_code != COLVARS_OK) return error_code;

  std::vector<int> const &map = atoms->sorted_ids_map();
  for (size_t i = 0; i < atoms->size(); i++) {
    (*pos)[map[i]] = sorted_pos[i];
  }

  return error_code;
}

int colvarmodule::load_coords_xyz(char const *filename,
                         std::vector<rvector> *pos,
                         cvm::atom_group *atoms,
                         bool keep_open)
{
  std::istream &xyz_is = proxy->input_stream(filename, "XYZ file");
  size_t natoms;
  char symbol[256];
  std::string line;
  cvm::real x = 0.0, y = 0.0, z = 0.0;

  std::string const error_msg("Error: cannot parse XYZ file \""+
                              std::string(filename)+"\".\n");

  if ( ! (xyz_is >> natoms) ) {
      // Return silent error when reaching the end of multi-frame files
      return keep_open ? COLVARS_NO_SUCH_FRAME : this->error(error_msg, COLVARS_INPUT_ERROR);
  }

  ++xyz_reader_use_count;
  if (xyz_reader_use_count < 2) {
    this->log("Warning: beginning from 2019-11-26 the XYZ file reader assumes Angstrom units.\n");
  }

  if (xyz_is.good()) {
    // skip comment line
    this->getline(xyz_is, line);
    this->getline(xyz_is, line);
    xyz_is.width(255);
  } else {
    proxy->close_input_stream(filename);
    return this->error(error_msg, COLVARS_INPUT_ERROR);
  }

  if (pos->size() > natoms) {
    proxy->close_input_stream(filename);
    return this->error("File \"" + std::string(filename) + "\" contains fewer atoms (" + this->to_str(natoms)
      + ") than expected (" + this->to_str(pos->size()) + ").", COLVARS_INPUT_ERROR);
  }

  std::vector<atom_pos>::iterator pos_i = pos->begin();
  size_t xyz_natoms = 0;
  if (pos->size() < natoms) { // Use specified indices
    int next = 0; // indices are zero-based
    if (!atoms) {
      // In the other branch of this test, reading all positions from the file,
      // a valid atom group pointer is not necessary
      return this->error("Trying to read partial positions with invalid atom group pointer",
                        COLVARS_BUG_ERROR);
    }

    if (static_cast<unsigned int>(atoms->sorted_ids().back()) > natoms) {
      proxy->close_input_stream(filename);
      return this->error("File \"" + std::string(filename) + "\" contains fewer atoms (" + this->to_str(natoms)
        + ") than expected (" + this->to_str(atoms->sorted_ids().back()) + ").", COLVARS_INPUT_ERROR);
    }

    std::vector<int>::const_iterator index = atoms->sorted_ids().begin();

    for ( ; pos_i != pos->end() ; pos_i++, index++) {
      while ( next < *index ) {
        this->getline(xyz_is, line);
        next++;
      }
      if (xyz_is.good()) {
        xyz_is >> symbol;
        xyz_is >> x >> y >> z;
        // XYZ files are assumed to be in Angstrom (as eg. VMD will)
        (*pos_i)[0] = proxy->angstrom_to_internal(x);
        (*pos_i)[1] = proxy->angstrom_to_internal(y);
        (*pos_i)[2] = proxy->angstrom_to_internal(z);
        xyz_natoms++;
      } else {
        proxy->close_input_stream(filename);
        return this->error(error_msg, COLVARS_INPUT_ERROR);
      }
    }

  } else {          // Use all positions

    for ( ; pos_i != pos->end() ; pos_i++) {
      if (xyz_is.good()) {
        xyz_is >> symbol;
        xyz_is >> x >> y >> z;
        (*pos_i)[0] = proxy->angstrom_to_internal(x);
        (*pos_i)[1] = proxy->angstrom_to_internal(y);
        (*pos_i)[2] = proxy->angstrom_to_internal(z);
        xyz_natoms++;
      } else {
        proxy->close_input_stream(filename);
        return this->error(error_msg, COLVARS_INPUT_ERROR);
      }
    }
  }

  if (xyz_natoms != pos->size()) {
    proxy->close_input_stream(filename);
    return this->error("Error: The number of positions read from file \""+
                      std::string(filename)+"\" does not match the number of "+
                      "positions required: "+this->to_str(xyz_natoms)+" vs. "+
                      this->to_str(pos->size())+".\n", COLVARS_INPUT_ERROR);
  }

  if (keep_open) {
    return COLVARS_OK;
  } else {
    return proxy->close_input_stream(filename);
  }
}

// Wrappers to proxy functions: these may go in the future


cvm::real colvarmodule::dt()
{
  return proxy->dt();
}


void colvarmodule::request_total_force()
{
  proxy->request_total_force(true);
}


cvm::real colvarmodule::rand_gaussian(void)
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


template<typename T> std::string _to_str_vector(T const &x,
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

std::string colvarmodule::to_str(const void* ptr) {
  return _to_str<const void*>(ptr, 0, 0);
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
  return _to_str_vector(x, width, prec);
}

std::string colvarmodule::to_str(std::vector<size_t> const &x,
                                 size_t width, size_t prec)
{
  return _to_str_vector(x, width, prec);
}

std::string colvarmodule::to_str(std::vector<long int> const &x,
                                 size_t width, size_t prec)
{
  return _to_str_vector(x, width, prec);
}

std::string colvarmodule::to_str(std::vector<cvm::real> const &x,
                                 size_t width, size_t prec)
{
  return _to_str_vector(x, width, prec);
}

std::string colvarmodule::to_str(std::vector<cvm::rvector> const &x,
                                 size_t width, size_t prec)
{
  return _to_str_vector(x, width, prec);
}

std::string colvarmodule::to_str(std::vector<cvm::quaternion> const &x,
                                 size_t width, size_t prec)
{
  return _to_str_vector(x, width, prec);
}

std::string colvarmodule::to_str(std::vector<colvarvalue> const &x,
                                 size_t width, size_t prec)
{
  return _to_str_vector(x, width, prec);
}

std::string colvarmodule::to_str(std::vector<std::string> const &x,
                                 size_t width, size_t prec)
{
  return _to_str_vector(x, width, prec);
}

#if ( defined(COLVARS_CUDA) || defined(COLVARS_HIP) )
std::string colvarmodule::to_str(std::vector<cvm::real, colvars_gpu::CudaHostAllocator<cvm::real>> const &x,
                                 size_t width, size_t prec)
{
  return _to_str_vector(x, width, prec);
}

std::string colvarmodule::to_str(std::vector<cvm::rvector, colvars_gpu::CudaHostAllocator<cvm::rvector>> const &x,
                                 size_t width, size_t prec)
{
  return _to_str_vector(x, width, prec);
}
#endif


std::string colvarmodule::wrap_string(std::string const &s, size_t nchars)
{
  if (!s.size()) {
    return std::string(nchars, ' ');
  } else {
    return ( (s.size() <= nchars) ?
             (s+std::string(nchars-s.size(), ' ')) :
             (std::string(s, 0, nchars)) );
  }
}



int colvarmodule::cite_feature(std::string const &feature)
{
  return usage_->cite_feature(feature);
}

std::string colvarmodule::feature_report(int flag)
{
  return usage_->report(flag);
}


colvarmodule::usage::usage(colvarmodule *cvmodule_in)
  : cvmodule(cvmodule_in)
{
#include "colvarmodule_refs.h"
}

int colvarmodule::usage::cite_feature(std::string const &feature)
{
  if (feature_count_.count(feature) > 0) {
    feature_count_[feature] += 1;
    return cite_paper(feature_paper_map_[feature]);
  }
  cvmodule->log("Warning: cannot cite unknown feature \""+feature+"\"\n");
  return COLVARS_OK;
}

int colvarmodule::usage::cite_paper(std::string const &paper)
{
  if (paper_count_.count(paper) > 0) {
    paper_count_[paper] += 1;
    return COLVARS_OK;
  }
  cvmodule->log("Warning: cannot cite unknown paper \""+paper+"\"\n");
  return COLVARS_OK;
}

std::string colvarmodule::usage::report(int flag)
{
  std::string result;
  if (flag == 0) {
    // Text
    result += "SUMMARY OF COLVARS FEATURES USED SO FAR AND THEIR CITATIONS:\n";
  }
  if (flag == 1) {
    // LAMMPS log friendly (one-line summary, lowercase message)
    result += "Colvars module (Fiorin2013, plus other works listed for specific features)\n\n";
  }

  std::map<std::string, int>::iterator p_iter = paper_count_.begin();
  for ( ; p_iter != paper_count_.end(); p_iter++) {
    std::string const paper = p_iter->first;
    int const count = p_iter->second;
    if (count > 0) {
      result += "\n";
      std::map<std::string, std::string>::iterator f_iter =
        feature_paper_map_.begin();
      for ( ; f_iter != feature_paper_map_.end(); f_iter++) {
        if ((f_iter->second == paper) &&
            (feature_count_[f_iter->first] > 0)) {
          if (flag == 0) {
            // URL
            result += "- " + f_iter->first + ":\n";
          }
          if (flag == 1) {
            // BibTeX
            result += "% " + f_iter->first + ":\n";
          }
        }
      }
      if (flag == 0) {
        result += "  " + paper + " " + paper_url_[paper] + "\n";
      }
      if (flag == 1) {
        result += paper_bibtex_[paper] + "\n";
      }
    }
  }

  return result;
}

// Static pointer to the proxy object
colvarproxy *colvarmodule::proxy_static = nullptr;

cvm::real colvarmodule::debug_gradients_step_size = 1.0e-07;
