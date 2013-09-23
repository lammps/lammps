#include "colvarmodule.h"
#include "colvarparse.h"
#include "colvarproxy.h"
#include "colvar.h"
#include "colvarbias.h"
#include "colvarbias_meta.h"
#include "colvarbias_abf.h"


colvarmodule::colvarmodule (char const  *config_filename,
                            colvarproxy *proxy_in)
{
  // pointer to the proxy object
  if (proxy == NULL) {
    proxy = proxy_in;
    parse = new colvarparse();
  } else {
    cvm::fatal_error ("Error: trying to allocate the collective "
                      "variable module twice.\n");
  }

  cvm::log (cvm::line_marker);
  cvm::log ("Initializing the collective variables module, version "+
            cvm::to_str(COLVARS_VERSION)+".\n");

  // "it_restart" will be set by the input state file, if any;
  // "it" should be updated by the proxy
  it = it_restart = 0;
  it_restart_from_state_file = true;

  // open the configfile
  config_s.open (config_filename);
  if (!config_s)
    cvm::fatal_error ("Error: in opening configuration file \""+
                      std::string (config_filename)+"\".\n");

  // read the config file into a string
  std::string conf = "";
  {
    std::string line;
    while (colvarparse::getline_nocomments (config_s, line))
      conf.append (line+"\n");
    // don't need the stream any more
    config_s.close();
  }

  std::string index_file_name;
  if (parse->get_keyval (conf, "indexFile", index_file_name)) {
    read_index_file (index_file_name.c_str());
  }

  parse->get_keyval (conf, "analysis", b_analysis, false);

  parse->get_keyval (conf, "debugGradientsStepSize", debug_gradients_step_size, 1.0e-07,
                     colvarparse::parse_silent);

  parse->get_keyval (conf, "eigenvalueCrossingThreshold",
                     colvarmodule::rotation::crossing_threshold, 1.0e-02,
                     colvarparse::parse_silent);

  parse->get_keyval (conf, "colvarsTrajFrequency", cv_traj_freq, 100);
  parse->get_keyval (conf, "colvarsRestartFrequency", restart_out_freq,
                     proxy->restart_frequency());

  // by default overwrite the existing trajectory file
  parse->get_keyval (conf, "colvarsTrajAppend", cv_traj_append, false);

  // input restart file
  restart_in_name = proxy->input_prefix().size() ?
    std::string (proxy->input_prefix()+".colvars.state") :
    std::string ("") ;

  // output restart file
  restart_out_name = proxy->restart_output_prefix().size() ?
    std::string (proxy->restart_output_prefix()+".colvars.state") :
    std::string ("");

  if (restart_out_name.size())
    cvm::log ("The restart output state file will be \""+restart_out_name+"\".\n");

  output_prefix = proxy->output_prefix();

  cvm::log ("The final output state file will be \""+
            (output_prefix.size() ?
             std::string (output_prefix+".colvars.state") :
             std::string ("colvars.state"))+"\".\n");

  cv_traj_name =
    (output_prefix.size() ?
     std::string (output_prefix+".colvars.traj") :
     std::string ("colvars.traj"));

  if (cv_traj_freq) {
    // open trajectory file
    if (cv_traj_append) {
      cvm::log ("Appending to colvar trajectory file \""+cv_traj_name+
                "\".\n");
      cv_traj_os.open (cv_traj_name.c_str(), std::ios::app);
    } else {
      cvm::log ("Writing to colvar trajectory file \""+cv_traj_name+
                "\".\n");
      proxy->backup_file (cv_traj_name.c_str());
      cv_traj_os.open (cv_traj_name.c_str(), std::ios::out);
    }
    cv_traj_os.setf (std::ios::scientific, std::ios::floatfield);
  }

  // parse the options for collective variables
  init_colvars (conf);

  // parse the options for biases
  init_biases (conf);

  // done with the parsing, check that all keywords are valid
  parse->check_keywords (conf, "colvarmodule");
  cvm::log (cvm::line_marker);

  // read the restart configuration, if available
  if (restart_in_name.size()) {
    // read the restart file
    std::ifstream input_is (restart_in_name.c_str());
    if (!input_is.good())
      fatal_error ("Error: in opening restart file \""+
                   std::string (restart_in_name)+"\".\n");
    else {
      cvm::log ("Restarting from file \""+restart_in_name+"\".\n");
      read_restart (input_is);
      cvm::log (cvm::line_marker);
    }
  }

  // check if it is possible to save output configuration
  if ((!output_prefix.size()) && (!restart_out_name.size())) {
    cvm::fatal_error ("Error: neither the final output state file or "
                      "the output restart file could be defined, exiting.\n");
  }

  cvm::log ("Collective variables module initialized.\n");
  cvm::log (cvm::line_marker);
}


std::istream & colvarmodule::read_restart (std::istream &is)
{
  {
    // read global restart information
    std::string restart_conf;
    if (is >> colvarparse::read_block ("configuration", restart_conf)) {
      if (it_restart_from_state_file) {
        parse->get_keyval (restart_conf, "step",
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
    if ( !((*cvi)->read_restart (is)) )
      cvm::fatal_error ("Error: in reading restart configuration for collective variable \""+
                        (*cvi)->name+"\".\n");
  }

  // biases restart
  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    if (!((*bi)->read_restart (is)))
      fatal_error ("Error: in reading restart configuration for bias \""+
                   (*bi)->name+"\".\n");
  }
  cvm::decrease_depth();

  return is;
}



void colvarmodule::init_colvars (std::string const &conf)
{
  if (cvm::debug())
    cvm::log ("Initializing the collective variables.\n");

  std::string colvar_conf = "";
  size_t pos = 0;
  while (parse->key_lookup (conf, "colvar", colvar_conf, pos)) {

    if (colvar_conf.size()) {
      cvm::log (cvm::line_marker);
      cvm::increase_depth();
      colvars.push_back (new colvar (colvar_conf));
      (colvars.back())->check_keywords (colvar_conf, "colvar");
      cvm::decrease_depth();
    } else {
      cvm::log ("Warning: \"colvar\" keyword found without any configuration.\n");
    }
    colvar_conf = "";
  }


  if (!colvars.size())
    cvm::fatal_error ("Error: no collective variables defined.\n");

  if (colvars.size())
    cvm::log (cvm::line_marker);
  cvm::log ("Collective variables initialized, "+
            cvm::to_str (colvars.size())+
            " in total.\n");
}


void colvarmodule::init_biases (std::string const &conf)
{
  if (cvm::debug())
    cvm::log ("Initializing the collective variables biases.\n");

  {
    /// initialize ABF instances
    std::string abf_conf = "";
    size_t abf_pos = 0;
    while (parse->key_lookup (conf, "abf", abf_conf, abf_pos)) {
      if (abf_conf.size()) {
        cvm::log (cvm::line_marker);
        cvm::increase_depth();
        biases.push_back (new colvarbias_abf (abf_conf, "abf"));
        (biases.back())->check_keywords (abf_conf, "abf");
        cvm::decrease_depth();
        n_abf_biases++;
      } else {
        cvm::log ("Warning: \"abf\" keyword found without configuration.\n");
      }
      abf_conf = "";
    }
  }

  {
    /// initialize harmonic restraints
    std::string harm_conf = "";
    size_t harm_pos = 0;
    while (parse->key_lookup (conf, "harmonic", harm_conf, harm_pos)) {
      if (harm_conf.size()) {
        cvm::log (cvm::line_marker);
        cvm::increase_depth();
        biases.push_back (new colvarbias_harmonic (harm_conf, "harmonic"));
        (biases.back())->check_keywords (harm_conf, "harmonic");
        cvm::decrease_depth();
        n_harm_biases++;
      } else {
        cvm::log ("Warning: \"harmonic\" keyword found without configuration.\n");
      }
      harm_conf = "";
    }
  }

  {
    /// initialize histograms
    std::string histo_conf = "";
    size_t histo_pos = 0;
    while (parse->key_lookup (conf, "histogram", histo_conf, histo_pos)) {
      if (histo_conf.size()) {
        cvm::log (cvm::line_marker);
        cvm::increase_depth();
        biases.push_back (new colvarbias_histogram (histo_conf, "histogram"));
        (biases.back())->check_keywords (histo_conf, "histogram");
        cvm::decrease_depth();
        n_histo_biases++;
      } else {
        cvm::log ("Warning: \"histogram\" keyword found without configuration.\n");
      }
      histo_conf = "";
    }
  }

  {
    /// initialize metadynamics instances
    std::string meta_conf = "";
    size_t meta_pos = 0;
    while (parse->key_lookup (conf, "metadynamics", meta_conf, meta_pos)) {
      if (meta_conf.size()) {
        cvm::log (cvm::line_marker);
        cvm::increase_depth();
        biases.push_back (new colvarbias_meta (meta_conf, "metadynamics"));
        (biases.back())->check_keywords (meta_conf, "metadynamics");
        cvm::decrease_depth();
        n_meta_biases++;
      } else {
        cvm::log ("Warning: \"metadynamics\" keyword found without configuration.\n");
      }
      meta_conf = "";
    }
  }

  if (biases.size())
    cvm::log (cvm::line_marker);
  cvm::log ("Collective variables biases initialized, "+
            cvm::to_str (biases.size())+" in total.\n");
}


void colvarmodule::change_configuration (std::string const &bias_name,
                                         std::string const &conf)
{
  cvm::increase_depth();
  int found = 0;
  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    if ( (*bi)->name == bias_name ) {
      ++found;
      (*bi)->change_configuration (conf);
    }
  }
  if (found < 1) cvm::fatal_error ("Error: bias not found.\n");
  if (found > 1) cvm::fatal_error ("Error: duplicate bias name.\n");
  cvm::decrease_depth();
}

cvm::real colvarmodule::energy_difference (std::string const &bias_name,
                                           std::string const &conf)
{
  cvm::increase_depth();
  cvm::real energy_diff = 0.;
  int found = 0;
  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    if ( (*bi)->name == bias_name ) {
      ++found;
      energy_diff = (*bi)->energy_difference (conf);
    }
  }
  if (found < 1) cvm::fatal_error ("Error: bias not found.\n");
  if (found > 1) cvm::fatal_error ("Error: duplicate bias name.\n");
  cvm::decrease_depth();
  return energy_diff;
}


void colvarmodule::calc() {
  cvm::real total_bias_energy = 0.0;
  cvm::real total_colvar_energy = 0.0;

  if (cvm::debug()) {
    cvm::log (cvm::line_marker);
    cvm::log ("Collective variables module, step no. "+
              cvm::to_str (cvm::step_absolute())+"\n");
  }

  // calculate collective variables and their gradients
  if (cvm::debug())
    cvm::log ("Calculating collective variables.\n");
  cvm::increase_depth();
  for (std::vector<colvar *>::iterator cvi = colvars.begin();
       cvi != colvars.end();
       cvi++) {
    (*cvi)->calc();
  }
  cvm::decrease_depth();

  // update the biases and communicate their forces to the collective
  // variables
  if (cvm::debug() && biases.size())
    cvm::log ("Updating collective variable biases.\n");
  cvm::increase_depth();
  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    total_bias_energy += (*bi)->update();
  }
  cvm::decrease_depth();

  // sum the forces from all biases for each collective variable
  if (cvm::debug() && biases.size())
    cvm::log ("Collecting forces from all biases.\n");
  cvm::increase_depth();
  for (std::vector<colvar *>::iterator cvi = colvars.begin();
       cvi != colvars.end();
       cvi++) {
    (*cvi)->reset_bias_force();
  }
  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    (*bi)->communicate_forces();
  }
  cvm::decrease_depth();

  if (cvm::b_analysis) {
    // perform runtime analysis of colvars and biases
    if (cvm::debug() && biases.size())
      cvm::log ("Perform runtime analyses.\n");
    cvm::increase_depth();
    for (std::vector<colvar *>::iterator cvi = colvars.begin();
         cvi != colvars.end();
         cvi++) {
      (*cvi)->analyse();
    }
    for (std::vector<colvarbias *>::iterator bi = biases.begin();
         bi != biases.end();
         bi++) {
      (*bi)->analyse();
    }
    cvm::decrease_depth();
  }

  // sum up the forces for each colvar and integrate any internal
  // equation of motion
  if (cvm::debug())
    cvm::log ("Updating the internal degrees of freedom "
              "of colvars (if they have any).\n");
  cvm::increase_depth();
  for (std::vector<colvar *>::iterator cvi = colvars.begin();
       cvi != colvars.end();
       cvi++) {
    total_colvar_energy += (*cvi)->update();
  }
  cvm::decrease_depth();
  proxy->add_energy (total_bias_energy + total_colvar_energy);

  // make collective variables communicate their forces to their
  // coupled degrees of freedom (i.e. atoms)
  if (cvm::debug())
    cvm::log ("Communicating forces from the colvars to the atoms.\n");
  cvm::increase_depth();
  for (std::vector<colvar *>::iterator cvi = colvars.begin();
       cvi != colvars.end();
       cvi++) {
    if ((*cvi)->tasks[colvar::task_gradients])
      (*cvi)->communicate_forces();
  }
  cvm::decrease_depth();

  // write restart file, if needed
  if (restart_out_freq && restart_out_name.size() && !cvm::b_analysis) {
    if ( (cvm::step_relative() > 0) &&
         ((cvm::step_absolute() % restart_out_freq) == 0) ) {
      cvm::log ("Writing the state file \""+
                restart_out_name+"\".\n");
      proxy->backup_file (restart_out_name.c_str());
      restart_out_os.open (restart_out_name.c_str());
      restart_out_os.setf (std::ios::scientific, std::ios::floatfield);
      if (!write_restart (restart_out_os))
        cvm::fatal_error ("Error: in writing restart file.\n");
      restart_out_os.close();
    }
  }

  // write trajectory file, if needed
  if (cv_traj_freq) {

    if (cvm::debug())
      cvm::log ("Writing trajectory file.\n");

    // (re)open trajectory file
    if (!cv_traj_os.good()) {
      if (cv_traj_append) {
        cvm::log ("Appending to colvar trajectory file \""+cv_traj_name+
                  "\".\n");
        cv_traj_os.open (cv_traj_name.c_str(), std::ios::app);
      } else {
        cvm::log ("Overwriting colvar trajectory file \""+cv_traj_name+
                  "\".\n");
        proxy->backup_file (cv_traj_name.c_str());
        cv_traj_os.open (cv_traj_name.c_str(), std::ios::out);
      }
      cv_traj_os.setf (std::ios::scientific, std::ios::floatfield);
    }

    // write labels in the traj file every 1000 lines and at first ts
    cvm::increase_depth();
    if ((cvm::step_absolute() % (cv_traj_freq * 1000)) == 0 || cvm::step_relative() == 0) {
      cv_traj_os << "# " << cvm::wrap_string ("step", cvm::it_width-2)
                 << " ";
      if (cvm::debug())
        cv_traj_os.flush();
      for (std::vector<colvar *>::iterator cvi = colvars.begin();
	   cvi != colvars.end();
	   cvi++) {
	(*cvi)->write_traj_label (cv_traj_os);
      }
      for (std::vector<colvarbias *>::iterator bi = biases.begin();
	   bi != biases.end();
	   bi++) {
	(*bi)->write_traj_label (cv_traj_os);
      }
      cv_traj_os << "\n";
      if (cvm::debug())
        cv_traj_os.flush();
    }
    cvm::decrease_depth();

    // write collective variable values to the traj file
    cvm::increase_depth();
    if ((cvm::step_absolute() % cv_traj_freq) == 0) {
      cv_traj_os << std::setw (cvm::it_width) << it
                 << " ";
      for (std::vector<colvar *>::iterator cvi = colvars.begin();
	   cvi != colvars.end();
	   cvi++) {
        (*cvi)->write_traj (cv_traj_os);
      }
      for (std::vector<colvarbias *>::iterator bi = biases.begin();
	   bi != biases.end();
	   bi++) {
	(*bi)->write_traj (cv_traj_os);
      }
      cv_traj_os << "\n";
      if (cvm::debug())
        cv_traj_os.flush();
    }
    cvm::decrease_depth();

    if (restart_out_freq) {
      // flush the trajectory file if we are at the restart frequency
      if ( (cvm::step_relative() > 0) &&
           ((cvm::step_absolute() % restart_out_freq) == 0) ) {
        cvm::log ("Synchronizing (emptying the buffer of) trajectory file \""+
                  cv_traj_name+"\".\n");
        cv_traj_os.flush();
      }
    }
  } // end if (cv_traj_freq)
}


void colvarmodule::analyze()
{
  if (cvm::debug()) {
    cvm::log ("colvarmodule::analyze(), step = "+cvm::to_str (it)+".\n");
  }

  if (cvm::step_relative() == 0)
    cvm::log ("Performing analysis.\n");

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

}

void colvarmodule::setup()
{
  // loop over all components of all colvars to reset masses of all groups
  for (std::vector<colvar *>::iterator cvi = colvars.begin();
       cvi != colvars.end();  cvi++) {
    (*cvi)->setup();
  }
}

colvarmodule::~colvarmodule()
{
  for (std::vector<colvar *>::iterator cvi = colvars.begin();
       cvi != colvars.end();
       cvi++) {
    delete *cvi;
  }
  colvars.clear();

  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    delete *bi;
  }
  biases.clear();

  if (cv_traj_os.good()) {
    cv_traj_os.close();
  }

  delete parse;
  proxy = NULL;
}


void colvarmodule::write_output_files()
{
  // if this is a simulation run (i.e. not a postprocessing), output data
  // must be written to be able to restart the simulation
  std::string const out_name =
    (output_prefix.size() ?
     std::string (output_prefix+".colvars.state") :
     std::string ("colvars.state"));
  cvm::log ("Saving collective variables state to \""+out_name+"\".\n");
  proxy->backup_file (out_name.c_str());
  std::ofstream out (out_name.c_str());
  out.setf (std::ios::scientific, std::ios::floatfield);
  this->write_restart (out);
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
}



bool colvarmodule::read_traj (char const *traj_filename,
                              size_t      traj_read_begin,
                              size_t      traj_read_end)
{
  cvm::log ("Opening trajectory file \""+
            std::string (traj_filename)+"\".\n");
  std::ifstream traj_is (traj_filename);

  while (true) {
    while (true) {

      std::string line ("");

      do {
        if (!colvarparse::getline_nocomments (traj_is, line)) {
          cvm::log ("End of file \""+std::string (traj_filename)+
                    "\" reached, or corrupted file.\n");
          traj_is.close();
          return false;
        }
      } while (line.find_first_not_of (colvarparse::white_space) == std::string::npos);

      std::istringstream is (line);

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
          cvm::log ("Reached the end of the trajectory, "
                    "read_end = "+cvm::to_str (traj_read_end)+"\n");
          return false;
        }

        for (std::vector<colvar *>::iterator cvi = colvars.begin();
             cvi != colvars.end();
             cvi++) {
          if (!(*cvi)->read_traj (is)) {
            cvm::log ("Error: in reading colvar \""+(*cvi)->name+
                      "\" from trajectory file \""+
                      std::string (traj_filename)+"\".\n");
            return false;
          }
        }

        break;
      }
    }
  }

  return true;
}


std::ostream & colvarmodule::write_restart (std::ostream &os)
{
  os << "configuration {\n"
     << "  step " << std::setw (it_width)
     << it << "\n"
     << "  dt " << dt() << "\n"
     << "}\n\n";

  cvm::increase_depth();
  for (std::vector<colvar *>::iterator cvi = colvars.begin();
       cvi != colvars.end();
       cvi++) {
    (*cvi)->write_restart (os);
  }

  for (std::vector<colvarbias *>::iterator bi = biases.begin();
       bi != biases.end();
       bi++) {
    (*bi)->write_restart (os);
  }
  cvm::decrease_depth();

  return os;
}



void cvm::log (std::string const &message)
{
  if (depth > 0)
    proxy->log ((std::string (2*depth, ' '))+message);
  else
    proxy->log (message);
}

void cvm::increase_depth()
{
  depth++;
}

void cvm::decrease_depth()
{
  if (depth) depth--;
}

void cvm::fatal_error (std::string const &message)
{
  proxy->fatal_error (message);
}

void cvm::exit (std::string const &message)
{
  proxy->exit (message);
}


void cvm::read_index_file (char const *filename)
{
  std::ifstream is (filename, std::ios::binary);
  if (!is.good())
    fatal_error ("Error: in opening index file \""+
                 std::string (filename)+"\".\n");
  // std::list<std::string>::iterator names_i = cvm::index_group_names.begin();
  // std::list<std::vector<int> >::iterator lists_i = cvm::index_groups.begin();
  while (is.good()) {
    char open, close;
    std::string group_name;
    if ( (is >> open) && (open == '[') &&
         (is >> group_name) &&
         (is >> close) && (close == ']') ) {
      cvm::index_group_names.push_back (group_name);
      cvm::index_groups.push_back (std::vector<int> ());
    } else {
      cvm::fatal_error ("Error: in parsing index file \""+
                 std::string (filename)+"\".\n");
    }

    int atom_number = 1;
    size_t pos = is.tellg();
    while ( (is >> atom_number) && (atom_number > 0) ) {
      (cvm::index_groups.back()).push_back (atom_number);
      pos = is.tellg();
    }
    is.clear();
    is.seekg (pos, std::ios::beg);
    std::string delim;
    if ( (is >> delim) && (delim == "[") ) {
      // new group
      is.clear();
      is.seekg (pos, std::ios::beg);
    } else {
      break;
    }
  }

  cvm::log ("The following index groups were read from the index file \""+
            std::string (filename)+"\":\n");
  std::list<std::string>::iterator names_i = cvm::index_group_names.begin();
  std::list<std::vector<int> >::iterator lists_i = cvm::index_groups.begin();
  for ( ; names_i != cvm::index_group_names.end() ; names_i++, lists_i++) {
    cvm::log ("  "+(*names_i)+" ("+cvm::to_str (lists_i->size())+" atoms).\n");
  }

}


// static pointers
std::vector<colvar *>     colvarmodule::colvars;
std::vector<colvarbias *> colvarmodule::biases;
size_t                    colvarmodule::n_abf_biases = 0;
size_t                    colvarmodule::n_harm_biases = 0;
size_t                    colvarmodule::n_histo_biases = 0;
size_t                    colvarmodule::n_meta_biases = 0;
colvarproxy              *colvarmodule::proxy = NULL;


// static runtime data
cvm::real colvarmodule::debug_gradients_step_size = 1.0e-03;
size_t    colvarmodule::it = 0;
size_t    colvarmodule::it_restart = 0;
size_t    colvarmodule::restart_out_freq = 0;
size_t    colvarmodule::cv_traj_freq = 0;
size_t    colvarmodule::depth = 0;
bool      colvarmodule::b_analysis = false;
cvm::real colvarmodule::rotation::crossing_threshold = 1.0E-04;
std::list<std::string> colvarmodule::index_group_names;
std::list<std::vector<int> > colvarmodule::index_groups;


// file name prefixes
std::string colvarmodule::output_prefix = "";
std::string colvarmodule::input_prefix = "";
std::string colvarmodule::restart_in_name = "";


// i/o constants
size_t const colvarmodule::it_width = 12;
size_t const colvarmodule::cv_prec  = 14;
size_t const colvarmodule::cv_width = 21;
size_t const colvarmodule::en_prec  = 14;
size_t const colvarmodule::en_width = 21;
std::string const colvarmodule::line_marker =
  "----------------------------------------------------------------------\n";






std::ostream & operator << (std::ostream &os, colvarmodule::rvector const &v)
{
  std::streamsize const w = os.width();
  std::streamsize const p = os.precision();

  os.width (2);
  os << "( ";
  os.width (w); os.precision (p);
  os << v.x << " , ";
  os.width (w); os.precision (p);
  os << v.y << " , ";
  os.width (w); os.precision (p);
  os << v.z << " )";
  return os;
}


std::istream & operator >> (std::istream &is, colvarmodule::rvector &v)
{
  size_t const start_pos = is.tellg();
  char sep;
  if ( !(is >> sep) || !(sep == '(') ||
       !(is >> v.x) || !(is >> sep)  || !(sep == ',') ||
       !(is >> v.y) || !(is >> sep)  || !(sep == ',') ||
       !(is >> v.z) || !(is >> sep)  || !(sep == ')') ) {
    is.clear();
    is.seekg (start_pos, std::ios::beg);
    is.setstate (std::ios::failbit);
    return is;
  }
  return is;
}



std::ostream & operator << (std::ostream &os, colvarmodule::quaternion const &q)
{
  std::streamsize const w = os.width();
  std::streamsize const p = os.precision();

  os.width (2);
  os << "( ";
  os.width (w); os.precision (p);
  os << q.q0 << " , ";
  os.width (w); os.precision (p);
  os << q.q1 << " , ";
  os.width (w); os.precision (p);
  os << q.q2 << " , ";
  os.width (w); os.precision (p);
  os << q.q3 << " )";
  return os;
}


std::istream & operator >> (std::istream &is, colvarmodule::quaternion &q)
{
  size_t const start_pos = is.tellg();

  std::string euler ("");

  if ( (is >> euler) && (colvarparse::to_lower_cppstr (euler) ==
                         std::string ("euler")) ) {

    // parse the Euler angles

    char sep;
    cvm::real phi, theta, psi;
    if ( !(is >> sep)   || !(sep == '(') ||
         !(is >> phi)   || !(is >> sep)  || !(sep == ',') ||
         !(is >> theta) || !(is >> sep)  || !(sep == ',') ||
         !(is >> psi)   || !(is >> sep)  || !(sep == ')') ) {
      is.clear();
      is.seekg (start_pos, std::ios::beg);
      is.setstate (std::ios::failbit);
      return is;
    }

    q = colvarmodule::quaternion (phi, theta, psi);

  } else {

    // parse the quaternion components

    is.seekg (start_pos, std::ios::beg);
    char sep;
    if ( !(is >> sep)  || !(sep == '(') ||
         !(is >> q.q0) || !(is >> sep)  || !(sep == ',') ||
         !(is >> q.q1) || !(is >> sep)  || !(sep == ',') ||
         !(is >> q.q2) || !(is >> sep)  || !(sep == ',') ||
         !(is >> q.q3) || !(is >> sep)  || !(sep == ')') ) {
      is.clear();
      is.seekg (start_pos, std::ios::beg);
      is.setstate (std::ios::failbit);
      return is;
    }
  }

  return is;
}


cvm::quaternion
cvm::quaternion::position_derivative_inner (cvm::rvector const &pos,
                                            cvm::rvector const &vec) const
{
  cvm::quaternion result (0.0, 0.0, 0.0, 0.0);


  result.q0 =   2.0 * pos.x * q0 * vec.x
               +2.0 * pos.y * q0 * vec.y
               +2.0 * pos.z * q0 * vec.z

               -2.0 * pos.y * q3 * vec.x
               +2.0 * pos.z * q2 * vec.x

               +2.0 * pos.x * q3 * vec.y
               -2.0 * pos.z * q1 * vec.y

               -2.0 * pos.x * q2 * vec.z
               +2.0 * pos.y * q1 * vec.z;


  result.q1 =  +2.0 * pos.x * q1 * vec.x
               -2.0 * pos.y * q1 * vec.y
               -2.0 * pos.z * q1 * vec.z

               +2.0 * pos.y * q2 * vec.x
               +2.0 * pos.z * q3 * vec.x

               +2.0 * pos.x * q2 * vec.y
               -2.0 * pos.z * q0 * vec.y

               +2.0 * pos.x * q3 * vec.z
               +2.0 * pos.y * q0 * vec.z;


  result.q2 =  -2.0 * pos.x * q2 * vec.x
               +2.0 * pos.y * q2 * vec.y
               -2.0 * pos.z * q2 * vec.z

               +2.0 * pos.y * q1 * vec.x
               +2.0 * pos.z * q0 * vec.x

               +2.0 * pos.x * q1 * vec.y
               +2.0 * pos.z * q3 * vec.y

               -2.0 * pos.x * q0 * vec.z
               +2.0 * pos.y * q3 * vec.z;


  result.q3 =  -2.0 * pos.x * q3 * vec.x
               -2.0 * pos.y * q3 * vec.y
               +2.0 * pos.z * q3 * vec.z

               -2.0 * pos.y * q0 * vec.x
               +2.0 * pos.z * q1 * vec.x

               +2.0 * pos.x * q0 * vec.y
               +2.0 * pos.z * q2 * vec.y

               +2.0 * pos.x * q1 * vec.z
               +2.0 * pos.y * q2 * vec.z;

  return result;
}






// Calculate the optimal rotation between two groups, and implement it
// as a quaternion.  The method is the one documented in: Coutsias EA,
// Seok C, Dill KA.  Using quaternions to calculate RMSD.  J Comput
// Chem. 25(15):1849-57 (2004) DOI: 10.1002/jcc.20110 PubMed: 15376254

void colvarmodule::rotation::build_matrix (std::vector<cvm::atom_pos> const &pos1,
                                           std::vector<cvm::atom_pos> const &pos2,
                                           matrix2d<cvm::real, 4, 4>        &S)
{
  cvm::rmatrix C;

  // build the correlation matrix
  C.reset();
  for (size_t i = 0; i < pos1.size(); i++) {
    C.xx() += pos1[i].x * pos2[i].x;
    C.xy() += pos1[i].x * pos2[i].y;
    C.xz() += pos1[i].x * pos2[i].z;
    C.yx() += pos1[i].y * pos2[i].x;
    C.yy() += pos1[i].y * pos2[i].y;
    C.yz() += pos1[i].y * pos2[i].z;
    C.zx() += pos1[i].z * pos2[i].x;
    C.zy() += pos1[i].z * pos2[i].y;
    C.zz() += pos1[i].z * pos2[i].z;
  }

  // build the "overlap" matrix, whose eigenvectors are stationary
  // points of the RMSD in the space of rotations
  S[0][0] =    C.xx() + C.yy() + C.zz();
  S[1][0] =    C.yz() - C.zy();
  S[0][1] = S[1][0];
  S[2][0] =  - C.xz() + C.zx() ;
  S[0][2] = S[2][0];
  S[3][0] =    C.xy() - C.yx();
  S[0][3] = S[3][0];
  S[1][1] =    C.xx() - C.yy() - C.zz();
  S[2][1] =    C.xy() + C.yx();
  S[1][2] = S[2][1];
  S[3][1] =    C.xz() + C.zx();
  S[1][3] = S[3][1];
  S[2][2] = - C.xx() + C.yy() - C.zz();
  S[3][2] =   C.yz() + C.zy();
  S[2][3] = S[3][2];
  S[3][3] = - C.xx() - C.yy() + C.zz();

  //   if (cvm::debug()) {
  //     for (size_t i = 0; i < 4; i++) {
  //       std::string line ("");
  //       for (size_t j = 0; j < 4; j++) {
  //         line += std::string (" S["+cvm::to_str (i)+
  //                              "]["+cvm::to_str (j)+"] ="+cvm::to_str (S[i][j]));
  //       }
  //       cvm::log (line+"\n");
  //     }
  //   }
}


void colvarmodule::rotation::diagonalize_matrix (matrix2d<cvm::real, 4, 4> &S,
                                                 cvm::real                  S_eigval[4],
                                                 matrix2d<cvm::real, 4, 4> &S_eigvec)
{
  // diagonalize
  int jac_nrot = 0;
  jacobi (S, 4, S_eigval, S_eigvec, &jac_nrot);
  eigsrt (S_eigval, S_eigvec, 4);
  // jacobi saves eigenvectors by columns
  transpose (S_eigvec, 4);

  // normalize eigenvectors
  for (size_t ie = 0; ie < 4; ie++) {
    cvm::real norm2 = 0.0;
    for (size_t i = 0; i < 4; i++) norm2 += std::pow (S_eigvec[ie][i], int (2));
    cvm::real const norm = std::sqrt (norm2);
    for (size_t i = 0; i < 4; i++) S_eigvec[ie][i] /= norm;
  }
}


// Calculate the rotation, plus its derivatives

void colvarmodule::rotation::calc_optimal_rotation
(std::vector<cvm::atom_pos> const &pos1,
 std::vector<cvm::atom_pos> const &pos2)
{
  matrix2d<cvm::real, 4, 4> S;
  matrix2d<cvm::real, 4, 4> S_backup;
  cvm::real                 S_eigval[4];
  matrix2d<cvm::real, 4, 4> S_eigvec;

//   if (cvm::debug()) {
//     cvm::atom_pos cog1 (0.0, 0.0, 0.0);
//     for (size_t i = 0; i < pos1.size(); i++) {
//       cog1 += pos1[i];
//     }
//     cog1 /= cvm::real (pos1.size());
//     cvm::atom_pos cog2 (0.0, 0.0, 0.0);
//     for (size_t i = 0; i < pos2.size(); i++) {
//       cog2 += pos2[i];
//     }
//     cog2 /= cvm::real (pos1.size());
//     cvm::log ("calc_optimal_rotation: centers of geometry are: "+
//               cvm::to_str (cog1, cvm::cv_width, cvm::cv_prec)+
//               " and "+cvm::to_str (cog2, cvm::cv_width, cvm::cv_prec)+".\n");
//   }

  build_matrix (pos1, pos2, S);
  S_backup = S;

  if (cvm::debug()) {
    if (b_debug_gradients) {
      cvm::log ("S     = "+cvm::to_str (cvm::to_str (S_backup), cvm::cv_width, cvm::cv_prec)+"\n");
    }
  }

  diagonalize_matrix (S, S_eigval, S_eigvec);

  // eigenvalues and eigenvectors
  cvm::real const &L0 = S_eigval[0];
  cvm::real const &L1 = S_eigval[1];
  cvm::real const &L2 = S_eigval[2];
  cvm::real const &L3 = S_eigval[3];
  cvm::real const *Q0 = S_eigvec[0];
  cvm::real const *Q1 = S_eigvec[1];
  cvm::real const *Q2 = S_eigvec[2];
  cvm::real const *Q3 = S_eigvec[3];

  lambda = L0;
  q = cvm::quaternion (Q0);

  if (q_old.norm2() > 0.0) {
    q.match (q_old);
    if (q_old.inner (q) < (1.0 - crossing_threshold)) {
      cvm::log ("Warning: one molecular orientation has changed by more than "+
                cvm::to_str (crossing_threshold)+": discontinuous rotation ?\n");
    }
  }
  q_old = q;

  if (cvm::debug()) {
    if (b_debug_gradients) {
      cvm::log ("L0 = "+cvm::to_str (L0, cvm::cv_width, cvm::cv_prec)+
                ", Q0 = "+cvm::to_str (cvm::quaternion (Q0), cvm::cv_width, cvm::cv_prec)+
                ", Q0*Q0 = "+cvm::to_str (cvm::quaternion (Q0).inner (cvm::quaternion (Q0)), cvm::cv_width, cvm::cv_prec)+
                "\n");
      cvm::log ("L1 = "+cvm::to_str (L1, cvm::cv_width, cvm::cv_prec)+
                ", Q1 = "+cvm::to_str (cvm::quaternion (Q1), cvm::cv_width, cvm::cv_prec)+
                ", Q0*Q1 = "+cvm::to_str (cvm::quaternion (Q0).inner (cvm::quaternion (Q1)), cvm::cv_width, cvm::cv_prec)+
                "\n");
      cvm::log ("L2 = "+cvm::to_str (L2, cvm::cv_width, cvm::cv_prec)+
                ", Q2 = "+cvm::to_str (cvm::quaternion (Q2), cvm::cv_width, cvm::cv_prec)+
                ", Q0*Q2 = "+cvm::to_str (cvm::quaternion (Q0).inner (cvm::quaternion (Q2)), cvm::cv_width, cvm::cv_prec)+
                "\n");
      cvm::log ("L3 = "+cvm::to_str (L3, cvm::cv_width, cvm::cv_prec)+
                ", Q3 = "+cvm::to_str (cvm::quaternion (Q3), cvm::cv_width, cvm::cv_prec)+
                ", Q0*Q3 = "+cvm::to_str (cvm::quaternion (Q0).inner (cvm::quaternion (Q3)), cvm::cv_width, cvm::cv_prec)+
                "\n");
    }
  }

  // calculate derivatives of L0 and Q0 with respect to each atom in
  // either group; note: if dS_1 is a null vector, nothing will be
  // calculated
  for (size_t ia = 0; ia < dS_1.size(); ia++) {

    cvm::real const &a2x = pos2[ia].x;
    cvm::real const &a2y = pos2[ia].y;
    cvm::real const &a2z = pos2[ia].z;

    matrix2d<cvm::rvector, 4, 4>    &ds_1 =  dS_1[ia];

    // derivative of the S matrix
    ds_1.reset();
    ds_1[0][0] = cvm::rvector ( a2x,  a2y,  a2z);
    ds_1[1][0] = cvm::rvector ( 0.0,  a2z, -a2y);
    ds_1[0][1] = ds_1[1][0];
    ds_1[2][0] = cvm::rvector (-a2z,  0.0,  a2x);
    ds_1[0][2] = ds_1[2][0];
    ds_1[3][0] = cvm::rvector ( a2y, -a2x,  0.0);
    ds_1[0][3] = ds_1[3][0];
    ds_1[1][1] = cvm::rvector ( a2x, -a2y, -a2z);
    ds_1[2][1] = cvm::rvector ( a2y,  a2x,  0.0);
    ds_1[1][2] = ds_1[2][1];
    ds_1[3][1] = cvm::rvector ( a2z,  0.0,  a2x);
    ds_1[1][3] = ds_1[3][1];
    ds_1[2][2] = cvm::rvector (-a2x,  a2y, -a2z);
    ds_1[3][2] = cvm::rvector ( 0.0,  a2z,  a2y);
    ds_1[2][3] = ds_1[3][2];
    ds_1[3][3] = cvm::rvector (-a2x, -a2y,  a2z);

    cvm::rvector              &dl0_1 = dL0_1[ia];
    vector1d<cvm::rvector, 4> &dq0_1 = dQ0_1[ia];

    // matrix multiplications; derivatives of L_0 and Q_0 are
    // calculated using Hellmann-Feynman theorem (i.e. exploiting the
    // fact that the eigenvectors Q_i form an orthonormal basis)

    dl0_1.reset();
    for (size_t i = 0; i < 4; i++) {
      for (size_t j = 0; j < 4; j++) {
        dl0_1 += Q0[i] * ds_1[i][j] * Q0[j];
      }
    }

    dq0_1.reset();
    for (size_t p = 0; p < 4; p++) {
      for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++) {
          dq0_1[p] +=
            (Q1[i] * ds_1[i][j] * Q0[j]) / (L0-L1) * Q1[p] +
            (Q2[i] * ds_1[i][j] * Q0[j]) / (L0-L2) * Q2[p] +
            (Q3[i] * ds_1[i][j] * Q0[j]) / (L0-L3) * Q3[p];
        }
      }
    }
  }

  // do the same for the second group
  for (size_t ia = 0; ia < dS_2.size(); ia++) {

    cvm::real const &a1x = pos1[ia].x;
    cvm::real const &a1y = pos1[ia].y;
    cvm::real const &a1z = pos1[ia].z;

    matrix2d<cvm::rvector, 4, 4> &ds_2 =  dS_2[ia];

    ds_2.reset();
    ds_2[0][0] = cvm::rvector ( a1x,  a1y,  a1z);
    ds_2[1][0] = cvm::rvector ( 0.0, -a1z,  a1y);
    ds_2[0][1] = ds_2[1][0];
    ds_2[2][0] = cvm::rvector ( a1z,  0.0, -a1x);
    ds_2[0][2] = ds_2[2][0];
    ds_2[3][0] = cvm::rvector (-a1y,  a1x,  0.0);
    ds_2[0][3] = ds_2[3][0];
    ds_2[1][1] = cvm::rvector ( a1x, -a1y, -a1z);
    ds_2[2][1] = cvm::rvector ( a1y,  a1x,  0.0);
    ds_2[1][2] = ds_2[2][1];
    ds_2[3][1] = cvm::rvector ( a1z,  0.0,  a1x);
    ds_2[1][3] = ds_2[3][1];
    ds_2[2][2] = cvm::rvector (-a1x,  a1y, -a1z);
    ds_2[3][2] = cvm::rvector ( 0.0,  a1z,  a1y);
    ds_2[2][3] = ds_2[3][2];
    ds_2[3][3] = cvm::rvector (-a1x, -a1y,  a1z);

    cvm::rvector              &dl0_2 = dL0_2[ia];
    vector1d<cvm::rvector, 4> &dq0_2 = dQ0_2[ia];

    dl0_2.reset();
    for (size_t i = 0; i < 4; i++) {
      for (size_t j = 0; j < 4; j++) {
        dl0_2 += Q0[i] * ds_2[i][j] * Q0[j];
      }
    }

    dq0_2.reset();
    for (size_t p = 0; p < 4; p++) {
      for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++) {
          dq0_2[p] +=
            (Q1[i] * ds_2[i][j] * Q0[j]) / (L0-L1) * Q1[p] +
            (Q2[i] * ds_2[i][j] * Q0[j]) / (L0-L2) * Q2[p] +
            (Q3[i] * ds_2[i][j] * Q0[j]) / (L0-L3) * Q3[p];
        }
      }
    }

    if (cvm::debug()) {

      if (b_debug_gradients) {

      matrix2d<cvm::real, 4, 4> S_new;
      cvm::real                 S_new_eigval[4];
      matrix2d<cvm::real, 4, 4> S_new_eigvec;

        // make an infitesimal move along each cartesian coordinate of
        // this atom, and solve again the eigenvector problem
        for (size_t comp = 0; comp < 3; comp++) {

          S_new = S_backup;
          // diagonalize the new overlap matrix
          for (size_t i = 0; i < 4; i++) {
            for (size_t j = 0; j < 4; j++) {
              S_new[i][j] +=
                colvarmodule::debug_gradients_step_size * ds_2[i][j][comp];
            }
          }

//           cvm::log ("S_new = "+cvm::to_str (cvm::to_str (S_new), cvm::cv_width, cvm::cv_prec)+"\n");

          diagonalize_matrix (S_new, S_new_eigval, S_new_eigvec);

          cvm::real const &L0_new = S_new_eigval[0];
          cvm::real const *Q0_new = S_new_eigvec[0];

          cvm::real const DL0 = (dl0_2[comp]) * colvarmodule::debug_gradients_step_size;
          cvm::quaternion const q0 (Q0);
          cvm::quaternion const DQ0 (dq0_2[0][comp] * colvarmodule::debug_gradients_step_size,
                                     dq0_2[1][comp] * colvarmodule::debug_gradients_step_size,
                                     dq0_2[2][comp] * colvarmodule::debug_gradients_step_size,
                                     dq0_2[3][comp] * colvarmodule::debug_gradients_step_size);

          cvm::log (  "|(l_0+dl_0) - l_0^new|/l_0 = "+
                      cvm::to_str (std::fabs (L0+DL0 - L0_new)/L0, cvm::cv_width, cvm::cv_prec)+
                      ", |(q_0+dq_0) - q_0^new| = "+
                      cvm::to_str ((Q0+DQ0 - Q0_new).norm(), cvm::cv_width, cvm::cv_prec)+
                      "\n");
        }
      }
    }
  }
}



// Numerical Recipes routine for diagonalization

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);    \
  a[k][l]=h+s*(g-h*tau);
void jacobi(cvm::real **a, int n, cvm::real d[], cvm::real **v, int *nrot)
{
  int j,iq,ip,i;
  cvm::real tresh,theta,tau,t,sm,s,h,g,c;

  std::vector<cvm::real> b (n, 0.0);
  std::vector<cvm::real> z (n, 0.0);

  for (ip=0;ip<n;ip++) {
    for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }
  for (ip=0;ip<n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=0;i<=50;i++) {
    sm=0.0;
    for (ip=0;ip<n-1;ip++) {
      for (iq=ip+1;iq<n;iq++)
        sm += std::fabs(a[ip][iq]);
    }
    if (sm == 0.0) {
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=0;ip<n-1;ip++) {
      for (iq=ip+1;iq<n;iq++) {
        g=100.0*std::fabs(a[ip][iq]);
        if (i > 4 && (cvm::real)(std::fabs(d[ip])+g) == (cvm::real)std::fabs(d[ip])
            && (cvm::real)(std::fabs(d[iq])+g) == (cvm::real)std::fabs(d[iq]))
          a[ip][iq]=0.0;
        else if (std::fabs(a[ip][iq]) > tresh) {
          h=d[iq]-d[ip];
          if ((cvm::real)(std::fabs(h)+g) == (cvm::real)std::fabs(h))
            t=(a[ip][iq])/h;
          else {
            theta=0.5*h/(a[ip][iq]);
            t=1.0/(std::fabs(theta)+std::sqrt(1.0+theta*theta));
            if (theta < 0.0) t = -t;
          }
          c=1.0/std::sqrt(1+t*t);
          s=t*c;
          tau=s/(1.0+c);
          h=t*a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[ip][iq]=0.0;
          for (j=0;j<=ip-1;j++) {
            ROTATE(a,j,ip,j,iq)
              }
          for (j=ip+1;j<=iq-1;j++) {
            ROTATE(a,ip,j,j,iq)
              }
          for (j=iq+1;j<n;j++) {
            ROTATE(a,ip,j,iq,j)
              }
          for (j=0;j<n;j++) {
            ROTATE(v,j,ip,j,iq)
              }
          ++(*nrot);
        }
      }
    }
    for (ip=0;ip<n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  cvm::fatal_error ("Too many iterations in routine jacobi.\n");
}
#undef ROTATE


void eigsrt(cvm::real d[], cvm::real **v, int n)
{
  int k,j,i;
  cvm::real p;

  for (i=0;i<n;i++) {
    p=d[k=i];
    for (j=i+1;j<n;j++)
      if (d[j] >= p) p=d[k=j];
    if (k != i) {
      d[k]=d[i];
      d[i]=p;
      for (j=0;j<n;j++) {
        p=v[j][i];
        v[j][i]=v[j][k];
        v[j][k]=p;
      }
    }
  }
}


void transpose(cvm::real **v, int n)
{
  cvm::real p;
  for (int i=0;i<n;i++) {
    for (int j=i+1;j<n;j++) {
      p=v[i][j];
      v[i][j]=v[j][i];
      v[j][i]=p;
    }
  }
}

