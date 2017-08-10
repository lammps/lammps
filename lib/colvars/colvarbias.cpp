// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarbias.h"


colvarbias::colvarbias(char const *key)
  : bias_type(to_lower_cppstr(key))
{
  init_cvb_requires();

  rank = 1;

  has_data = false;
  b_output_energy = false;
  reset();
  state_file_step = 0;
  description = "uninitialized " + cvm::to_str(key) + " bias";
}


int colvarbias::init(std::string const &conf)
{
  colvarparse::init(conf);

  if (name.size() == 0) {

    // first initialization

    cvm::log("Initializing a new \""+bias_type+"\" instance.\n");
    rank = cvm::num_biases_type(bias_type);
    get_keyval(conf, "name", name, bias_type+cvm::to_str(rank));

    {
      colvarbias *bias_with_name = cvm::bias_by_name(this->name);
      if (bias_with_name != NULL) {
        if ((bias_with_name->rank != this->rank) ||
            (bias_with_name->bias_type != this->bias_type)) {
          cvm::error("Error: this bias cannot have the same name, \""+this->name+
                     "\", as another bias.\n", INPUT_ERROR);
          return INPUT_ERROR;
        }
      }
    }

    description = "bias " + name;

    {
      // lookup the associated colvars
      std::vector<std::string> colvar_names;
      if (get_keyval(conf, "colvars", colvar_names)) {
        if (num_variables()) {
          cvm::error("Error: cannot redefine the colvars that a bias was already defined on.\n",
                     INPUT_ERROR);
          return INPUT_ERROR;
        }
        for (size_t i = 0; i < colvar_names.size(); i++) {
          add_colvar(colvar_names[i]);
        }
      }
    }

    if (!num_variables()) {
      cvm::error("Error: no collective variables specified.\n", INPUT_ERROR);
      return INPUT_ERROR;
    }
  } else {
    cvm::log("Reinitializing bias \""+name+"\".\n");
  }

  output_prefix = cvm::output_prefix();

  get_keyval(conf, "outputEnergy", b_output_energy, b_output_energy);

  get_keyval(conf, "timeStepFactor", time_step_factor, 1);
  if (time_step_factor < 1) {
    cvm::error("Error: timeStepFactor must be 1 or greater.\n");
    return COLVARS_ERROR;
  }

  // Now that children are defined, we can solve dependencies
  enable(f_cvb_active);
  if (cvm::debug()) print_state();

  return COLVARS_OK;
}


int colvarbias::reset()
{
  bias_energy = 0.0;
  for (size_t i = 0; i < num_variables(); i++) {
    colvar_forces[i].reset();
  }
  return COLVARS_OK;
}


colvarbias::colvarbias()
  : colvarparse(), has_data(false)
{}


colvarbias::~colvarbias()
{
  colvarbias::clear();
}


int colvarbias::clear()
{
  free_children_deps();

  // Remove references to this bias from colvars
  for (std::vector<colvar *>::iterator cvi = colvars.begin();
       cvi != colvars.end();
       ++cvi) {
    for (std::vector<colvarbias *>::iterator bi = (*cvi)->biases.begin();
         bi != (*cvi)->biases.end();
         ++bi) {
      if ( *bi == this) {
        (*cvi)->biases.erase(bi);
        break;
      }
    }
  }

  colvarmodule *cv = cvm::main();
  // ...and from the colvars module
  for (std::vector<colvarbias *>::iterator bi = cv->biases.begin();
       bi != cv->biases.end();
       ++bi) {
    if ( *bi == this) {
      cv->biases.erase(bi);
      break;
    }
  }

  return COLVARS_OK;
}


int colvarbias::add_colvar(std::string const &cv_name)
{
  if (colvar *cv = cvm::colvar_by_name(cv_name)) {

    if (cvm::debug()) {
      cvm::log("Applying this bias to collective variable \""+
               cv->name+"\".\n");
    }

    colvars.push_back(cv);

    colvar_forces.push_back(colvarvalue());
    colvar_forces.back().type(cv->value()); // make sure each force is initialized to zero
    colvar_forces.back().is_derivative(); // colvar constraints are not applied to the force
    colvar_forces.back().reset();

    cv->biases.push_back(this); // add back-reference to this bias to colvar

    if (is_enabled(f_cvb_apply_force)) {
      cv->enable(f_cv_gradient);
    }

    // Add dependency link.
    // All biases need at least the value of each colvar
    // although possibly not at all timesteps
    add_child(cv);

  } else {
    cvm::error("Error: cannot find a colvar named \""+
               cv_name+"\".\n", INPUT_ERROR);
    return INPUT_ERROR;
  }

  return COLVARS_OK;
}


int colvarbias::update()
{
  if (cvm::debug()) {
    cvm::log("Updating the "+bias_type+" bias \""+this->name+"\".\n");
  }

  has_data = true;

  bias_energy = 0.0;
  for (size_t ir = 0; ir < num_variables(); ir++) {
    colvar_forces[ir].reset();
  }

  return COLVARS_OK;
}


void colvarbias::communicate_forces()
{
  for (size_t i = 0; i < num_variables(); i++) {
    if (cvm::debug()) {
      cvm::log("Communicating a force to colvar \""+
               variables(i)->name+"\".\n");
    }
    // Impulse-style multiple timestep
    // Note that biases with different values of time_step_factor
    // may send forces to the same colvar
    // which is why rescaling has to happen now: the colvar is not
    // aware of this bias' time_step_factor
    variables(i)->add_bias_force(cvm::real(time_step_factor) * colvar_forces[i]);
  }
}


int colvarbias::change_configuration(std::string const &conf)
{
  cvm::error("Error: change_configuration() not implemented.\n",
             COLVARS_NOT_IMPLEMENTED);
  return COLVARS_NOT_IMPLEMENTED;
}


cvm::real colvarbias::energy_difference(std::string const &conf)
{
  cvm::error("Error: energy_difference() not implemented.\n",
             COLVARS_NOT_IMPLEMENTED);
  return 0.0;
}


// So far, these are only implemented in colvarbias_abf
int colvarbias::bin_num()
{
  cvm::error("Error: bin_num() not implemented.\n");
  return COLVARS_NOT_IMPLEMENTED;
}
int colvarbias::current_bin()
{
  cvm::error("Error: current_bin() not implemented.\n");
  return COLVARS_NOT_IMPLEMENTED;
}
int colvarbias::bin_count(int bin_index)
{
  cvm::error("Error: bin_count() not implemented.\n");
  return COLVARS_NOT_IMPLEMENTED;
}
int colvarbias::replica_share()
{
  cvm::error("Error: replica_share() not implemented.\n");
  return COLVARS_NOT_IMPLEMENTED;
}


std::string const colvarbias::get_state_params() const
{
  std::ostringstream os;
  os << "step " << cvm::step_absolute() << "\n"
     << "name " << this->name << "\n";
  return os.str();
}


int colvarbias::set_state_params(std::string const &conf)
{
  std::string new_name = "";
  if (colvarparse::get_keyval(conf, "name", new_name,
                              std::string(""), colvarparse::parse_silent) &&
      (new_name != this->name)) {
    cvm::error("Error: in the state file, the "
               "\""+bias_type+"\" block has a different name, \""+new_name+
               "\": different system?\n", INPUT_ERROR);
  }

  if (name.size() == 0) {
    cvm::error("Error: \""+bias_type+"\" block within the restart file "
               "has no identifiers.\n", INPUT_ERROR);
  }

  colvarparse::get_keyval(conf, "step", state_file_step,
                          cvm::step_absolute(), colvarparse::parse_silent);

  return COLVARS_OK;
}


std::ostream & colvarbias::write_state(std::ostream &os)
{
  if (cvm::debug()) {
    cvm::log("Writing state file for bias \""+name+"\"\n");
  }
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(cvm::cv_prec);
  os << bias_type << " {\n"
     << "  configuration {\n";
  std::istringstream is(get_state_params());
  std::string line;
  while (std::getline(is, line)) {
    os << "    " << line << "\n";
  }
  os << "  }\n";
  write_state_data(os);
  os << "}\n\n";
  return os;
}


std::istream & colvarbias::read_state(std::istream &is)
{
  size_t const start_pos = is.tellg();

  std::string key, brace, conf;
  if ( !(is >> key)   || !(key == bias_type) ||
       !(is >> brace) || !(brace == "{") ||
       !(is >> colvarparse::read_block("configuration", conf)) ||
       (set_state_params(conf) != COLVARS_OK) ) {
    cvm::error("Error: in reading state configuration for \""+bias_type+"\" bias \""+
               this->name+"\" at position "+
               cvm::to_str(is.tellg())+" in stream.\n", INPUT_ERROR);
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
    return is;
  }

  if (!read_state_data(is)) {
    cvm::error("Error: in reading state data for \""+bias_type+"\" bias \""+
               this->name+"\" at position "+
               cvm::to_str(is.tellg())+" in stream.\n", INPUT_ERROR);
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
  }

  is >> brace;
  if (brace != "}") {
    cvm::error("Error: corrupt restart information for \""+bias_type+"\" bias \""+
               this->name+"\": no matching brace at position "+
               cvm::to_str(is.tellg())+" in stream.\n");
    is.setstate(std::ios::failbit);
  }

  return is;
}


std::istream & colvarbias::read_state_data_key(std::istream &is, char const *key)
{
  size_t const start_pos = is.tellg();
  std::string key_in;
  if ( !(is >> key_in) ||
       !(key_in == to_lower_cppstr(std::string(key))) ) {
    cvm::error("Error: in reading restart configuration for "+
               bias_type+" bias \""+this->name+"\" at position "+
               cvm::to_str(is.tellg())+" in stream.\n", INPUT_ERROR);
    is.clear();
    is.seekg(start_pos, std::ios::beg);
    is.setstate(std::ios::failbit);
    return is;
  }
  return is;
}



std::ostream & colvarbias::write_traj_label(std::ostream &os)
{
  os << " ";
  if (b_output_energy)
    os << " E_"
       << cvm::wrap_string(this->name, cvm::en_width-2);
  return os;
}


std::ostream & colvarbias::write_traj(std::ostream &os)
{
  os << " ";
  if (b_output_energy)
    os << " "
       << bias_energy;
  return os;
}

// Static members

std::vector<colvardeps::feature *> colvarbias::cvb_features;
