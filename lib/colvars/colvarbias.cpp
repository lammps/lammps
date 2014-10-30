/// -*- c++ -*-

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarbias.h"


colvarbias::colvarbias(std::string const &conf, char const *key)
  : colvarparse(), has_data(false)
{
  cvm::log("Initializing a new \""+std::string(key)+"\" instance.\n");

  size_t rank = 1;
  std::string const key_str(key);

  if (to_lower_cppstr(key_str) == std::string("abf")) {
    rank = cvm::n_abf_biases+1;
  }
  if (to_lower_cppstr(key_str) == std::string("harmonic") ||
      to_lower_cppstr(key_str) == std::string("linear")) {
    rank = cvm::n_rest_biases+1;
  }
  if (to_lower_cppstr(key_str) == std::string("histogram")) {
    rank = cvm::n_histo_biases+1;
  }
  if (to_lower_cppstr(key_str) == std::string("metadynamics")) {
    rank = cvm::n_meta_biases+1;
  }

  get_keyval(conf, "name", name, key_str+cvm::to_str(rank));

  if (cvm::bias_by_name(this->name) != NULL) {
    cvm::error("Error: this bias cannot have the same name, \""+this->name+
                "\", as another bias.\n", INPUT_ERROR);
    return;
  }

  // lookup the associated colvars
  std::vector<std::string> colvars_str;
  if (get_keyval(conf, "colvars", colvars_str)) {
    for (size_t i = 0; i < colvars_str.size(); i++) {
      add_colvar(colvars_str[i]);
    }
  }
  if (!colvars.size()) {
    cvm::error("Error: no collective variables specified.\n");
    return;
  }

  get_keyval(conf, "outputEnergy", b_output_energy, false);
}


colvarbias::colvarbias()
  : colvarparse(), has_data(false)
{}

colvarbias::~colvarbias()
{
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
  // ...and from the colvars module
  for (std::vector<colvarbias *>::iterator bi = cvm::biases.begin();
       bi != cvm::biases.end();
       ++bi) {
    if ( *bi == this) {
      cvm::biases.erase(bi);
      break;
    }
  }
}

void colvarbias::add_colvar(std::string const &cv_name)
{
  if (colvar *cv = cvm::colvar_by_name(cv_name)) {
    cv->enable(colvar::task_gradients);
    if (cvm::debug())
      cvm::log("Applying this bias to collective variable \""+
                cv->name+"\".\n");
    colvars.push_back(cv);
    colvar_forces.push_back(colvarvalue());
    colvar_forces.back().type(cv->value()); // make sure each forces is initialized to zero
    colvar_forces.back().reset();
    cv->biases.push_back(this); // add back-reference to this bias to colvar
  } else {
    cvm::error("Error: cannot find a colvar named \""+
               cv_name+"\".\n");
  }
}


void colvarbias::communicate_forces()
{
  for (size_t i = 0; i < colvars.size(); i++) {
    if (cvm::debug()) {
      cvm::log("Communicating a force to colvar \""+
               colvars[i]->name+"\".\n");
    }
    colvars[i]->add_bias_force(colvar_forces[i]);
  }
}


void colvarbias::change_configuration(std::string const &conf)
{
  cvm::error("Error: change_configuration() not implemented.\n");
}


cvm::real colvarbias::energy_difference(std::string const &conf)
{
  cvm::error("Error: energy_difference() not implemented.\n");
  return 0.;
}


// So far, these are only implemented in colvarsbias_abf
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
