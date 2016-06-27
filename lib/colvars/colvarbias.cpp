// -*- c++ -*-

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarbias.h"


colvarbias::colvarbias(char const *key)
  : bias_type(to_lower_cppstr(key))
{
  init_cvb_requires();

  rank = 1;

  if (bias_type == std::string("abf")) {
    rank = cvm::n_abf_biases+1;
  }
  if (bias_type == std::string("harmonic") ||
      bias_type == std::string("linear")) {
    rank = cvm::n_rest_biases+1;
  }
  if (bias_type == std::string("histogram")) {
    rank = cvm::n_histo_biases+1;
  }
  if (bias_type == std::string("metadynamics")) {
    rank = cvm::n_meta_biases+1;
  }

  has_data = false;
  b_output_energy = false;
  reset();

  // Start in active state by default
  enable(f_cvb_active);
}


int colvarbias::init(std::string const &conf)
{
  colvarparse::init(conf);

  if (name.size() == 0) {
    cvm::log("Initializing a new \""+bias_type+"\" instance.\n");
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
        if (colvars.size()) {
          cvm::error("Error: cannot redefine the colvars that a bias was already defined on.\n",
                     INPUT_ERROR);
          return INPUT_ERROR;
        }
        for (size_t i = 0; i < colvar_names.size(); i++) {
          add_colvar(colvar_names[i]);
        }
      }
    }

    if (!colvars.size()) {
      cvm::error("Error: no collective variables specified.\n", INPUT_ERROR);
      return INPUT_ERROR;
    }

  } else {
    cvm::log("Reinitializing bias \""+name+"\".\n");
  }

  get_keyval(conf, "outputEnergy", b_output_energy, b_output_energy);

  return COLVARS_OK;
}


int colvarbias::reset()
{
  bias_energy = 0.0;
  for (size_t i = 0; i < colvars.size(); i++) {
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

  return COLVARS_OK;
}


int colvarbias::add_colvar(std::string const &cv_name)
{
  if (colvar *cv = cvm::colvar_by_name(cv_name)) {
    // Removed this as nor all biases apply forces eg histogram
    // cv->enable(colvar::task_gradients);
    if (cvm::debug()) {
      cvm::log("Applying this bias to collective variable \""+
               cv->name+"\".\n");
    }
    colvars.push_back(cv);

    colvar_forces.push_back(colvarvalue());
    colvar_forces.back().type(cv->value()); // make sure each force is initialized to zero
    colvar_forces.back().reset();

    cv->biases.push_back(this); // add back-reference to this bias to colvar

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
  // Note: if anything is added here, it should be added also in the SMP block of calc_biases()
  has_data = true;
  return COLVARS_OK;
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

// Static members

std::vector<cvm::deps::feature *> colvarbias::cvb_features;
