// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <algorithm>

#include "colvarmodule.h"
#include "colvartypes.h"
#include "colvarproxy_system.h"



colvarproxy_system::colvarproxy_system()
{
  angstrom_value_ = 0.0;
  kcal_mol_value_ = 0.0;
  timestep_ = 1.0;
  target_temperature_ = 0.0;
  boltzmann_ = 0.001987191; // Default: kcal/mol/K
  total_force_requested = false;
  indirect_lambda_biasing_force = 0.0;
  cached_alch_lambda_changed = false;
  cached_alch_lambda = -1.0;
}


colvarproxy_system::~colvarproxy_system() {}


int colvarproxy_system::set_unit_system(std::string const & /* units */,
                                        bool /* check_only */)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_system::set_target_temperature(cvm::real T)
{
  target_temperature_ = T;
  return COLVARS_OK;
}


int colvarproxy_system::set_integration_timestep(cvm::real dt)
{
  timestep_ = dt;
  return COLVARS_OK;
}

int colvarproxy_system::set_time_step_factor(int fact)
{
  time_step_factor_ = fact;
  return COLVARS_OK;
}

cvm::real colvarproxy_system::rand_gaussian()
{
  // TODO define, document and implement a user method to set the value of this
  return 0.0;
}


void colvarproxy_system::add_energy(cvm::real /* energy */) {}


void colvarproxy_system::request_total_force(bool yesno)
{
  if (yesno == true)
    cvm::error_static("Error: total forces are currently not implemented.\n",
               COLVARS_NOT_IMPLEMENTED);
}


bool colvarproxy_system::total_forces_enabled() const
{
  return false;
}


bool colvarproxy_system::total_forces_same_step() const
{
  return false;
}


cvm::rvector colvarproxy_system::position_distance(cvm::atom_pos const &pos1,
                                                   cvm::atom_pos const &pos2) const
{
  return boundaries_.position_distance(pos1, pos2);
}


int colvarproxy_system::get_molid(int &)
{
  cvm::error_static("Error: only VMD allows the use of multiple \"molecules\", "
             "i.e. multiple molecular systems.", COLVARS_NOT_IMPLEMENTED);
  return -1;
}


int colvarproxy_system::get_alch_lambda(cvm::real * /* lambda */)
{
  return cvm::error_static("Error in get_alch_lambda: alchemical lambda dynamics is not supported by this build.",
    COLVARS_NOT_IMPLEMENTED);
}


void colvarproxy_system::set_alch_lambda(cvm::real lambda)
{
  cached_alch_lambda = lambda;
  cached_alch_lambda_changed = true;
}


int colvarproxy_system::send_alch_lambda()
{
  return cvm::error_static("Error in set_alch_lambda: alchemical lambda dynamics is not supported by this build.",
    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy_system::get_dE_dlambda(cvm::real * /* force */)
{
  return cvm::error_static("Error in get_dE_dlambda: alchemical lambda dynamics is not supported by this build.",
    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy_system::apply_force_dE_dlambda(cvm::real* /* force */)
{
  return cvm::error_static("Error in apply_force_dE_dlambda: function is not implemented by this build.",
    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy_system::get_d2E_dlambda2(cvm::real*)
{
  return cvm::error_static("Error in get_d2E_dlambda2: function is not implemented by this build.",
    COLVARS_NOT_IMPLEMENTED);
}
