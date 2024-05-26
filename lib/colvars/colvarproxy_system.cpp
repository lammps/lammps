// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.


#include "colvarmodule.h"
#include "colvartypes.h"
#include "colvarproxy_system.h"



colvarproxy_system::colvarproxy_system()
{
  angstrom_value_ = 0.0;
  kcal_mol_value_ = 0.0;
  target_temperature_ = 0.0;
  boltzmann_ = 0.001987191; // Default: kcal/mol/K
  boundaries_type = boundaries_unsupported;
  total_force_requested = false;
  indirect_lambda_biasing_force = 0.0;
  cached_alch_lambda_changed = false;
  cached_alch_lambda = -1.0;
  reset_pbc_lattice();
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


cvm::real colvarproxy_system::dt()
{
  // TODO define, document and implement a user method to set the value of this
  return 1.0;
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
    cvm::error("Error: total forces are currently not implemented.\n",
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


inline int round_to_integer(cvm::real x)
{
  return int(cvm::floor(x+0.5));
}


void colvarproxy_system::update_pbc_lattice()
{
  // Periodicity is assumed in all directions

  if (boundaries_type == boundaries_unsupported ||
      boundaries_type == boundaries_non_periodic) {
    cvm::error("Error: setting PBC lattice with unsupported boundaries.\n",
               COLVARS_BUG_ERROR);
    return;
  }

  {
    cvm::rvector const v = cvm::rvector::outer(unit_cell_y, unit_cell_z);
    reciprocal_cell_x = v/(v*unit_cell_x);
  }
  {
    cvm::rvector const v = cvm::rvector::outer(unit_cell_z, unit_cell_x);
    reciprocal_cell_y = v/(v*unit_cell_y);
  }
  {
    cvm::rvector const v = cvm::rvector::outer(unit_cell_x, unit_cell_y);
    reciprocal_cell_z = v/(v*unit_cell_z);
  }
}


void colvarproxy_system::reset_pbc_lattice()
{
  unit_cell_x.reset();
  unit_cell_y.reset();
  unit_cell_z.reset();
  reciprocal_cell_x.reset();
  reciprocal_cell_y.reset();
  reciprocal_cell_z.reset();
}


cvm::rvector colvarproxy_system::position_distance(cvm::atom_pos const &pos1,
                                                   cvm::atom_pos const &pos2)
  const
{
  if (boundaries_type == boundaries_unsupported) {
    cvm::error("Error: unsupported boundary conditions.\n", COLVARS_INPUT_ERROR);
  }

  cvm::rvector diff = (pos2 - pos1);

  if (boundaries_type == boundaries_non_periodic) return diff;

  cvm::real const x_shift = round_to_integer(reciprocal_cell_x*diff);
  cvm::real const y_shift = round_to_integer(reciprocal_cell_y*diff);
  cvm::real const z_shift = round_to_integer(reciprocal_cell_z*diff);

  diff.x -= x_shift*unit_cell_x.x + y_shift*unit_cell_y.x +
    z_shift*unit_cell_z.x;
  diff.y -= x_shift*unit_cell_x.y + y_shift*unit_cell_y.y +
    z_shift*unit_cell_z.y;
  diff.z -= x_shift*unit_cell_x.z + y_shift*unit_cell_y.z +
    z_shift*unit_cell_z.z;

  return diff;
}


int colvarproxy_system::get_molid(int &)
{
  cvm::error("Error: only VMD allows the use of multiple \"molecules\", "
             "i.e. multiple molecular systems.", COLVARS_NOT_IMPLEMENTED);
  return -1;
}


int colvarproxy_system::get_alch_lambda(cvm::real * /* lambda */)
{
  return cvm::error("Error in get_alch_lambda: alchemical lambda dynamics is not supported by this build.",
    COLVARS_NOT_IMPLEMENTED);
}


void colvarproxy_system::set_alch_lambda(cvm::real lambda)
{
  cached_alch_lambda = lambda;
  cached_alch_lambda_changed = true;
}


int colvarproxy_system::send_alch_lambda()
{
  return cvm::error("Error in set_alch_lambda: alchemical lambda dynamics is not supported by this build.",
    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy_system::get_dE_dlambda(cvm::real * /* force */)
{
  return cvm::error("Error in get_dE_dlambda: alchemical lambda dynamics is not supported by this build.",
    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy_system::apply_force_dE_dlambda(cvm::real* /* force */)
{
  return cvm::error("Error in apply_force_dE_dlambda: function is not implemented by this build.",
    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy_system::get_d2E_dlambda2(cvm::real*)
{
  return cvm::error("Error in get_d2E_dlambda2: function is not implemented by this build.",
    COLVARS_NOT_IMPLEMENTED);
}
