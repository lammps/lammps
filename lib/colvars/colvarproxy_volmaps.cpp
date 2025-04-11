// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarmodule.h"
#include "colvarproxy_volmaps.h"
#include "colvarmodule_utils.h"


colvarproxy_volmaps::colvarproxy_volmaps()
{
  volmaps_rms_applied_force_ = volmaps_max_applied_force_ = 0.0;
}


colvarproxy_volmaps::~colvarproxy_volmaps() {}


int colvarproxy_volmaps::check_volmaps_available()
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_volmaps::reset()
{
  for (size_t i = 0; i < volmaps_ids.size(); i++) {
    clear_volmap(i);
  }
  volmaps_ids.clear();
  volmaps_refcount.clear();
  volmaps_values.clear();
  volmaps_new_colvar_forces.clear();
  return COLVARS_OK;
}


int colvarproxy_volmaps::add_volmap_slot(int volmap_id)
{
  volmaps_ids.push_back(volmap_id);
  volmaps_refcount.push_back(1);
  volmaps_values.push_back(0.0);
  volmaps_new_colvar_forces.push_back(0.0);
  return (volmaps_ids.size() - 1);
}


int colvarproxy_volmaps::check_volmap_by_id(int /* volmap_id */)
{
  return cvm::error("Error: selecting volumetric maps is not available.\n",
                    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy_volmaps::check_volmap_by_name(const char * /* volmap_name */)
{
  return cvm::error("Error: selecting volumetric maps by name is not "
                    "available.\n", COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy_volmaps::init_volmap_by_name(char const * /* volmap_name */)
{
  return -1;
}


int colvarproxy_volmaps::init_volmap_by_id(int /* volmap_id */)
{
  return -1;
}


int colvarproxy_volmaps::init_volmap_by_name(std::string const &volmap_name)
{
  return init_volmap_by_name(volmap_name.c_str());
}


int colvarproxy_volmaps::check_volmap_by_name(std::string const &volmap_name)
{
  return check_volmap_by_name(volmap_name.c_str());
}


void colvarproxy_volmaps::clear_volmap(int index)
{
  if (((size_t) index) >= volmaps_ids.size()) {
    cvm::error("Error: trying to unrequest a volumetric map that was not "
               "previously requested.\n", COLVARS_INPUT_ERROR);
  }

  if (volmaps_refcount[index] > 0) {
    volmaps_refcount[index] -= 1;
  }
}


int colvarproxy_volmaps::get_volmap_id_from_name(char const *volmap_name)
{
  // Raise error
  colvarproxy_volmaps::check_volmap_by_name(volmap_name);
  return -1;
}


int colvarproxy_volmaps::compute_volmap(int /* flags */,
                                        int /* volmap_id */,
                                        cvm::atom_iter /* atom_begin */,
                                        cvm::atom_iter /* atom_end */,
                                        cvm::real * /* value */,
                                        cvm::real * /* atom_field */)
{
  return COLVARS_NOT_IMPLEMENTED;
}


void colvarproxy_volmaps::compute_rms_volmaps_applied_force()
{
  volmaps_rms_applied_force_ =
    compute_norm2_stats<cvm::real, 0, false>(volmaps_new_colvar_forces);
}


void colvarproxy_volmaps::compute_max_volmaps_applied_force()
{
  volmaps_max_applied_force_ =
    compute_norm2_stats<cvm::real, 1, false>(volmaps_new_colvar_forces);
}
