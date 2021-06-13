// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarmodule.h"
#include "colvarproxy_volmaps.h"


colvarproxy_volmaps::colvarproxy_volmaps() {}


colvarproxy_volmaps::~colvarproxy_volmaps() {}


int colvarproxy_volmaps::volmaps_available()
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_volmaps::reset()
{
  for (size_t i = 0; i < volmaps_ids.size(); i++) {
    clear_volmap(i);
  }
  volmaps_ids.clear();
  volmaps_ncopies.clear();
  volmaps_values.clear();
  volmaps_new_colvar_forces.clear();
  return COLVARS_OK;
}


int colvarproxy_volmaps::add_volmap_slot(int volmap_id)
{
  volmaps_ids.push_back(volmap_id);
  volmaps_ncopies.push_back(1);
  volmaps_values.push_back(0.0);
  volmaps_new_colvar_forces.push_back(0.0);
  return (volmaps_ids.size() - 1);
}


int colvarproxy_volmaps::init_volmap(int volmap_id)
{
  return cvm::error("Error: access to volumetric maps is unavailable "
                    "in this build.\n",
                    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy_volmaps::init_volmap(const char *volmap_name)
{
  return cvm::error("Error: access to volumetric maps is unavailable "
                    "in this build.\n",
                    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy_volmaps::init_volmap(const std::string &volmap_name)
{
  return init_volmap(volmap_name.c_str());
}


void colvarproxy_volmaps::clear_volmap(int index)
{
  if (((size_t) index) >= volmaps_ids.size()) {
    cvm::error("Error: trying to unrequest a volumetric map that was not "
               "previously requested.\n", INPUT_ERROR);
  }

  if (volmaps_ncopies[index] > 0) {
    volmaps_ncopies[index] -= 1;
  }
}
