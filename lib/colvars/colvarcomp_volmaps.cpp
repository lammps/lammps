// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvar.h"
#include "colvarcomp.h"



colvar::map_total::map_total()
  : cvc(), volmap_index(-1)
{
  function_type = "map_total";
  x.type(colvarvalue::type_scalar);
}


colvar::map_total::map_total(std::string const &conf)
  : cvc(), volmap_index(-1)
{
  function_type = "map_total";
  x.type(colvarvalue::type_scalar);
  map_total::init(conf);
}


int colvar::map_total::init(std::string const &conf)
{
  int error_code = cvc::init(conf);
  get_keyval(conf, "mapName", map_name, map_name);
  volmap_index = (cvm::proxy)->init_volmap(map_name);
  error_code |= volmap_index > 0 ? COLVARS_OK : INPUT_ERROR;
  return error_code;
}


void colvar::map_total::calc_value()
{
  x.real_value = (cvm::proxy)->get_volmap_value(volmap_index);
}


void colvar::map_total::calc_gradients()
{
  // Atomic coordinates are not available here
}


void colvar::map_total::apply_force(colvarvalue const &force)
{
  (cvm::proxy)->apply_volmap_force(volmap_index, force.real_value);
}
