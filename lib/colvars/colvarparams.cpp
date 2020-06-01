// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <sstream>
#include <iostream>
#include <algorithm>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparams.h"



colvarparams::colvarparams()
{}


colvarparams::~colvarparams()
{}


void colvarparams::register_param(std::string const &param_name,
                                  void *param_ptr)
{
  param_map[param_name] = param_ptr;
}


void colvarparams::register_param_grad(std::string const &param_name,
                                       colvarvalue *param_grad_ptr)
{
  param_grad_map[param_name] = param_grad_ptr;
}


int colvarparams::param_exists(std::string const &param_name)
{
  if (param_map.count(param_name) > 0) {
    return COLVARS_OK;
  }
  return INPUT_ERROR;
}


std::vector<std::string> colvarparams::get_param_names()
{
  std::vector<std::string> result;
  for (std::map<std::string, void const *>::const_iterator elem =
         param_map.begin(); elem != param_map.end(); elem++) {
    result.push_back(elem->first);
  }
  return result;
}


std::vector<std::string> colvarparams::get_param_grad_names()
{
  std::vector<std::string> result;
  for (std::map<std::string, colvarvalue const *>::const_iterator elem =
         param_grad_map.begin(); elem != param_grad_map.end(); elem++) {
    result.push_back(elem->first);
  }
  return result;
}


void const *colvarparams::get_param_ptr(std::string const &param_name)
{
  if (param_map.count(param_name) > 0) {
    return param_map[param_name];
  }
  cvm::error("Error: parameter \""+param_name+"\" not found.\n", INPUT_ERROR);
  return NULL;
}


void const *colvarparams::get_param_grad_ptr(std::string const &param_name)
{
  if (param_grad_map.count(param_name) > 0) {
    return param_grad_map[param_name];
  }
  cvm::error("Error: gradient of parameter \""+param_name+"\" not found.\n",
             INPUT_ERROR);
  return NULL;
}


cvm::real colvarparams::get_param(std::string const &param_name)
{
  cvm::real const *ptr =
    reinterpret_cast<cvm::real const *>(get_param_ptr(param_name));
  return ptr != NULL ? *ptr : 0.0;
}


int colvarparams::set_param(std::string const &param_name,
                            void const *new_value)
{
  if (param_map.count(param_name) > 0) {
    return cvm::error("Error: parameter \""+param_name+"\" cannot be "
                      "modified.\n", COLVARS_NOT_IMPLEMENTED);
  }
  return cvm::error("Error: parameter \""+param_name+"\" not found.\n",
                    INPUT_ERROR);
}
