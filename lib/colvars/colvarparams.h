// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARPARAMS_H
#define COLVARPARAMS_H

#include <map>

/// \file colvarparams.h Functions to handle scalar parameters used in objects


class colvarparams {

 public:

  /// Whether the parameter param_name exists
  int param_exists(std::string const &param_name);

  /// Get a copy of the names of registered parameters
  virtual std::vector<std::string> get_param_names();

  /// Get a copy of the names of registered parameter gradients
  virtual std::vector<std::string> get_param_grad_names();

  /// Pointer to the parameter param_name
  virtual void const *get_param_ptr(std::string const &param_name);

  /// Pointer to the gradient of parameter param_name
  virtual void const *get_param_grad_ptr(std::string const &param_name);

  /// Value of the parameter param_name (must be a scalar)
  virtual cvm::real get_param(std::string const &param_name);

  /// Set the named parameter to the given value
  virtual int set_param(std::string const &param_name, void const *new_value);

 protected:

  /// Default constructor
  colvarparams();

  /// Default destructor
  virtual ~colvarparams();

  /// Pointers to relevant parameters that may be accessed by other objects
  std::map<std::string, void const *> param_map;

  /// Derivatives of the object with respect to internal parameters
  std::map<std::string, colvarvalue const *> param_grad_map;

  /// Register the given parameter
  void register_param(std::string const &param_name, void *param_ptr);

  /// Register the gradient of the given parameter
  void register_param_grad(std::string const &param_name,
                           colvarvalue *param_grad_ptr);

};

#endif
